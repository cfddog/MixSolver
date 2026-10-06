!===============================================================================
! main.f90 -- mixed structured/unstructured weak-coupling driver (phase 6).
!
! Architecture:
!   - All ranks share MPI_COMM_WORLD.
!   - Ranks are split into two sub-communicators by colour:
!       STRUCT_GROUP (colour 0) -> runs the structured solver
!       UNS_GROUP    (colour 1) -> runs the unstructured solver
!   - The split is controlled by the first n_struct ranks (from mix.control or
!     command line).  Default: rank 0 = struct, everyone else = uns.
!   - Each group runs its own solver init + coupling iteration loop.  Between
!     iterations the interface state is exchanged across the groups (phase 5
!     units + exchange modules).
!
! Phase-6 scope:
!   - Unstructured side: fully wired (init / set-interface-BC / step / extract).
!   - Structured side: fully wired (init / step / extract / set-ghost-cell).
!   - MPI exchange: rank 0 of each group exchanges interface SI state via
!     blocking Send/Recv across COMM_WORLD (mod_coupling_exchange).
!
! Usage: mpirun -np N bin/mixsolver_mpi  [mix.control]  [struct.cas/.control]
!                                                     [uns.cas uns.control]
!===============================================================================
program mixsolver
   use mpi
   use mod_precision, only: dp
   use mod_reference_state, only: read_mix_control
   implicit none

   integer, parameter :: STRUCT_GROUP = 0, UNS_GROUP = 1

   integer :: ierr, rank, nproc, colour, key
   integer :: comm_struct, comm_uns
   integer :: n_struct_ranks
   character(len=512) :: mixfile, scas, sctl, ucas, uctl

   ! ---- MPI init -------------------------------------------------------------
   call MPI_Init( ierr )
   call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
   call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )

   ! ---- arguments ------------------------------------------------------------
   mixfile = 'mix.control'
   scas    = 'grid_BC/Mesh3d.x'
   sctl    = 'grid_BC/control.ec'
   ucas    = 'grid_BC/unMesh.cas'
   uctl    = 'grid_BC/unMesh.control'
   if ( command_argument_count() >= 1 ) call get_command_argument(1, mixfile)
   if ( command_argument_count() >= 2 ) call get_command_argument(2, scas)
   if ( command_argument_count() >= 3 ) call get_command_argument(3, sctl)
   if ( command_argument_count() >= 4 ) call get_command_argument(4, ucas)
   if ( command_argument_count() >= 5 ) call get_command_argument(5, uctl)

   ! ---- reference state (all ranks read the same mix.control) ----------------
   call read_mix_control( trim(mixfile), ierr )

   ! ---- split communicator ---------------------------------------------------
   ! Default: rank 0 -> structured, rest -> unstructured.
   ! Special case: nproc==1 -> single rank runs both (struct first, then uns).
   n_struct_ranks = 1
   if ( nproc == 1 ) then
      n_struct_ranks = 0   ! all ranks -> uns (no struct side)
   end if
   if ( rank < n_struct_ranks ) then
      colour = STRUCT_GROUP
   else
      colour = UNS_GROUP
   end if
   key = rank

   call MPI_Comm_split( MPI_COMM_WORLD, colour, key, comm_struct, ierr )
   comm_uns = comm_struct

   if ( rank == 0 ) then
      write(*,'(a)') ''
      write(*,'(a)') '=========================================================='
      write(*,'(a)') '  MixNSSolver -- weak-coupling driver (phase 6)'
      write(*,'(a)') '=========================================================='
      write(*,'(a,i0)') '  total ranks   : ', nproc
      write(*,'(a,i0)') '  struct ranks  : ', n_struct_ranks
      write(*,'(a,i0)') '  uns ranks     : ', nproc - n_struct_ranks
      write(*,'(a,a)')    '  mix.control   : ', trim(mixfile)
      write(*,'(a)') '=========================================================='
   end if

   ! ---- dispatch by group ----------------------------------------------------
   select case ( colour )
   case ( STRUCT_GROUP )
      call struct_group_driver( comm_struct, scas, sctl, rank, nproc )
   case ( UNS_GROUP )
      call uns_group_driver( comm_uns, ucas, uctl, rank, nproc )
   end select

   call MPI_Finalize( ierr )

contains

   !===========================================================================
   ! Structured-group driver.
   ! Initialises the structured solver, then runs the coupling loop:
   !   1. Extract interface state (non-dim), convert to SI, send to uns root
   !   2. Receive uns SI state from uns root, convert to non-dim, set ghost cells
   !   3. Advance structured solver by one time step
   !===========================================================================
   subroutine struct_group_driver( comm, casfile, ctlfile, rank, nproc )
      use mod_struct_driver, only: struct_solver_init, struct_solver_step, &
                                    struct_extract_iface, struct_set_iface_bc, &
                                    struct_solver_save
      use mod_coupling_exchange, only: exchange_struct_to_uns, recv_uns_iface_state
      use mod_interface_units, only: struct_to_SI, SI_to_struct
      use mod_reference_state, only: get_ref_state, reference_state_t, &
                                     get_coupling_params
      implicit none

      integer,          intent(in) :: comm, rank, nproc
      character(len=*), intent(in) :: casfile, ctlfile

      integer :: iter, iter0, n_couple, n_uns_root, ierr2
      integer :: nfaces, n_uns_recv, n_uns_steps_d, iface_ramp_d
      integer :: n_struct_steps_d, isub, save_interval, couple_restart
      real(dp) :: iface_relax_d
      real(dp), allocatable :: rho_nd(:), u_nd(:), v_nd(:), w_nd(:), T_nd(:), p_nd(:)
      real(dp), allocatable :: urho(:), uu(:,:), uT(:), up(:)
      real(dp), allocatable :: rho_nd2(:), u_nd2(:), v_nd2(:), w_nd2(:), T_nd2(:), p_nd2(:)
      real(dp), allocatable :: rho_si(:), u_si(:), v_si(:), w_si(:), T_si(:), p_si(:)
      type(reference_state_t) :: r

      call get_coupling_params( n_couple, n_uns_steps_d, iface_relax_d, iface_ramp_d, &
                                n_struct_steps = n_struct_steps_d, &
                                save_interval = save_interval, &
                                couple_restart = couple_restart )
      n_uns_root = 1   ! uns root is global rank 1 (struct is rank 0)
      write(*,'(a,i0)') '[struct rank 0] struct pseudo-time substeps/coupling iter: ', &
                        n_struct_steps_d

      ! joint restart: recover the last saved coupling iteration count
      iter0 = 0
      if ( couple_restart == 1 ) call read_couple_state( iter0 )

      write(*,'(a,i0,a)') '[struct rank ', rank, '] initialising structured solver'
      call struct_solver_init(comm, ctlfile, force_restart = (iter0 > 0) )
      write(*,'(a,i0,a)') '[struct rank ', rank, '] init OK'
      if ( iter0 > 0 ) &
         write(*,'(a,i0,a,i0)') '[struct rank ', rank, '] restart from coupling iter ', &
                                iter0

      ! ---- coupling loop ------------------------------------------------------
      do iter = iter0+1, n_couple
         ! extract interface state (non-dim)
         call struct_extract_iface(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, nfaces)
         if (iter == 1) &
            write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, '] interface faces: ', nfaces

         ! send struct SI state to uns root (count + payload, works for nfaces=0)
         if (nfaces > 0) then
            call exchange_struct_to_uns(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                        nfaces, n_uns_root, ierr2)
         else
            call exchange_struct_to_uns(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                        0, n_uns_root, ierr2)
         end if

         ! receive uns SI state from uns root
         call recv_uns_iface_state(n_uns_recv, urho, uu, uT, up, n_uns_root, ierr2)
         if (iter == 1) &
            write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, '] received uns faces: ', n_uns_recv

         ! periodic SI diagnostics (interface means, both peers)
         if ( mod(iter,25)==0 .or. iter==1 .or. iter==n_couple ) then
            if (nfaces > 0) then
               allocate(rho_si(nfaces), u_si(nfaces), v_si(nfaces), &
                        w_si(nfaces), T_si(nfaces), p_si(nfaces))
               call convert_to_si(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                  rho_si, u_si, v_si, w_si, T_si, p_si, nfaces)
               write(*,'(a,i0,a,i0,a,3es12.4)') &
                  '[struct rank ', rank, '] diag iter ', iter, &
                  ' struct SI: mean u_x, mean p(abs), mean T =', &
                  sum(u_si)/nfaces, sum(p_si)/nfaces, sum(T_si)/nfaces
               deallocate(rho_si, u_si, v_si, w_si, T_si, p_si)
            end if
            if (n_uns_recv > 0) then
               write(*,'(a,i0,a,i0,a,3es12.4)') &
                  '[struct rank ', rank, '] diag iter ', iter, &
                  ' uns    SI: mean u_x, mean p(gauge), mean T =', &
                  sum(uu(1,:))/n_uns_recv, sum(up)/n_uns_recv, &
                  sum(uT)/n_uns_recv
            end if
         end if

         ! advance structured solver (subcycled pseudo-time steps to keep
         ! pace with the uns side, which takes several SIMPLE outer steps)
         do isub = 1, n_struct_steps_d
            call struct_solver_step()
         end do

         ! set ghost cells from the received uns first-cell state.  Done
         ! after the substeps.  The characteristic back-pressure strength is
         ! alpha = ALPHA_MAX * min(1, iter/iface_ramp): the ramp cushions the
         ! cold start, and ALPHA_MAX<1 permanently damps the partitioned
         ! coupling loop (full alpha=1 was empirically unstable at Ma=0.1:
         ! the acoustic gain 1/(rho*c) makes the exchange loop gain > 1).
         ! At convergence p1 -> p_peer for any alpha>0 because the struct
         ! inner pressure relaxes toward pb every coupling iteration.
         if (nfaces > 0 .and. n_uns_recv > 0) then
            block
               real(dp), parameter :: ALPHA_MAX = 0.3_dp
               real(dp) :: alpha_p
            allocate(rho_nd2(nfaces), u_nd2(nfaces), v_nd2(nfaces), &
                     w_nd2(nfaces), T_nd2(nfaces), p_nd2(nfaces))
            r = get_ref_state()
            call convert_uns_to_struct_nd(urho, uu, uT, up, n_uns_recv, &
                                          rho_nd2, u_nd2, v_nd2, w_nd2, T_nd2, p_nd2, &
                                          nfaces, r)
            alpha_p = ALPHA_MAX * min(1.0_dp, &
                        real(iter,dp)/real(max(iface_ramp_d,1),dp))
            call struct_set_iface_bc(rho_nd2, u_nd2, v_nd2, w_nd2, T_nd2, p_nd2, nfaces, &
                                     alpha = alpha_p)
            deallocate(rho_nd2, u_nd2, v_nd2, w_nd2, T_nd2, p_nd2)
            end block
         else if (nfaces > 0) then
            ! no uns data received; keep current state as BC (extrapolation)
            call struct_set_iface_bc(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, nfaces, &
                                     alpha = 0.0_dp)
         end if
         if ( mod(iter,100)==0 ) &
            write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, '] coupling iter ', iter, ' done'

         ! phase-10 auto-save: aligned with the uns side at the same coupling
         ! iter (after both sides have exchanged and stepped).
         if ( save_interval > 0 ) then
            if ( mod(iter, save_interval) == 0 ) then
               call struct_solver_save
               if ( rank == 0 ) call write_couple_state( iter )
               write(*,'(a,i0,a,i0,a)') '[struct rank ', rank, &
                     '] saved flow3d.dat at coupling iter ', iter
            end if
         end if

         if (allocated(rho_nd)) deallocate(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd)
         if (allocated(urho)) deallocate(urho, uu, uT, up)
      end do

      write(*,'(a,i0,a)') '[struct rank ', rank, '] coupling loop done'
   end subroutine struct_group_driver

   !===========================================================================
   ! Unstructured-group driver.
   ! Initialises the unstructured solver, then runs the coupling loop:
   !   1. Receive struct SI state from struct root, set interface BC
   !   2. Advance SIMPLE by n_uns_steps
   !   3. Extract interface SI state, send to struct root
   !===========================================================================
   subroutine uns_group_driver( comm, casfile, ctlfile, rank, nproc )
      use mod_uns_driver, only: uns_solver_init, uns_solver_step, &
                                uns_solver_extract_iface
      use mod_uns_mesh, only: mesh_t
      use mod_uns_connectivity, only: conn_t
      use mod_uns_geometry, only: geom_t
      use mod_uns_control, only: ctrl_t, BC_INTERFACE, BC_POUTLET, BC_FARFIELD
      use mod_uns_bc, only: bc_t, set_interface_vel, &
                            set_interface_T, bc_face_vel
      use mod_uns_fields, only: fields_t
#ifdef HAVE_MPI
      ! The uns field dump/restart lives in the MPI-only mod_uns_restart stack
      ! (-> mod_uns_mpi_core / mod_uns_partition / mod_uns_local_mesh).  The
      ! serial 'single-process MPI' build (make all) omits that stack, so the
      ! save path below is compiled out there; the MPI build enables it.
      use mod_uns_restart, only: write_field_dump
#endif
      use mod_coupling_exchange, only: exchange_uns_to_struct
      use mod_reference_state, only: get_coupling_params
      implicit none

      integer,          intent(in) :: comm, rank, nproc
      character(len=*), intent(in) :: casfile, ctlfile

      type(mesh_t)   :: m
      type(conn_t)   :: c
      type(geom_t)   :: g
      type(ctrl_t)   :: ctrl
      type(bc_t)     :: bcs
      type(fields_t) :: fld
      integer        :: ier, iter, iter0, n_couple, n_uns_steps, n_struct_root
      integer, allocatable :: ifaces(:)
      real(dp), allocatable :: irho(:), iT(:), ip(:), iu(:,:)
      real(dp), allocatable :: srho(:), su(:,:), sT(:), sp(:)
      real(dp), allocatable :: su_bc(:,:)
      integer :: n_uns_faces, iface_ramp, save_interval, couple_restart
      real(dp) :: iface_relax, omega
      logical :: has_struct
      character(len=*), parameter :: uns_dump = 'unMesh_restart.dat'

      call get_coupling_params( n_couple, n_uns_steps, iface_relax, iface_ramp, &
                                save_interval = save_interval, &
                                couple_restart = couple_restart )
      n_struct_root = 0   ! struct root is global rank 0
      ! When there is no struct group (nproc==1, n_struct_ranks==0) the
      ! exchange is skipped entirely.
      has_struct = ( nproc > 1 )

      ! joint restart: recover the last saved coupling iteration count
      iter0 = 0
      if ( couple_restart == 1 ) call read_couple_state( iter0 )

      write(*,'(a,i0,a)') '[uns rank ', rank, '] initialising unstructured solver'
      if ( iter0 > 0 ) then
         call uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier, &
                               restart_file = uns_dump )
      else
         call uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier )
      end if
      if ( ier /= 0 ) then
         write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] init FAILED ier=', ier
         call MPI_Abort( MPI_COMM_WORLD, ier, ierr )
         return
      end if
      write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] init OK, ncells=', m%ncells
      if ( iter0 > 0 ) &
         write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] restart from coupling iter ', iter0

      ! ---- extract initial interface state (for exchange sizing) ---------------
      call uns_solver_extract_iface( m, g, bcs, fld, ifaces, irho, iu, iT, ip, ier, ctrl%rho )
      n_uns_faces = size(ifaces)
      write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] interface faces: ', n_uns_faces

      ! ---- coupling loop ------------------------------------------------------
      do iter = iter0+1, n_couple

         ! receive struct SI state, send uns SI state back
         ! (skipped when there is no struct group, e.g. nproc==1)
         if ( has_struct ) then
            call exchange_uns_to_struct(irho, iu, iT, ip, n_uns_faces, &
                                        n_struct_root, srho, su, sT, sp, ierr)

            ! set interface BC from struct state with under-relaxation and
            ! linear ramp over the first iface_ramp coupling iterations:
            !   omega  = iface_relax * min(1, iter/iface_ramp)
            !   BC     = omega*struct + (1-omega)*uns_current
            ! (iu/ip hold the uns interface state from the previous extract)
            if (n_uns_faces > 0) then
               omega = iface_relax * min( 1.0_dp, &
                            real(iter,dp) / real(max(iface_ramp,1),dp) )
               allocate( su_bc(3,n_uns_faces) )
               su_bc = omega * su + (1.0_dp - omega) * iu
               ! Interface partition (phase-11 flow B, Dirichlet-Neumann):
               ! the uns side imposes only the peer face VELOCITY here;
               ! interface pressure is zero-gradient (bc_face_p for
               ! BC_INTERFACE returns pP), while the struct side imposes the
               ! uns back pressure via its subsonic-outflow characteristic
               ! ghost state.  Setting pressure as well over-constrains both
               ! sides and freezes a kPa-scale pressure jump.  The received
               ! sp (absolute Pa) is retained for diagnostics only.
               ! mass-conservation fix: with all-Dirichlet velocity
               ! boundaries (mass-flow-inlet + interface) the net boundary
               ! volume flux must vanish or the pure-Neumann PPE has no
               ! solution.  Sum the flux over all non-interface boundaries
               ! with the CURRENT bc state (includes the inlet ramp factor),
               ! then add a uniform normal correction to the interface face
               ! velocities so the total boundary flux is zero.
               block
                  integer  :: kf, gi, fi
                  logical  :: has_p_bc
                  real(dp) :: uf(3), F_other, F_iface, A_iface, du_n, nvec(3)
                  ! Uniform-normal flux correction is required only when every
                  ! non-interface boundary is velocity-Dirichlet (pure-Neumann
                  ! PPE).  With a pressure-specified boundary (pressure-outlet
                  ! / far-field) the PPE is well posed and the interface inflow
                  ! must NOT be cancelled to balance the not-yet-developed
                  ! outlet flux: that blocks the inlet at startup, the struct
                  ! side compresses against a closed wall, and the interface
                  ! pressure gradient explodes.
                  has_p_bc = any( bcs%gb(:)%btype == BC_POUTLET ) .or. &
                             any( bcs%gb(:)%btype == BC_FARFIELD )
                  F_other = 0.0_dp
                  do kf = 1, m%nfaces
                     if ( m%f(kf)%c1 /= 0 ) cycle
                     gi = bcs%fgrp(kf)
                     if ( gi == 0 ) cycle
                     if ( bcs%gb(gi)%btype == BC_INTERFACE ) cycle
                     call bc_face_vel( bcs, kf, g%xf(:,kf), g%sf(:,kf), &
                                       fld%u(:,m%f(kf)%c0), uf )
                     F_other = F_other + dot_product( uf, g%sf(:,kf) )
                  end do
                  F_iface = 0.0_dp
                  A_iface = 0.0_dp
                  do fi = 1, n_uns_faces
                     kf = ifaces(fi)
                     F_iface = F_iface + dot_product( su_bc(:,fi), g%sf(:,kf) )
                     A_iface = A_iface + g%area(kf)
                  end do
                  du_n = 0.0_dp
                  if ( .not. has_p_bc .and. A_iface > 0.0_dp ) then
                     du_n = -( F_other + F_iface ) / A_iface
                     do fi = 1, n_uns_faces
                        kf = ifaces(fi)
                        nvec = g%sf(:,kf) / g%area(kf)
                        su_bc(:,fi) = su_bc(:,fi) + du_n * nvec
                     end do
                  end if
                  if ( iter <= 8 .or. mod(iter,25)==0 ) then
                     write(*,'(a,i0,a,l1,a,3es11.3)') '[uns rank ', rank, &
                           '] iface has_p_bc=', has_p_bc, &
                           ' F_other/F_iface/du_n =', F_other, F_iface, du_n
                     write(*,'(a,i0,a,6es11.3)') '[uns rank ', rank, &
                           '] recv su x min/max, sT min/max, sp min/max =', &
                           minval(su(1,:)), maxval(su(1,:)), &
                           minval(sT), maxval(sT), minval(sp), maxval(sp)
                     write(*,'(a,i0,a,4es11.3)') '[uns rank ', rank, &
                           '] BC: su_bc x min/max =', &
                           minval(su_bc(1,:)), maxval(su_bc(1,:))
                  end if
               end block
               if ( mod(iter,100)==0 .or. iter==1 ) &
                  write(*,'(a,i0,a,i0,a,f6.3)') '[uns rank ', rank, &
                        '] iter ', iter, ' iface omega = ', omega
               call set_interface_vel( bcs, ifaces, su_bc )
               call set_interface_T( bcs, ifaces, sT )
               deallocate( su_bc )
            end if
         end if

         ! advance unstructured solver
         call uns_solver_step( m, c, g, ctrl, bcs, fld, n_uns_steps, ier )
         if ( ier /= 0 ) then
            write(*,'(a,i0,a,i0)') '[uns rank ', rank, '] step FAILED ier=', ier
            call MPI_Abort( MPI_COMM_WORLD, ier, ierr )
            return
         end if

         ! extract interface state for next exchange
         call uns_solver_extract_iface( m, g, bcs, fld, ifaces, irho, iu, iT, ip, ier, ctrl%rho )
         if ( size(ifaces) > 0 .and. (mod(iter,25)==0 .or. iter==1 .or. iter==n_couple) ) then
            write(*,'(a,i0,a,i0,a,3es12.4)') '[uns rank ', rank, '] diag iter ', iter, &
                  ' uns SI: mean u_x, mean p(gauge), |u|max =', &
                  sum(iu(1,:))/size(ifaces), sum(ip)/size(ifaces), &
                  maxval( sqrt( sum(iu*iu, dim=1) ) )
         end if

         ! No reflection: pass the raw uns cell-centre state back to the
         ! struct solver.  The struct side will compute a self-consistent
         ! ghost cell (ghost = 2*desired_face - inner) so that its own
         ! extract_iface arithmetic-mean returns exactly the received value.

         ! phase-10 auto-save: write the unstructured dump aligned with the
         ! struct side.  Only the uns root writes (global rank n_struct_ranks:
         ! 1 for a coupled run, 0 for a standalone nproc==1 run); the dump is
         ! a global field, serial-read on restart.
         if ( save_interval > 0 ) then
            if ( mod(iter, save_interval) == 0 ) then
               if ( rank == merge(0, 1, nproc == 1) ) then
#ifdef HAVE_MPI
                  call write_field_dump( uns_dump, m, fld, trim(casfile), ier )
                  if ( ier == 0 ) write(*,'(a,a,a,i0)') '[uns] saved dump: ', &
                        trim(uns_dump), ' at coupling iter ', iter
#endif
                  ! standalone run (no struct group): the uns root owns
                  ! couple_state.dat so a restart can recover the iter count.
                  if ( .not. has_struct ) call write_couple_state( iter )
               end if
            end if
         end if

         if (allocated(srho)) deallocate(srho, su, sT, sp)
      end do

      write(*,'(a,i0,a)') '[uns rank ', rank, '] coupling loop done'

      ! final field snapshot for archiving (ASCII VTU) + dump
      if ( .not. has_struct ) then
         write(*,'(a)') '[uns] standalone run: skipping coupled VTU snapshot'
      else
         block
            use mod_uns_output, only: vtk_write
            integer :: iov
            call vtk_write( m, c, g, fld, 'unMesh_coupled.vtu', iov )
            if ( iov == 0 ) write(*,'(a)') '[uns] wrote unMesh_coupled.vtu'
         end block
      end if
   end subroutine uns_group_driver

   !---------------------------------------------------------------------------
   ! Joint restart state helpers (phase 10).  couple_state.dat stores the last
   ! saved coupling iteration so that both sides restart from the same iter.
   ! Only rank 0 (struct root) reads/writes; the other side reads the file
   ! independently after a barrier-free filesystem sync.
   !---------------------------------------------------------------------------
   subroutine write_couple_state( iter )
      integer, intent(in) :: iter
      integer :: u, ios
      open( newunit=u, file='couple_state.dat', status='replace', action='write', &
            iostat=ios )
      if ( ios == 0 ) then
         write(u,*) iter
         close(u)
      end if
   end subroutine write_couple_state

   subroutine read_couple_state( iter )
      integer, intent(out) :: iter
      logical :: ex
      integer :: u, ios
      iter = 0
      inquire( file='couple_state.dat', exist=ex )
      if ( .not. ex ) return
      open( newunit=u, file='couple_state.dat', status='old', action='read', &
            iostat=ios )
      if ( ios == 0 ) then
         read(u,*,iostat=ios) iter
         close(u)
      end if
   end subroutine read_couple_state

   !===========================================================================
   ! convert_to_si -- array wrapper for struct_to_SI (non-dim -> SI).
   !===========================================================================
   subroutine convert_to_si(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                            rho, u, v, w, T, p, n)
      use mod_interface_units, only: struct_to_SI
      use mod_reference_state, only: get_ref_state, reference_state_t
      implicit none
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: rho_nd(n), u_nd(n), v_nd(n), w_nd(n), T_nd(n), p_nd(n)
      real(dp), intent(out) :: rho(n), u(n), v(n), w(n), T(n), p(n)
      integer :: i
      type(reference_state_t) :: r
      r = get_ref_state()
      do i = 1, n
         call struct_to_SI(rho_nd(i), u_nd(i), v_nd(i), w_nd(i), T_nd(i), p_nd(i), &
                           rho(i), u(i), v(i), w(i), T(i), p(i), r)
      end do
   end subroutine convert_to_si

   !===========================================================================
   ! convert_uns_to_struct_nd -- map the uns SI interface state onto the
   ! struct interface faces, then SI -> struct non-dim.
   !
   ! Two mappings are supported:
   !   * n_struct == n_uns: direct 1:1 mapping in face-list order.  Both
   !     patches are generated as matching (j outer, k inner) quad arrays,
   !     so no geometric interpolation is needed in this case.
   !   * otherwise: uniform area-weighted average as a placeholder until
   !     the full uns->struct area-weighted interpolation
   !     (mod_interface_exchange) is wired into the driver.
   ! NB: the uns pressure is GAUGE; the struct non-dim pressure is absolute
   ! (p_inf* = 1/(gamma*Ma^2)), so p_ref is added before conversion.
   !===========================================================================
   subroutine convert_uns_to_struct_nd(urho, uu, uT, up, n_uns, &
                                       rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, &
                                       n_struct, r)
      use mod_interface_units, only: SI_to_struct
      use mod_reference_state, only: reference_state_t
      implicit none
      integer,  intent(in)  :: n_uns, n_struct
      real(dp), intent(in)  :: urho(n_uns), uu(3,n_uns), uT(n_uns), up(n_uns)
      real(dp), intent(out) :: rho_nd(n_struct), u_nd(n_struct), v_nd(n_struct), &
                               w_nd(n_struct), T_nd(n_struct), p_nd(n_struct)
      type(reference_state_t), intent(in) :: r

      integer :: i, j
      real(dp) :: avg_rho, avg_u(3), avg_T, avg_p
      real(dp) :: rho_s, u_s, v_s, w_s, T_s, p_s, p_abs

      if ( n_struct == n_uns ) then
         ! direct 1:1 mapping (matched patches)
         do i = 1, n_struct
            p_abs = up(i) + r%p_ref
            call SI_to_struct(urho(i), uu(1,i), uu(2,i), uu(3,i), uT(i), p_abs, &
                              rho_s, u_s, v_s, w_s, T_s, p_s, r)
            rho_nd(i) = rho_s
            u_nd(i)   = u_s
            v_nd(i)   = v_s
            w_nd(i)   = w_s
            T_nd(i)   = T_s
            p_nd(i)   = p_s
         end do
      else
         ! uniform area-weighted average (equal weights for now)
         avg_rho = sum(urho) / real(n_uns, dp)
         avg_u   = sum(uu, dim=2) / real(n_uns, dp)
         avg_T   = sum(uT) / real(n_uns, dp)
         avg_p   = sum(up) / real(n_uns, dp) + r%p_ref

         do i = 1, n_struct
            call SI_to_struct(avg_rho, avg_u(1), avg_u(2), avg_u(3), avg_T, avg_p, &
                              rho_s, u_s, v_s, w_s, T_s, p_s, r)
            rho_nd(i) = rho_s
            u_nd(i)   = u_s
            v_nd(i)   = v_s
            w_nd(i)   = w_s
            T_nd(i)   = T_s
            p_nd(i)   = p_s
         end do
      end if
   end subroutine convert_uns_to_struct_nd

end program mixsolver
