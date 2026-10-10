!===============================================================================
! mod_struct_driver.f90 -- structured-solver coupling driver (phase 6)
!
! Wraps the structured solver (OpenCFD-EC) for the mixed coupling driver:
!   - init:  read_parameter + init + set_control_para + check_mesh_quality
!            + Init_flow  (mirrors main.f90 pre-loop)
!   - step:  advance the finest mesh by one time step (NS_Time_advance)
!   - extract: read interface inner-cell state (d, u, v, w, T, p) from B%U
!   - set:   write ghost-cell state B%U(:, ig,jg,kg) from SI exchange data
!
! All variables in B%U are non-dimensional (OpenCFD-EC convention):
!   U(1)=rho*, U(2:4)=rho*u* (momentum), U(5)=rho*E* (total energy density)
! The ghost-cell state is written as primitive variables converted to
! conservative form here (rho, rho*u, rho*v, rho*w, rho*E).
!
! Unit convention: struct solver is non-dimensional; the exchange layer
! (mod_interface_units) converts SI <-> struct non-dim before/after
! calling these routines.
!===============================================================================
module mod_struct_driver
   use precision_EC
   use mod_precision, only: dp
   implicit none
   private

   public :: struct_solver_init
   public :: struct_solver_step
   public :: struct_extract_iface
   public :: struct_set_iface_bc
   public :: struct_solver_save
   public :: struct_solver_run

contains

   !---------------------------------------------------------------------------
   ! struct_solver_init -- full structured-solver initialisation.
   ! Mirrors the pre-loop part of src/structured/main.f90:
   !   read_parameter -> Init -> set_control_para -> check_mesh_quality -> Init_flow
   ! Must be called AFTER MPI_Init and rank splitting (Global_Var my_id set).
   !
   ! comm:    the STRUCT_GROUP sub-communicator.  The pristine solver hardcodes
   !          MPI_COMM_WORLD everywhere; because the coupling driver assigns the
   !          struct group to global ranks 0..n_struct-1, group-local ids equal
   !          global ids, so those COMM_WORLD calls stay correct as long as
   !          my_id/Total_proc are re-initialised from the sub-communicator
   !          (otherwise partition logic waits for struct ranks that do not
   !          exist and the run hangs).
   ! ctlfile: path to the control.ec namelist file.  read_parameter accepts an
   !          optional filename, so no symlink/copy hack is needed.
   !---------------------------------------------------------------------------
   subroutine struct_solver_init(comm, ctlfile, force_restart)
      use mpi
      use Global_Var
      use mod_struct_grid, only: check_mesh_multigrid, set_control_para, check_mesh_quality
      use mod_struct_init, only: read_parameter, init, Init_flow, Iflag_init
      implicit none
      integer,          intent(in)           :: comm
      character(len=*), intent(in), optional :: ctlfile
      logical,          intent(in), optional :: force_restart
      integer :: ierr
      integer, parameter :: IBuffer_Size = 10000000
      real(PRE_EC), allocatable, save :: buf_mpi(:)   ! Bsend buffer (stays attached)

      ! group-local rank/size (equal to global ids by the rank-split convention)
      call MPI_Comm_rank(comm, my_id, ierr)
      call MPI_Comm_size(comm, Total_proc, ierr)
      Struct_Comm = comm   ! all struct-internal collectives use the group comm

      ! attach the Bsend buffer that Init_mpi would normally attach
      if ( .not. allocated(buf_mpi) ) then
         allocate( buf_mpi(IBuffer_Size) )
         call MPI_Buffer_attach(buf_mpi, 8*IBuffer_Size, ierr)
      end if

      if ( present(ctlfile) ) then
         call read_parameter(ctlfile)
      else
         call read_parameter
      end if
      if (my_id == 0) call check_mesh_multigrid
      call Init
      call set_control_para
      call check_mesh_quality
      ! Coupled restart: force continuation from flow3d.dat regardless of the
      ! control.ec setting (phase 10 joint restart).
      if ( present(force_restart) ) then
         if ( force_restart ) then
            Iflag_init = 1
            if ( my_id == 0 ) print*, '[struct] joint restart: Iflag_init=1 forced'
         end if
      end if
      call Init_flow

      if (my_id == 0) print*, "[struct] initialisation done, Kstep=", Mesh(1)%Kstep
   end subroutine struct_solver_init

   !---------------------------------------------------------------------------
   ! struct_solver_save -- write flow3d.dat + Step_mess.dat via the native
   ! output_flow (mod_struct_io), plus the node-interpolated Plot3D function
   ! file flow3d_node.dat (phase 10).  Collective over Struct_Comm.
   !---------------------------------------------------------------------------
   subroutine struct_solver_save
      use mod_struct_io, only: output_flow, output_flow_nodes
      implicit none
      call output_flow
      call output_flow_nodes
   end subroutine struct_solver_save

   !---------------------------------------------------------------------------
   ! struct_solver_step -- advance the structured solver by one time step.
   ! Only the finest mesh (nMesh=1) is advanced; multigrid is not used in
   ! the coupling loop (consistent with the uns side single-level SIMPLE).
   !---------------------------------------------------------------------------
   subroutine struct_solver_step()
      use mod_struct_solver, only: NS_Time_advance
      implicit none
      call NS_Time_advance(1)
   end subroutine struct_solver_step

   !---------------------------------------------------------------------------
   ! struct_extract_iface -- extract the inner-cell primitive state at every
   ! registered interface face (struct side, non-dimensional).
   !
   ! Output arrays are sized Num_Interface and ordered to match Interface_List.
   ! The caller converts non-dim -> SI via mod_interface_units.
   !
   ! Primitive state: rho, u, v, w, T, p  (non-dimensional, OpenCFD-EC units)
   !
   ! The exported state is the INTERFACE FACE state, i.e. the arithmetic
   ! mean of the last interior cell (ic,jc,kc) and the ghost cell
   ! (ig,jg,kg), mirroring the first-order face reconstruction used by the
   ! structured flux between those two cells.  The ghost cell holds the
   ! peer (unstructured) first-cell state.  Exporting the interior cell
   ! value alone while the peer applies it as a face DIRICHLET value
   ! introduces a systematic half-cell bias (observed as a ~10% velocity
   ! deficit frozen in the coupling fixed point).
   !---------------------------------------------------------------------------
   subroutine struct_extract_iface(rho, u, v, w, T, p, nfaces)
      use Global_Var
      use mod_interface, only: Interface_List, Num_Interface, PEER_STRUCT
      implicit none
      integer, intent(out) :: nfaces
      real(dp), allocatable, intent(out) :: rho(:), u(:), v(:), w(:), T(:), p(:)
      integer :: i, ic, jc, kc, ig, jg, kg, blk
      real(PRE_EC) :: d1, uu1, v1, w1, p1, T1
      real(PRE_EC) :: d2, uu2, v2, w2, p2, T2
      Type(Block_TYPE), pointer :: B

      nfaces = 0
      do i = 1, Num_Interface
         if (Interface_List(i)%solver == PEER_STRUCT) nfaces = nfaces + 1
      end do
      if (nfaces == 0) return

      allocate(rho(nfaces), u(nfaces), v(nfaces), w(nfaces), T(nfaces), p(nfaces))

      nfaces = 0
      do i = 1, Num_Interface
         if (Interface_List(i)%solver /= PEER_STRUCT) cycle
         nfaces = nfaces + 1
         blk = Interface_List(i)%block_no
         B => Mesh(1)%Block(blk)
         ic = Interface_List(i)%ic; jc = Interface_List(i)%jc; kc = Interface_List(i)%kc
         ig = Interface_List(i)%ig; jg = Interface_List(i)%jg; kg = Interface_List(i)%kg

         ! interior cell: conservative -> primitive
         d1  = B%U(1, ic, jc, kc)
         uu1 = B%U(2, ic, jc, kc) / d1
         v1  = B%U(3, ic, jc, kc) / d1
         w1  = B%U(4, ic, jc, kc) / d1
         p1  = (gamma - 1.d0) * (B%U(5, ic, jc, kc) - &
               0.5d0 * d1 * (uu1*uu1 + v1*v1 + w1*w1))
         T1  = p1 * gamma * Ma * Ma / d1

         ! ghost cell (peer state): conservative -> primitive
         d2  = B%U(1, ig, jg, kg)
         uu2 = B%U(2, ig, jg, kg) / d2
         v2  = B%U(3, ig, jg, kg) / d2
         w2  = B%U(4, ig, jg, kg) / d2
         p2  = (gamma - 1.d0) * (B%U(5, ig, jg, kg) - &
               0.5d0 * d2 * (uu2*uu2 + v2*v2 + w2*w2))
         T2  = p2 * gamma * Ma * Ma / d2

         ! face state = arithmetic mean of the two cell states
         rho(nfaces) = real(0.5d0*(d1+d2),    dp)
         u(nfaces)   = real(0.5d0*(uu1+uu2),  dp)
         v(nfaces)   = real(0.5d0*(v1+v2),    dp)
         w(nfaces)   = real(0.5d0*(w1+w2),    dp)
         T(nfaces)   = real(0.5d0*(T1+T2),    dp)
         p(nfaces)   = real(0.5d0*(p1+p2),    dp)
      end do
   end subroutine struct_extract_iface

   !---------------------------------------------------------------------------
   ! struct_set_iface_bc -- characteristic subsonic-outflow ghost state.
   !
   ! Only the peer face PRESSURE p_nd is imposed (Dirichlet-Neumann partition):
   ! the interface is a subsonic OUTFLOW for the compressible struct side, so
   ! exactly one physical condition (back pressure) is prescribed; all other
   ! quantities leave the domain along characteristics.  The ghost cell is
   ! built with the same linearised Riemann reflection used by
   ! boundary_Farfield for a subsonic outlet:
   !
   !   c1 = sqrt(gamma*p1/d1)                 (inner-cell sound speed)
   !   db = d1 + (pb-p1)/c1**2                (face density)
   !   ub = u1 + (p1-pb)/(d1*c1) * n_out      (face velocity, normal wave)
   !   p2 = 2*pb - p1 ; d2 = 2*db - d1 ;
   !   u2 = 2*ub - u1                         (ghost)
   !
   ! so the struct interface FACE pressure equals the unstructured cell-centre
   ! pressure p_nd, and the face velocity (arithmetic mean used by
   ! struct_extract_iface) is ub -- the mass flux dynamically relaxes until
   ! the inner pressure matches the peer level, instead of freezing a
   ! pressure jump.  Peer velocity/temperature arrays are accepted for
   ! interface signature compatibility but intentionally NOT imposed:
   ! imposing velocity as well over-specifies a subsonic outflow.
   !
   ! alpha (0..1) under-relaxes the imposed back pressure: pb is blended with
   ! the inner pressure, so alpha=0 degenerates to pure extrapolation
   ! (ghost = inner state, fully passive).  At low Mach number the acoustic
   ! gain 1/(rho*c) is large and the cold-start peer pressure swings would
   ! otherwise feed back as huge face-velocity jumps; the driver ramps alpha
   ! over iface_ramp coupling iterations.
   !
   ! Valid while the normal face velocity stays outflow-directed (the only
   ! regime the coupling cases use).
   !
   ! The ghost cell indices (ig,jg,kg) were stored at registration time
   ! (phase 6 addition to mod_interface).
   !---------------------------------------------------------------------------
   subroutine struct_set_iface_bc(rho_nd, u_nd, v_nd, w_nd, T_nd, p_nd, nfaces, alpha)
      use Global_Var
      use mod_interface, only: Interface_List, Num_Interface, PEER_STRUCT
      implicit none
      integer, intent(in) :: nfaces
      real(dp), intent(in) :: rho_nd(nfaces), u_nd(nfaces), v_nd(nfaces), &
                              w_nd(nfaces), T_nd(nfaces), p_nd(nfaces)
      real(dp), intent(in), optional :: alpha     ! back-pressure ramp 0..1 (def 1)
      integer :: i, ig, jg, kg, ic, jc, kc, blk, cnt
      real(PRE_EC) :: d1, u1, v1, w1, p1, c1, pb, db, ub(3), nvec(3)
      real(PRE_EC) :: d2, u2, v2, w2, p2, E2, a_ramp, p_peer
      Type(Block_TYPE), pointer :: B

      a_ramp = 1.0d0
      if ( present(alpha) ) a_ramp = min(1.0d0, max(0.0d0, real(alpha, PRE_EC)))

      cnt = 0
      do i = 1, Num_Interface
         if (Interface_List(i)%solver /= PEER_STRUCT) cycle
         cnt = cnt + 1
         if (cnt > nfaces) exit

         blk = Interface_List(i)%block_no
         B => Mesh(1)%Block(blk)
         ig = Interface_List(i)%ig; jg = Interface_List(i)%jg; kg = Interface_List(i)%kg
         ic = Interface_List(i)%ic; jc = Interface_List(i)%jc; kc = Interface_List(i)%kc

         ! inner cell primitive state (non-dimensional)
         d1 = B%U(1, ic, jc, kc)
         u1 = B%U(2, ic, jc, kc) / d1
         v1 = B%U(3, ic, jc, kc) / d1
         w1 = B%U(4, ic, jc, kc) / d1
         p1 = (gamma - 1.d0) * (B%U(5, ic, jc, kc) - &
              0.5d0 * d1 * (u1*u1 + v1*v1 + w1*w1))
         c1 = sqrt(gamma * p1 / d1)

         ! outward unit normal at the interface: take it from the block's
         ! own face-normal tables (same convention as boundary_Farfield),
         ! keyed by the stored face direction 1..6.
         select case ( Interface_List(i)%face )
         case (1)   ! i-
            nvec = (/ -B%ni1(ic,jc,kc), -B%ni2(ic,jc,kc), -B%ni3(ic,jc,kc) /)
         case (4)   ! i+
            nvec = (/  B%ni1(ic+1,jc,kc),  B%ni2(ic+1,jc,kc),  B%ni3(ic+1,jc,kc) /)
         case (2)   ! j-
            nvec = (/ -B%nj1(ic,jc,kc), -B%nj2(ic,jc,kc), -B%nj3(ic,jc,kc) /)
         case (5)   ! j+
            nvec = (/  B%nj1(ic,jc+1,kc),  B%nj2(ic,jc+1,kc),  B%nj3(ic,jc+1,kc) /)
         case (3)   ! k-
            nvec = (/ -B%nk1(ic,jc,kc), -B%nk2(ic,jc,kc), -B%nk3(ic,jc,kc) /)
         case (6)   ! k+
            nvec = (/  B%nk1(ic,jc,kc+1),  B%nk2(ic,jc,kc+1),  B%nk3(ic,jc,kc+1) /)
         case default
            nvec = real(Interface_List(i)%normal, PRE_EC)
         end select

         ! prescribed face pressure: ramp peer pressure onto inner pressure
         p_peer = real(p_nd(cnt), PRE_EC)
         pb = a_ramp * p_peer + (1.0d0 - a_ramp) * p1

         ! linearised characteristic reflection
         db = d1 + (pb - p1) / (c1*c1)
         ub(1) = u1 + (p1 - pb)/(d1*c1) * nvec(1)
         ub(2) = v1 + (p1 - pb)/(d1*c1) * nvec(2)
         ub(3) = w1 + (p1 - pb)/(d1*c1) * nvec(3)

         p2 = 2.d0*pb - p1
         d2 = 2.d0*db - d1
         d2 = max(d2, 0.05d0)               ! positivity guard during transients
         u2 = 2.d0*ub(1) - u1
         v2 = 2.d0*ub(2) - v1
         w2 = 2.d0*ub(3) - w1
         E2 = p2 / (gamma - 1.d0) + 0.5d0 * d2 * (u2*u2 + v2*v2 + w2*w2)

         B%U(1, ig, jg, kg) = d2
         B%U(2, ig, jg, kg) = d2 * u2
         B%U(3, ig, jg, kg) = d2 * v2
         B%U(4, ig, jg, kg) = d2 * w2
         B%U(5, ig, jg, kg) = E2
      end do
   end subroutine struct_set_iface_bc

   !---------------------------------------------------------------------------
   ! struct_solver_run -- run the structured solver to the end of its
   ! pseudo-time horizon as a STANDALONE solve (struct-only mode of the
   ! self-dispatching bin/mixnsolver).
   !
   ! Additive wrapper: struct_solver_init (init) + the time loop mirrored from
   ! src/structured/main.f90's `program main`.  It deliberately does NOT call
   ! MPI_Finalize -- the caller owns the MPI lifecycle (the pristine standalone
   ! program main does finalise, but a library routine must not, or the caller
   ! would get a double-finalise).
   !---------------------------------------------------------------------------
   subroutine struct_solver_run(comm, ctlfile)
      use mpi
      use Global_Var
      use mod_struct_solver, only: NS_Time_advance, NS_2stge_multigrid, &
                                   NS_3stge_multigrid, Filtering_oneMesh, output_Res
      use mod_struct_io, only: comput_force, output_flow, output_vt, &
                               Time_average, output_flow_average
      implicit none
      integer,          intent(in) :: comm
      character(len=*), intent(in), optional :: ctlfile

      call struct_solver_init( comm, ctlfile )

      if ( my_id == 0 ) print*, " Start ......"

      ! time advancement: single / double / triple grid, Euler or RK3
      do while( Mesh(1)%tt < t_end )
         if ( Num_Mesh .eq. 1 ) then
            call NS_Time_advance(1)                ! single-grid, one step
         else if ( Num_Mesh .eq. 2 ) then
            call NS_2stge_multigrid                ! two-grid multigrid
         else
            call NS_3stge_multigrid                ! three-grid multigrid
         end if

         ! optional smoothing for stability
         if ( Kstep_smooth > 0 ) then
            if ( mod(Mesh(1)%Kstep, Kstep_smooth) == 0 ) call Filtering_oneMesh(1)
         end if

         ! periodic force + residual output
         if ( mod(Mesh(1)%Kstep, Kstep_show) == 0 ) then
            call comput_force
            call output_Res(1)
         end if
         ! periodic field dump (flow3d.dat, PLOT3D format)
         if ( mod(Mesh(1)%Kstep, Kstep_Save) == 0 ) then
            call output_flow
            if ( If_debug == 1 .and. If_viscous == 1 .and. &
                 Iflag_turbulence_model /= 0 ) call output_vt
         end if
         ! time averaging
         if ( Kstep_average > 0 ) then
            if ( mod(Mesh(1)%Kstep, Kstep_average) == 0 ) call Time_average
            if ( mod(Mesh(1)%Kstep, Kstep_Save) == 0 ) call output_flow_average
         end if
      end do
   end subroutine struct_solver_run

end module mod_struct_driver
