!===============================================================================
! mod_uns_driver.f90 -- callable wrapper for the unstructured solver (phase 6).
!
! Extracts the initialisation + step logic from main_uns.f90 so that the
! coupling driver (src/main.f90) can drive the unstructured solver inside a
! weak-coupling iteration loop:
!
!   uns_solver_init(casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier)
!   uns_solver_step(m, c, g, ctrl, bcs, fld, nsteps, ier)   ! SIMPLE N steps
!   uns_solver_extract_iface(m, g, bcs, fld, faces, rho, u, T, p, ier)
!
! The standalone main_uns.f90 keeps working unchanged (it is a thin wrapper
! around these routines plus the post-run reporting / output).
!===============================================================================
module mod_uns_driver
   use mod_precision, only: dp
   use mod_uns_mesh, only: mesh_t
   use mod_uns_connectivity, only: conn_t
   use mod_uns_geometry, only: geom_t
   use mod_uns_control, only: ctrl_t, BC_INTERFACE
   use mod_uns_bc, only: bc_t
   use mod_uns_fields, only: fields_t, init_fields, setup_porous_fields
   use mod_uns_cas_reader, only: read_cas
   use mod_uns_connectivity, only: build_connectivity
   use mod_uns_geometry, only: compute_geometry, register_interface_zones
   use mod_uns_control, only: read_control, resolve_cell_zones, read_mesh_scale
   use mod_uns_bc, only: build_bc
   use mod_uns_simple, only: simple_run
   ! The coupled restart reader lives in the MPI-only mod_uns_restart stack
   ! (-> mod_uns_mpi_core / mod_uns_partition / mod_uns_local_mesh), so it is
   ! pulled in only for the MPI tree; the serial tree falls back to the
   ! 'restart unavailable' warning at the call site below.
#ifdef HAVE_MPI
   use mod_uns_restart, only: read_field_dump
#endif
   implicit none
   private

   public :: uns_solver_init
   public :: uns_solver_step
   public :: uns_solver_extract_iface

contains

   !---------------------------------------------------------------------------
   ! Full initialisation: read mesh, build connectivity/geometry, register
   ! interface zones, read control, resolve cell zones, build BCs, init fields.
   ! Mirrors the pre-solver part of main_uns.f90.
   !---------------------------------------------------------------------------
   subroutine uns_solver_init( casfile, ctlfile, m, c, g, ctrl, bcs, fld, ier, &
                               restart_file )
      character(len=*), intent(in)  :: casfile, ctlfile
      type(mesh_t),     intent(out) :: m
      type(conn_t),     intent(out) :: c
      type(geom_t),     intent(out) :: g
      type(ctrl_t),     intent(out) :: ctrl
      type(bc_t),       intent(out) :: bcs
      type(fields_t),   intent(out) :: fld
      integer,          intent(out) :: ier
      character(len=*), intent(in), optional :: restart_file

      real(dp) :: mscale
      logical  :: rex
#ifdef HAVE_MPI
      integer  :: ier2
      character(len=512) :: src_file
#endif

      ier = 0

      call read_cas( trim(casfile), m, ier )
      if ( ier /= 0 ) return

      ! optional unit scaling (e.g. mm -> m) before connectivity/geometry
      call read_mesh_scale( trim(ctlfile), mscale )
      if ( mscale /= 1.0_dp ) then
         m%x = m%x * mscale
         write(*,'(a,es10.3)') '[uns] mesh coordinates scaled by ', mscale
      end if

      call build_connectivity( m, c, ier )
      if ( ier /= 0 ) return

      call compute_geometry( m, g, ier )
      if ( ier /= 0 ) return

      call register_interface_zones( m, g )

      call read_control( trim(ctlfile), ctrl, ier )
      if ( ier /= 0 ) return

      call resolve_cell_zones( m, ctrl, ier )
      if ( ier /= 0 ) return

      call build_bc( m, ctrl, bcs, ier )
      if ( ier /= 0 ) return

      call init_fields( m, g, fld )
      call setup_porous_fields( m, ctrl, fld )

      ! Optional uniform initial temperature (+ x-asymmetric perturbation).
      ! Mirrors main_uns.f90: the coupled driver needs a physically correct
      ! initial T because the interface temperature is exchanged (a T=0
      ! field would send vacuum-state data to the structured ghost cells).
      if ( ctrl%init_T /= 0.0_dp .or. ctrl%t_pert /= 0.0_dp ) then
         block
            real(dp) :: xmin, xmax, lx, pi_v
            integer  :: ic
            pi_v = acos( -1.0_dp )
            xmin = minval( m%x(1,:) )
            xmax = maxval( m%x(1,:) )
            lx   = max( xmax - xmin, epsilon(1.0_dp) )
            do ic = 1, m%ncells
               fld%T(ic) = ctrl%init_T + ctrl%t_pert &
                           * sin( pi_v * ( g%xc(1,ic) - xmin ) / lx )
               fld%T_s(ic) = fld%T(ic)
            end do
            fld%T_old      = fld%T
            fld%T_old_old  = fld%T
            fld%T_s_old    = fld%T_s
            fld%T_s_old_old = fld%T_s
         end block
      end if

      ! ---- coupled restart: overwrite the initialised fields with the dump --
      ! The coupled uns group holds the FULL mesh on every rank (the mixed
      ! driver does not partition the unstructured mesh), so local == global and
      ! the dump maps straight onto fld -- no scatter step is needed.
      ! Serial read (serial=.true., no MPI_Bcast): the coupled driver never
      ! calls mpi_bootstrap, so mpi_comm is MPI_COMM_WORLD and a collective read
      ! would deadlock against the struct ranks.  Each uns rank reads the same
      ! file independently.  T_s is not stored in dump v1 -- LTNE cases restart
      ! T_s from init.
      if ( present(restart_file) ) then
         inquire( file=trim(restart_file), exist=rex )
         if ( rex ) then
#ifdef HAVE_MPI
            ier2 = 0
            call read_field_dump( trim(restart_file), m%ncells, fld, src_file, &
                                  ier2, serial=.true. )
            if ( ier2 /= 0 ) then
               write(*,'(a,a,a,i0)') '[uns] restart read FAILED: ', &
                                     trim(restart_file), ' ier=', ier2
               ier = ier2
               return
            end if
            write(*,'(a,a)') '[uns] restarted from dump: ', trim(restart_file)
#else
            ! Serial tree: mod_uns_restart is MPI-only and not linked here, so a
            ! coupled restart cannot be honoured -- warn and keep the init state.
            write(*,'(a)') '[uns] WARNING: restart unavailable in the serial '// &
                           '(non-HAVE_MPI) build; keeping initialised fields'
#endif
         else
            write(*,'(a,a)') '[uns] restart requested but dump missing: ', &
                             trim(restart_file)
            ier = -1
            return
         end if
      end if

   end subroutine uns_solver_init

   !---------------------------------------------------------------------------
   ! Advance the SIMPLE solver by nsteps outer iterations (or until convergence
   ! if that happens first).  Wraps simple_run with the nsteps cap.
   !---------------------------------------------------------------------------
   subroutine uns_solver_step( m, c, g, ctrl, bcs, fld, nsteps, ier )
      type(mesh_t),    intent(in)    :: m
      type(conn_t),    intent(in)    :: c
      type(geom_t),    intent(in)    :: g
      type(ctrl_t),    intent(in)    :: ctrl
      type(bc_t),      intent(in)    :: bcs
      type(fields_t),  intent(inout) :: fld
      integer,         intent(in)    :: nsteps
      integer,         intent(out)   :: ier
      call simple_run( m, c, g, ctrl, bcs, fld, ier, nsteps )
   end subroutine uns_solver_step

   !---------------------------------------------------------------------------
   ! Extract SI state on the coupling-interface faces.
   !
   ! For each BC_INTERFACE face the owner-cell state is taken as the face
   ! state (zero-gradient extrapolation -- the unstructured side sends its
   ! near-interface cell state to the structured side).
   !
   ! Output arrays are sized to the number of interface faces and returned
   ! alongside the face index list so the exchange layer can map them.
   !---------------------------------------------------------------------------
   subroutine uns_solver_extract_iface( m, g, bcs, fld, faces, rho, u, T, p, ier, rho_in )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld
      integer, allocatable, intent(out) :: faces(:)
      real(dp), allocatable, intent(out) :: rho(:), T(:), p(:), u(:,:)
      integer,        intent(out)   :: ier
      real(dp),       intent(in), optional :: rho_in   ! SI density (ctrl%rho)

      integer :: i, nif, gi, c0
      real(dp) :: rho_const

      ier = 0
      rho_const = 1.0_dp   ! safe default (water-like); the driver passes
                           ! ctrl%rho for the actual case fluid.
      if ( present(rho_in) ) rho_const = rho_in

      ! count interface faces
      nif = 0
      do i = 1, bcs%nb
         if ( bcs%gb(i)%btype == BC_INTERFACE ) nif = nif + bcs%gb(i)%nf
      end do

      allocate( faces(nif), rho(nif), T(nif), p(nif), u(3, nif) )
      nif = 0
      do gi = 1, bcs%nb
         if ( bcs%gb(gi)%btype /= BC_INTERFACE ) cycle
         do i = 1, bcs%gb(gi)%nf
            nif = nif + 1
            faces(nif) = bcs%gb(gi)%faces(i)
            c0 = m%f(faces(nif))%c0
            rho(nif) = rho_const
            u(:,nif) = fld%u(:,c0)
            T(nif)   = fld%T(c0)
            p(nif)   = fld%p(c0)
         end do
      end do

   end subroutine uns_solver_extract_iface

end module mod_uns_driver
