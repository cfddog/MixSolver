!===============================================================================
! mod_control.f90 -- Control-parameter file reader (phase 3)
!
! Plain keyword = value format, one item per line; '#' starts a comment.
! Boundary conditions are given one per line:
!     bc = <zone> wall
!     bc = <zone> symmetry
!     bc = <zone> velocity-inlet <ux> <uy> <uz>
!     bc = <zone> velocity-inlet-parabolic <Umean> <span> <axis> <origin> [<T>]
!                  (fully-developed profile along the inward face normal:
!                   |u| = 6*Umean*xi*(1-xi), xi = (x_axis-origin)/span)
!     bc = <zone> pressure-outlet <p>
!     bc = <zone> outflow           (fully-developed: zero gradient, global
!                                    mass scaling enforces outflow = inflow)
! A moving-wall (lid) sub-region inside a wall zone is selected by a plane
! filter (faces whose centroid coordinate along <dir> equals <coord>):
!     lid = <zone> <dir 1|2|3> <coord> <ux> <uy> <uz>
!===============================================================================
module mod_uns_control
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh, only: mesh_t
   implicit none
   private
   public :: ctrl_t, read_control, read_mesh_scale, BC_NONE, BC_WALL, BC_SYMMETRY, &
             BC_VINLET, BC_VINLET_PARAB, BC_POUTLET, BC_FARFIELD, BC_SLIPWALL, &
             BC_MASSINLET, &
             BC_OUTFLOW, BC_INTERFACE, &
             bc_type_name, resolve_cell_zones, cell_zone_t, &
             CZ_AUTO, CZ_FLUID, CZ_POROUS

   integer, parameter :: BC_NONE     = 0
   integer, parameter :: BC_WALL     = 1
   integer, parameter :: BC_SYMMETRY = 2
   integer, parameter :: BC_VINLET   = 3
   integer, parameter :: BC_POUTLET  = 4
   integer, parameter :: BC_FARFIELD  = 5     ! pressure-far-field (Riemann-style)
   integer, parameter :: BC_SLIPWALL  = 6     ! slip wall (inviscid: u.n=0, free tangential)
   integer, parameter :: BC_MASSINLET = 7     ! mass-flux inlet: mdot kg/m^2/s,
                                              ! u_f = -(mdot/rho) n_outward
   integer, parameter :: BC_OUTFLOW   = 8     ! fully-developed outflow: zero
                                              ! normal gradient for u, p, T;
                                              ! face fluxes rescaled globally so
                                              ! total outflow = total inflow
   integer, parameter :: BC_INTERFACE = 9     ! coupling interface (phase 6):
                                              ! face velocity/pressure come from
                                              ! the peer solver via the exchange
                                              ! layer (stored in bc_t%iface_*).
   integer, parameter :: BC_VINLET_PARAB = 10 ! parabolic (fully-developed) inlet:
                                              ! |u_f| = 6*Umean*xi*(1-xi) applied
                                              ! along the INWARD face normal, with
                                              ! xi = (x(axis)-origin)/span clamped
                                              ! to [0,1].  The mean of the profile
                                              ! over the span is exactly Umean, so
                                              ! the volume flow equals that of the
                                              ! uniform 'velocity-inlet Umean'.

   integer, parameter :: MAXBC    = 32
   integer, parameter :: MAXCZ    = 32
   integer, parameter :: LINELEN  = 512

   ! cell-zone block type specification (phase 3 step B; porous props phase 12)
   ! Format in .control:
   !   cell_zone = <id|name> fluid
   !   cell_zone = <id|name> porous  [perm=..] [inertial=..] [porosity=..]
   !                          [k_s=..] [cp_s=..] [rho_s=..] [h_sf=..] [a_sf=..]
   type :: cell_zone_t
      character(len=32) :: name = ''     ! zone name (if matched by name)
      integer  :: id       = 0            ! zone id (if matched by id)
      integer  :: ztype    = 0            ! 1=fluid, 2=porous
      real(dp) :: perm     = 0.0_dp       ! isotropic permeability (m^2)
      real(dp) :: perm_xx  = 0.0_dp       ! diagonal tensor xx component (m^2);
      real(dp) :: perm_yy  = 0.0_dp       !   0 = fall back to scalar perm
      real(dp) :: perm_zz  = 0.0_dp       !   (axis-aligned anisotropy only;
      !                                        off-diagonal terms would need a
      !                                        block-coupled momentum assembly)
      real(dp) :: inertial = 0.0_dp       ! inertial (Forchheimer) coefficient (1/m)
      real(dp) :: porosity = 1.0_dp       ! void fraction (0..1)
      real(dp) :: disp_l   = 0.0_dp       ! longitudinal thermal dispersivity (m)
      real(dp) :: disp_t   = 0.0_dp       ! transverse thermal dispersivity (m)
      real(dp) :: bj_alpha = 0.0_dp       ! Beavers-Joseph slip coefficient at
                                          ! fluid/porous interface faces;
                                          ! 0 = disabled (stress continuity)
      ! --- solid-phase thermal properties (for porous zones, phase 12) ---
      real(dp) :: k_s      = 0.0_dp       ! solid thermal conductivity (W/m/K)
      real(dp) :: cp_s     = 0.0_dp       ! solid specific heat (J/kg/K)
      real(dp) :: rho_s    = 0.0_dp       ! solid density (kg/m^3)
      real(dp) :: h_sf     = 0.0_dp       ! fluid-solid heat transfer coeff (W/m^2/K)
      real(dp) :: a_sf     = 0.0_dp       ! specific surface area (m^2/m^3)
   end type cell_zone_t

   integer, parameter :: CZ_AUTO   = 0   ! cell_zone line supplies coefficients only;
                                         ! block type comes from a VC tag / fluid default
   integer, parameter :: CZ_FLUID  = 1
   integer, parameter :: CZ_POROUS = 2

   ! thermal model selection (phase 12)
   !   'lte'  = local thermal equilibrium (single T equation, effective props)
   !   'ltne' = local thermal non-equilibrium (fluid T_f + solid T_s)
   character(len=16), parameter :: THERM_LTE  = 'lte'
   character(len=16), parameter :: THERM_LTNE = 'ltne'

   ! one boundary-condition specification
   type :: bc_spec_t
      integer  :: zone     = 0          ! face zone id
      integer  :: btype    = BC_NONE
      real(dp) :: uvel(3)  = 0.0_dp     ! inlet velocity
      real(dp) :: pval     = 0.0_dp     ! outlet pressure
      real(dp) :: mdot     = 0.0_dp     ! mass-flux inlet (kg/m^2/s, positive in)
      ! --- parabolic (fully-developed) velocity inlet (BC_VINLET_PARAB) ---
      ! |u_f| = 6*umean*xi*(1-xi) along the inward face normal, with
      ! xi = (x(axis) - uorigin) / uspan clamped to [0,1] ("origin" is the
      ! coordinate where the profile vanishes, e.g. the channel wall y=0).
      real(dp) :: umean    = 0.0_dp     ! Umean of the profile (= uniform u)
      real(dp) :: uspan    = 0.0_dp     ! span over which the profile spans 0..1
      real(dp) :: uorigin   = 0.0_dp    ! coordinate where |u| = 0
      integer  :: uaxis    = 2          ! profile coordinate: 1=x, 2=y, 3=z
      logical  :: has_lid  = .false.    ! moving-wall plane filter present
      integer  :: lid_dir  = 0          ! 1=x, 2=y, 3=z
      real(dp) :: lid_coord = 0.0_dp
      real(dp) :: lid_vel(3) = 0.0_dp
      ! --- thermal boundary conditions (phase 10) ---
      ! ttype: 0=adiabatic (default), 1=fixed-temperature (Dirichlet),
      !        2=fixed-heat-flux (Neumann, q>0 heats the domain)
      integer  :: ttype    = 0          ! zone-level thermal bc type
      real(dp) :: tval     = 0.0_dp     ! fixed-temperature wall value
      real(dp) :: qval     = 0.0_dp     ! fixed-heat-flux wall value
      ! plane-filtered thermal bc (up to 4 planes per zone, mirrors lid)
      integer  :: ntbc_plane = 0        ! number of plane filters (0 = none)
      integer  :: tbc_pdir(4)   = 0     ! 1=x, 2=y, 3=z
      real(dp) :: tbc_pcoord(4) = 0.0_dp
      integer  :: tbc_pttype(4) = 0     ! plane thermal bc type
      real(dp) :: tbc_ptval(4) = 0.0_dp
      real(dp) :: tbc_pqval(4) = 0.0_dp
   end type bc_spec_t

   ! whole control container
   type :: ctrl_t
      real(dp) :: rho       = 1.0_dp
      real(dp) :: mu        = 0.01_dp
      real(dp) :: alpha_u   = 0.7_dp
      real(dp) :: alpha_p   = 0.3_dp
      integer  :: outer_max = 2000
      real(dp) :: outer_tol = 1.0e-6_dp
      real(dp) :: lin_tol   = 1.0e-8_dp
      integer  :: lin_max   = 200
      real(dp) :: conv_blend = 0.0_dp   ! 0 = 1st-order upwind, 1 = limited
                                        ! 2nd-order upwind (deferred correction)
      real(dp) :: nonorth_corr = 1.0_dp ! 0 = orthogonal diffusion only,
                                        ! 1 = full non-orthogonal correction
                                        ! (deferred, lagged; interior faces)
      character(len=16) :: ppe_precond = 'ic0'  ! PPE CG preconditioner:
                                        ! 'jacobi' or 'ic0' (incomplete Cholesky)
      character(len=16) :: out_format = 'vtu'  ! result file format:
                                        ! 'vtu' (VTK XML, default) or
                                        ! 'tecplot' (binary .plt via TecIO)
      ! --- transient PISO parameters (phase 8) ---
      ! transient = .false. selects steady SIMPLE (simple_run); when .true.,
      ! main dispatches to piso_run and the momentum assembly adds an implicit
      ! time term while skipping the steady under-relaxation (alpha_u and
      ! alpha_p are forced to 1 internally).
      ! time_scheme selects the temporal discretisation (phase 9):
      !   1 = implicit Euler  (1st order, default; unconditionally linear-stable
      !                          but CFL-limited here due to explicit F)
      !   2 = BDF2            (2nd order; (3u^{n+1}-4u^n+u^{n-1})/(2dt). Step 1
      !                          bootstraps with Euler since u^{n-1} is absent)
      logical  :: transient   = .false.
      integer  :: time_scheme = 1         ! 1=implicit Euler, 2=BDF2
      real(dp) :: dt          = 0.0_dp   ! time step size (must be > 0 if transient)
      integer  :: n_time_max  = 0        ! max number of time steps to advance
      integer  :: n_correct   = 2        ! PISO pressure-correction sweeps (>=1)
      integer  :: n_out_every = 0        ! write VTU every N steps (0 = final only)
      ! --- PIMPLE stabilisation (phase 11) ---
      ! PIMPLE = PISO + SIMPLE-like momentum under-relaxation + optional outer
      ! iterations per time step.  When transient=.true. and pimple=.true. the
      ! driver calls pimple_run instead of piso_run.  Momentum is relaxed with
      ! alpha_u (the same parameter as steady SIMPLE) while pressure is still
      ! corrected in full (alpha_p=1); n_outer_iter outer sweeps per time step
      ! give extra robustness for large dt / CFL > 1.  Default n_outer_iter=1
      ! recovers "PISO with momentum relaxation" (the most common PIMPLE form).
      logical  :: pimple      = .false.
      integer  :: n_outer_iter= 1        ! outer iterations per time step (>=1)
      ! --- energy / Boussinesq parameters (phase 10) ---
      ! Energy equation: rho*cp*(dT/dt + u.grad T) = div(k grad T).
      ! Passive scalar when boussinesq=.false.; buoyancy coupling when
      ! boussinesq=.true. adds -rho*beta*(T-Tref)*gravity(:)*V to momentum.
      real(dp) :: cp        = 1005.0_dp ! specific heat capacity
      real(dp) :: k_cond    = 0.026_dp  ! thermal conductivity
      real(dp) :: tref      = 0.0_dp    ! Boussinesq reference temperature
      real(dp) :: beta      = 0.0_dp    ! thermal expansion coefficient
      real(dp) :: gravity(3)= 0.0_dp    ! gravity vector (points toward earth)
      real(dp) :: body_force(3) = 0.0_dp ! uniform body force per unit volume
                                         ! (N/m^3), e.g. a mean pressure
                                         ! gradient surrogate in open channels
      logical  :: boussinesq = .false.  ! .true. enables buoyancy coupling
      ! thermal model (phase 12): 'lte' (default) or 'ltne'
      character(len=16) :: thermal_model = 'lte'
      ! initial temperature field (phase 12 natural-convection validation):
      !   init_T : uniform interior temperature at t=0 / iteration 0
      !   t_pert : amplitude of a smooth x-asymmetric perturbation used to
      !            break symmetry in Rayleigh-Benard (hot-bottom) setups:
      !            T(x) = init_T + t_pert*sin(pi*(x-xmin)/Lx)
      real(dp) :: init_T = 0.0_dp
      real(dp) :: t_pert = 0.0_dp
      ! --- restart parameters (MPI Step 7) ---
      ! restart = .true. enables restart mode: rank 0 reads the field dump
      !           via read_field_dump, MPI_Bcast-s it to all ranks, and
      !           scatter_global_to_local seeds each rank's local fld.
      ! restart_file   : path to the field dump file (required if restart)
      ! partition_file : optional pre-existing partition map.  If the map
      !                  matches the current nprocs and mesh size it is
      !                  reused (no METIS call); otherwise the mesh is
      !                  re-partitioned and the default <stem>.part.map
      !                  is overwritten.
      ! dump_file      : if non-empty, write the final field state to this
      !                  file at end of run (for chained restarts).
      logical  :: restart         = .false.
      character(len=256) :: restart_file   = ''
      character(len=256) :: partition_file = ''
      character(len=256) :: dump_file      = ''
      ! mesh_scale: multiplicative factor applied to the node coordinates right
      ! after read_cas (before connectivity/geometry).  Use 1.0e-3 for meshes
      ! exported in millimetres so the solver sees SI (metre) geometry.
      real(dp) :: mesh_scale = 1.0_dp
      ! inlet_ramp: number of SIMPLE outer iterations over which the
      ! mass-flow-inlet speed ramps linearly from 1/inlet_ramp to the full
      ! value (cold-start stabilisation; 0/1 = no ramp).
      integer  :: inlet_ramp = 1
      integer  :: nbc       = 0
      type(bc_spec_t) :: bc(MAXBC)
      ! --- cell-zone block types (phase 3 step B) ---
      ! Each cell zone may be flagged as 'fluid' (default) or 'porous'.
      integer  :: ncz = 0
      type(cell_zone_t) :: cz(MAXCZ)
   end type ctrl_t

contains

   !----------------------------------------------------------------------------
   ! read_mesh_scale -- lightweight pre-scan of the control file for the
   ! 'mesh_scale' key only.  Needed because the node-coordinate scaling must
   ! be applied right after read_cas (before connectivity/geometry), while the
   ! full read_control happens later in the init sequence.
   ! Returns scale=1.0 when the key is absent or the file cannot be opened.
   !----------------------------------------------------------------------------
   subroutine read_mesh_scale( filename, scale )
      character(len=*), intent(in)  :: filename
      real(dp),         intent(out) :: scale
      character(len=LINELEN) :: line, key, val
      integer :: u, ios, ie

      scale = 1.0_dp
      open( newunit = u, file = trim(filename), status = 'old', &
            action = 'read', iostat = ios )
      if ( ios /= 0 ) return

      do
         read( u, '(a)', iostat = ios ) line
         if ( ios /= 0 ) exit
         ie = scan( line, '#' )
         if ( ie > 0 ) line = line(:ie-1)
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle
         ie = scan( line, '=' )
         if ( ie == 0 ) cycle
         key = lowercase( adjustl( line(:ie-1) ) )
         val = adjustl( line(ie+1:) )
         if ( trim(key) == 'mesh_scale' ) then
            read( val, *, iostat = ios ) scale
            if ( ios /= 0 ) scale = 1.0_dp
            exit
         end if
      end do
      close( u )
   end subroutine read_mesh_scale

   !----------------------------------------------------------------------------
   ! Read and parse the control file
   !----------------------------------------------------------------------------
   subroutine read_control( filename, ctrl, ier )
      character(len=*), intent(in)  :: filename
      type(ctrl_t),     intent(out) :: ctrl
      integer,          intent(out) :: ier

      character(len=LINELEN) :: line, key, val
      integer :: u, ios, ie, ib

      ier  = 0
      ib   = 0
      open( newunit = u, file = trim(filename), status = 'old', &
            action = 'read', iostat = ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'ERROR: cannot open control file: ' // trim(filename)
         ier = 1
         return
      end if

      do
         read( u, '(a)', iostat = ios ) line
         if ( ios /= 0 ) exit

         ie = scan( line, '#' )              ! strip comments
         if ( ie > 0 ) line = line(:ie-1)
         line = adjustl( line )
         if ( len_trim(line) == 0 ) cycle

         ie = scan( line, '=' )
         if ( ie == 0 ) then
            write(*,'(a)') 'ERROR: control line without "=": ' // trim(line)
            ier = 2
            exit
         end if
         key = lowercase( adjustl( line(:ie-1) ) )
         val = adjustl( line(ie+1:) )

         select case ( trim(key) )
         case ( 'rho' );        read( val, *, iostat = ios ) ctrl%rho
         case ( 'mu' );         read( val, *, iostat = ios ) ctrl%mu
         case ( 'alpha_u' );    read( val, *, iostat = ios ) ctrl%alpha_u
         case ( 'alpha_p' );    read( val, *, iostat = ios ) ctrl%alpha_p
         case ( 'outer_max' );  read( val, *, iostat = ios ) ctrl%outer_max
         case ( 'outer_tol' );  read( val, *, iostat = ios ) ctrl%outer_tol
         case ( 'lin_tol' );    read( val, *, iostat = ios ) ctrl%lin_tol
         case ( 'lin_max' );    read( val, *, iostat = ios ) ctrl%lin_max
         case ( 'conv_blend' ); read( val, *, iostat = ios ) ctrl%conv_blend
         case ( 'nonorth_corr' ); read( val, *, iostat = ios ) ctrl%nonorth_corr
         case ( 'ppe_precond' )
            ctrl%ppe_precond = trim(adjustl( lowercase(val) ))
         case ( 'out_format', 'out_fmt' )
            ctrl%out_format = trim(adjustl( lowercase(val) ))
            if ( ctrl%out_format /= 'vtu' .and. &
                 ctrl%out_format /= 'tecplot' ) then
               write(*,'(a)') 'ERROR: out_format must be vtu or tecplot: ' &
                              // trim(val)
               ier = 2; exit
            end if

         ! --- transient PISO keywords (phase 8) ---
         case ( 'transient' )
            ! accept true/false, 1/0, yes/no
            select case ( trim(adjustl( lowercase(val) ) ) )
            case ( 'true', '1', 'yes', 'on' );  ctrl%transient = .true.
            case ( 'false', '0', 'no', 'off' ); ctrl%transient = .false.
            case default
               write(*,'(a)') 'ERROR: bad transient value (use true/false): ' &
                              // trim(val)
               ier = 2; exit
            end select
         case ( 'dt' );          read( val, *, iostat = ios ) ctrl%dt
         case ( 'n_time_max' ); read( val, *, iostat = ios ) ctrl%n_time_max
         case ( 'n_correct' );  read( val, *, iostat = ios ) ctrl%n_correct
         case ( 'n_out_every' );read( val, *, iostat = ios ) ctrl%n_out_every
         case ( 'pimple' )
            select case ( trim(adjustl( lowercase(val) ) ) )
            case ( 'true', '1', 'yes', 'on' );  ctrl%pimple = .true.
            case ( 'false', '0', 'no', 'off' ); ctrl%pimple = .false.
            case default
               write(*,'(a)') 'ERROR: bad pimple value (use true/false): ' &
                              // trim(val)
               ier = 2; exit
            end select
         case ( 'n_outer_iter' );read( val, *, iostat = ios ) ctrl%n_outer_iter
         case ( 'time_scheme' )
            read( val, *, iostat = ios ) ctrl%time_scheme
            if ( ios /= 0 ) then
               write(*,'(a)') 'ERROR: bad time_scheme value: ' // trim(val)
               ier = 2; exit
            end if
            if ( ctrl%time_scheme < 1 .or. ctrl%time_scheme > 2 ) then
               write(*,'(a,i0)') 'ERROR: time_scheme must be 1 (Euler) or 2 (BDF2), got ', &
                                 ctrl%time_scheme
               ier = 2; exit
            end if

         case ( 'bc' )
            if ( ctrl%nbc >= MAXBC ) then
               write(*,'(a,i0)') 'ERROR: too many bc lines (max ', MAXBC
               ier = 3; exit
            end if
            ctrl%nbc = ctrl%nbc + 1
            ib = ctrl%nbc
            call parse_bc( val, ctrl%bc(ib), ier )
            if ( ier /= 0 ) exit

         case ( 'lid' )
            call parse_lid( val, ctrl, ier )
            if ( ier /= 0 ) exit

         ! --- energy / Boussinesq keywords (phase 10) ---
         case ( 'cp' );         read( val, *, iostat = ios ) ctrl%cp
         case ( 'k_cond' );     read( val, *, iostat = ios ) ctrl%k_cond
         case ( 'tref' );       read( val, *, iostat = ios ) ctrl%tref
         case ( 'beta' );       read( val, *, iostat = ios ) ctrl%beta
         case ( 'gravity' );   read( val, *, iostat = ios ) ctrl%gravity
         case ( 'body_force' )
            read( val, *, iostat = ios ) ctrl%body_force
            if ( ios /= 0 ) then
               write(*,'(a)') 'ERROR: bad body_force value: ' // trim(val)
               ier = 9; return
            end if
         case ( 'boussinesq' )
            select case ( trim(adjustl( lowercase(val) ) ) )
            case ( 'true', '1', 'yes', 'on' );  ctrl%boussinesq = .true.
            case ( 'false', '0', 'no', 'off' ); ctrl%boussinesq = .false.
            case default
               write(*,'(a)') 'ERROR: bad boussinesq value (use true/false): ' &
                              // trim(val)
               ier = 2; exit
            end select

         case ( 'thermal_model' )
            select case ( trim(adjustl( lowercase(val) ) ) )
            case ( 'lte' );   ctrl%thermal_model = THERM_LTE
            case ( 'ltne' );  ctrl%thermal_model = THERM_LTNE
            case default
               write(*,'(a)') 'ERROR: thermal_model must be lte or ltne: ' &
                              // trim(val)
               ier = 2; exit
            end select

         case ( 'init_t' );     read( val, *, iostat = ios ) ctrl%init_T
         case ( 't_pert' );     read( val, *, iostat = ios ) ctrl%t_pert

         case ( 'tbc' )
            call parse_tbc( val, ctrl, ier )
            if ( ier /= 0 ) exit

         case ( 'tbc_plane' )
            call parse_tbc_plane( val, ctrl, ier )
            if ( ier /= 0 ) exit

         ! --- restart keywords (MPI Step 7) ---
         case ( 'restart' )
            select case ( trim(adjustl( lowercase(val) ) ) )
            case ( 'true', '1', 'yes', 'on' );  ctrl%restart = .true.
            case ( 'false', '0', 'no', 'off' ); ctrl%restart = .false.
            case default
               write(*,'(a)') 'ERROR: bad restart value (use true/false): ' &
                              // trim(val)
               ier = 2; exit
            end select
         case ( 'restart_file' );    ctrl%restart_file   = trim(adjustl(val))
         case ( 'partition_file' );  ctrl%partition_file = trim(adjustl(val))
         case ( 'dump_file' );       ctrl%dump_file      = trim(adjustl(val))
         case ( 'mesh_scale' );      read( val, *, iostat = ios ) ctrl%mesh_scale
         case ( 'inlet_ramp' );      read( val, *, iostat = ios ) ctrl%inlet_ramp

         ! --- cell-zone block types (phase 3 step B) ---
         ! cell_zone = <id|name> fluid
         ! cell_zone = <id|name> porous  [perm=..] [inertial=..] [porosity=..]
         case ( 'cell_zone' )
            if ( ctrl%ncz >= MAXCZ ) then
               write(*,'(a,i0)') 'ERROR: too many cell_zone lines (max ', MAXCZ
               ier = 8; exit
            end if
            ctrl%ncz = ctrl%ncz + 1
            call parse_cell_zone( val, ctrl%cz(ctrl%ncz), ier )
            if ( ier /= 0 ) exit

         case default
            write(*,'(a)') 'WARNING: unknown control keyword: ' // trim(key)
         end select

         if ( ier /= 0 ) exit
      end do
      close( u )

      if ( ier == 0 .and. ctrl%nbc == 0 ) then
         write(*,'(a)') 'ERROR: control file defines no boundary conditions'
         ier = 5
      end if

   end subroutine read_control

   !----------------------------------------------------------------------------
   ! Resolve "cell_zone" control entries against the CAS cell-zone table and
   ! mark every cell with its block type in m%cztype:
   !   CZ_FLUID (1)  -- ordinary fluid region
   !   CZ_POROUS (2) -- porous medium (coefficients are stored in ctrl%cz;
   !                    the momentum sink itself is implemented later)
   !
   ! Block-type resolution order (highest priority first):
   !   1. explicit "cell_zone = <id|name> fluid|porous" line in the control file
   !   2. "VC: porous" / "VC: fluid" tag embedded in the CAS zone name
   !      (e.g. 'Zone 2 ..., VC: porous Fluid = 1'); case-insensitive
   !   3. CZ_FLUID default
   ! A cell_zone line without a type token (key=value options only) carries
   ! CZ_AUTO: it supplies porous coefficients but never overrides the type.
   !
   ! Matching is by integer zone id when the control entry gave an id, or by
   ! zone name (cond_name/user_name, case-insensitive) when it gave a name.
   !
   ! verbose=.false. suppresses the report (MPI ranks other than rank 0).
   !----------------------------------------------------------------------------
   subroutine resolve_cell_zones( m, ctrl, ier, verbose )
      type(mesh_t),    intent(inout) :: m
      type(ctrl_t),    intent(in)    :: ctrl
      integer,         intent(out)   :: ier
      logical, optional, intent(in)  :: verbose

      integer :: i, ic, z, iz, nmatched, vctype, npor
      logical :: loud, found, typed
      logical, allocatable :: explicit_type(:)
      character(len=128) :: want

      ier = 0
      loud = .true.
      if ( present(verbose) ) loud = verbose

      if ( .not. allocated(m%czone) ) then
         write(*,'(a)') 'ERROR: mesh has no cell-zone map (old reader?)'
         ier = 10
         return
      end if

      ! default: every cell is fluid unless overridden below
      if ( allocated(m%cztype) ) deallocate( m%cztype )
      allocate( m%cztype(m%ncells), source = CZ_FLUID )
      allocate( explicit_type(m%ncells), source = .false. )

      if ( loud ) then
         write(*,'(a)') ''
         write(*,'(a)') '--- Cell (volume) zones ---'
         write(*,'(a)') '    id  condition          name                ncells  type    source'
      end if

      ! --- priority 1: explicit fluid|porous from the control file ----------
      do i = 1, ctrl%ncz
         found = .false.
         do iz = 1, m%nczone
            if ( ctrl%cz(i)%id /= 0 ) then
               found = ( m%czt(iz)%id == ctrl%cz(i)%id )
            else
               want = trim(adjustl( lowercase(ctrl%cz(i)%name) ))
               found = ( trim(adjustl(lowercase(m%czt(iz)%cond_name))) == trim(want) .or. &
                         trim(adjustl(lowercase(m%czt(iz)%user_name))) == trim(want) )
            end if
            if ( found ) then
               z = iz
               exit
            end if
         end do

         if ( .not. found ) then
            if ( ctrl%cz(i)%id /= 0 ) then
               write(*,'(a,i0)') 'ERROR: cell_zone refers to unknown zone id ', &
                                 ctrl%cz(i)%id
            else
               write(*,'(a)') 'ERROR: cell_zone refers to unknown zone name: ' // &
                              trim(ctrl%cz(i)%name)
            end if
            ier = 11
            return
         end if

         if ( ctrl%cz(i)%ztype /= CZ_AUTO ) then
            do ic = 1, m%ncells
               if ( m%czone(ic) == m%czt(z)%id ) then
                  m%cztype(ic)      = ctrl%cz(i)%ztype
                  explicit_type(ic) = .true.
               end if
            end do
         end if
      end do

      ! --- priority 2: VC: porous/fluid tags embedded in the CAS zone names --
      do iz = 1, m%nczone
         vctype = vc_tag_type( m%czt(iz)%cond_name, m%czt(iz)%user_name )
         if ( vctype == 0 ) cycle
         typed = .false.
         do ic = 1, m%ncells
            if ( m%czone(ic) == m%czt(iz)%id ) then
               if ( explicit_type(ic) ) then
                  typed = .true.
                  cycle
               end if
               m%cztype(ic) = vctype
            end if
         end do
         if ( typed .and. loud ) &
            write(*,'(a,i0,a)') '  note: zone ', m%czt(iz)%id, &
               ' carries a VC tag but its type was already set explicitly'
      end do

      if ( loud ) then
         do iz = 1, m%nczone
            nmatched = count_cells_type( m, m%czt(iz)%id, CZ_POROUS )
            npor     = nmatched
            if ( npor > 0 ) then
               call zone_source( iz, 'porous' )
            else
               call zone_source( iz, 'fluid'  )
            end if
         end do
         if ( ctrl%ncz == 0 ) &
            write(*,'(a)') '  (no cell_zone lines in control file; types from VC tags or fluid default)'
      end if

      ! --- sanity: porous zones without a coefficients entry ----------------
      ! perm defaults to 0, which silently switches the Darcy sink off.
      if ( loud ) then
         do iz = 1, m%nczone
            if ( count_cells_type( m, m%czt(iz)%id, CZ_POROUS ) == 0 ) cycle
            found = .false.
            do i = 1, ctrl%ncz
               if ( ctrl%cz(i)%id /= 0 ) then
                  found = ( m%czt(iz)%id == ctrl%cz(i)%id )
               else
                  want = trim(adjustl( lowercase(ctrl%cz(i)%name) ))
                  found = ( trim(adjustl(lowercase(m%czt(iz)%cond_name))) == trim(want) .or. &
                            trim(adjustl(lowercase(m%czt(iz)%user_name))) == trim(want) )
               end if
               if ( found ) exit
            end do
            if ( .not. found ) then
               write(*,'(a,i0,a)') 'WARNING: porous zone ', m%czt(iz)%id, &
                  ' has no cell_zone coefficients line (perm defaults to 0 -> no Darcy sink)'
            end if
         end do
      end if

      deallocate( explicit_type )

   contains

      ! Report one cell zone with the source of its block type.
      subroutine zone_source( zid, tname )
         integer, intent(in) :: zid
         character(len=*), intent(in) :: tname

         integer :: jc, nexpl, ntot
         character(len=10) :: src

         nexpl = 0; ntot = 0
         do jc = 1, m%ncells
            if ( m%czone(jc) == m%czt(zid)%id ) then
               ntot = ntot + 1
               if ( explicit_type(jc) ) nexpl = nexpl + 1
            end if
         end do
         if ( nexpl == ntot ) then
            src = 'control'
         else if ( vc_tag_type( m%czt(zid)%cond_name, m%czt(zid)%user_name ) /= 0 ) then
            src = 'VC tag'
         else
            src = 'default'
         end if
         write(*,'(i6,2x,a18,2x,a16,i8,2x,a6,2x,a)') &
            m%czt(zid)%id, trim(m%czt(zid)%cond_name), &
            trim(m%czt(zid)%user_name), ntot, trim(tname), trim(src)
      end subroutine zone_source

   end subroutine resolve_cell_zones

   !----------------------------------------------------------------------------
   ! Count cells of a given cell zone whose resolved type equals ztype.
   !----------------------------------------------------------------------------
   integer function count_cells_type( m, zid, ztype ) result( n )
      type(mesh_t), intent(in) :: m
      integer, intent(in) :: zid, ztype
      integer :: i
      n = 0
      do i = 1, m%ncells
         if ( m%czone(i) == zid .and. m%cztype(i) == ztype ) n = n + 1
      end do
   end function count_cells_type

   !----------------------------------------------------------------------------
   ! Detect a "VC: porous" / "VC: fluid" marker inside a CAS zone name.
   ! Both condition and user names are searched, case-insensitive.  The token
   ! following "vc" (separators ':' and blanks skipped) must be 'porous' or
   ! 'fluid' (trailing punctuation allowed); returns CZ_POROUS / CZ_FLUID,
   ! or 0 when no valid marker is present.
   !----------------------------------------------------------------------------
   integer function vc_tag_type( name1, name2 ) result( ztype )
      character(len=*), intent(in) :: name1, name2

      character(len=260) :: lc
      character(len=24) :: tok
      integer :: p, q, e, r, ll

      ztype = 0
      lc = trim(adjustl(lowercase(name1))) // ' ' // trim(adjustl(lowercase(name2)))
      ll = len_trim(lc)
      p = index( lc, 'vc' )
      do while ( p > 0 )
         ! word boundary before 'vc'
         if ( p == 1 .or. lc(p-1:p-1) < 'a' .or. lc(p-1:p-1) > 'z' ) then
            q = p + 2
            do while ( q <= ll .and. (lc(q:q) == ':' .or. lc(q:q) == ' ') )
               q = q + 1
            end do
            if ( q <= ll ) then
               e = q
               do while ( e <= ll .and. lc(e:e) /= ' ' )
                  e = e + 1
               end do
               tok = adjustl( lc(q:e-1) )
               do while ( len_trim(tok) > 0 .and. &
                          scan( tok(len_trim(tok):len_trim(tok)), ',;:=.' ) > 0 )
                  tok = tok(:len_trim(tok)-1)
               end do
               select case ( trim(tok) )
               case ( 'porous' )
                  ztype = CZ_POROUS
                  return
               case ( 'fluid' )
                  ztype = CZ_FLUID
                  return
               end select
            end if
         end if
         r = index( lc(p+1:), 'vc' )   ! next occurrence, global position
         p = merge( p + r, 0, r > 0 )
      end do
   end function vc_tag_type

   !----------------------------------------------------------------------------
   ! Parse a cell_zone line:
   !   "<id|name> [fluid|porous] [perm=..] [inertial=..] [porosity=..] ..."
   !
   ! The first token is either an integer zone id or a zone name.  The block
   ! type token ('fluid'/'porous') is OPTIONAL: when omitted (the next token is
   ! a key=value pair) ztype is CZ_AUTO and the line supplies coefficients only,
   ! leaving the type to a VC tag or the fluid default.  key=value pairs set the
   ! porous-medium coefficients (only meaningful for porous zones).
   !----------------------------------------------------------------------------
   subroutine parse_cell_zone( val, cz, ier )
      character(len=*), intent(in)  :: val
      type(cell_zone_t), intent(out) :: cz
      integer,          intent(out) :: ier

      character(len=32) :: zname
      character(len=LINELEN) :: rest, token, kkey, kval
      integer :: ios, ie, id, s0, s1

      ier = 0
      read( val, *, iostat = ios ) zname
      if ( ios /= 0 ) then
         write(*,'(a)') 'ERROR: bad cell_zone line: ' // trim(val)
         ier = 9
         return
      end if

      ! first token: integer -> zone id; otherwise -> zone name
      read( zname, *, iostat = ios ) id
      if ( ios == 0 ) then
         cz%id = id
         cz%name = ''
      else
         cz%name = trim(adjustl(zname))
         cz%id   = 0
      end if

      ! Defaults: block type CZ_AUTO (inferred from a VC tag / fluid default in
      ! resolve_cell_zones); coefficients default to "no porous medium".
      cz%ztype    = CZ_AUTO
      cz%perm     = 0.0_dp
      cz%inertial = 0.0_dp
      cz%porosity = 1.0_dp

      ! Locate the raw remainder after the first token (zone names with spaces
      ! are not supported on cell_zone lines).
      s0 = 1
      do while ( s0 <= len_trim(val) .and. val(s0:s0) == ' ' )
         s0 = s0 + 1
      end do
      s1 = s0
      do while ( s1 <= len_trim(val) .and. val(s1:s1) /= ' ' )
         s1 = s1 + 1
      end do
      rest = adjustl( val(s1:) )

      ! Optional block-type token: a bare word without '=' must be fluid|porous.
      ! A key=value token here means the type is omitted (CZ_AUTO).
      if ( len_trim(rest) > 0 ) then
         rest = adjustl(rest)
         ie = scan( rest, ' ' )
         if ( ie == 0 ) then
            token = rest
            rest  = ''
         else
            token = rest(:ie-1)
            rest  = adjustl(rest(ie+1:))
         end if
         if ( scan( token, '=' ) == 0 ) then
            select case ( trim(adjustl( lowercase(token) )) )
            case ( 'fluid' )
               cz%ztype = CZ_FLUID
            case ( 'porous' )
               cz%ztype = CZ_POROUS
            case default
               write(*,'(a)') 'ERROR: cell_zone type must be fluid or porous: ' // trim(token)
               ier = 9
               return
            end select
         else
            rest = trim(token) // ' ' // trim(rest)
         end if
      end if

      ! parse key=value coefficient pairs from the rest of the line
      do while ( len_trim(rest) > 0 )
            rest = adjustl( rest )
            ie = scan( rest, ' ' )
            if ( ie == 0 ) then
               token = rest
               rest = ''
            else
               token = rest(:ie-1)
               rest = rest(ie+1:)
            end if
            if ( len_trim(token) == 0 ) cycle

            ie = scan( token, '=' )
            if ( ie == 0 ) then
               write(*,'(a)') 'ERROR: bad cell_zone option (need key=value): ' // trim(token)
               ier = 9
               return
            end if
            kkey = trim(adjustl( lowercase(token(:ie-1)) ))
            kval = adjustl( token(ie+1:) )
            select case ( kkey )
            case ( 'perm' )
               read( kval, *, iostat = ios ) cz%perm
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad perm value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'perm_xx' )
               read( kval, *, iostat = ios ) cz%perm_xx
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad perm_xx value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'perm_yy' )
               read( kval, *, iostat = ios ) cz%perm_yy
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad perm_yy value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'perm_zz' )
               read( kval, *, iostat = ios ) cz%perm_zz
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad perm_zz value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'disp_l' )
               read( kval, *, iostat = ios ) cz%disp_l
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad disp_l value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'disp_t' )
               read( kval, *, iostat = ios ) cz%disp_t
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad disp_t value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'bj_alpha' )
               read( kval, *, iostat = ios ) cz%bj_alpha
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad bj_alpha value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'inertial' )
               read( kval, *, iostat = ios ) cz%inertial
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad inertial value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'porosity' )
               read( kval, *, iostat = ios ) cz%porosity
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad porosity value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'k_s' )
               read( kval, *, iostat = ios ) cz%k_s
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad k_s value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'cp_s' )
               read( kval, *, iostat = ios ) cz%cp_s
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad cp_s value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'rho_s' )
               read( kval, *, iostat = ios ) cz%rho_s
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad rho_s value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'h_sf' )
               read( kval, *, iostat = ios ) cz%h_sf
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad h_sf value: ' // trim(kval)
                  ier = 9; return
               end if
            case ( 'a_sf' )
               read( kval, *, iostat = ios ) cz%a_sf
               if ( ios /= 0 ) then
                  write(*,'(a)') 'ERROR: bad a_sf value: ' // trim(kval)
                  ier = 9; return
               end if
            case default
               write(*,'(a)') 'ERROR: unknown cell_zone option: ' // trim(kkey)
               ier = 9; return
            end select
      end do

   end subroutine parse_cell_zone

   !----------------------------------------------------------------------------
   ! Parse "zone type [value]" of a bc line
   !----------------------------------------------------------------------------
   subroutine parse_bc( val, spec, ier )
      character(len=*), intent(in)  :: val
      type(bc_spec_t),  intent(out) :: spec
      integer,          intent(out) :: ier

      character(len=32) :: bname
      integer :: ios

      ier = 0
      spec%btype = BC_NONE
      read( val, *, iostat = ios ) spec%zone, bname
      if ( ios /= 0 ) then
         write(*,'(a)') 'ERROR: bad bc line: ' // trim(val)
         ier = 6
         return
      end if

      select case ( trim(bname) )
      case ( 'wall' )
         spec%btype = BC_WALL
      case ( 'symmetry' )
         spec%btype = BC_SYMMETRY
      case ( 'velocity-inlet', 'velocity_inlet', 'inlet' )
         spec%btype = BC_VINLET
         read( val, *, iostat = ios ) spec%zone, bname, spec%uvel
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: velocity-inlet needs ux uy uz: ' // trim(val)
            ier = 6
         else
            ! optional 4th token: static temperature for the inlet Dirichlet T
            ! (bc_face_T anchors BC_VINLET to tval; leaving it at the default 0
            ! drains the domain towards T=0 in standalone runs)
            read( val, *, iostat = ios ) spec%zone, bname, spec%uvel, spec%tval
            ios = 0
         end if
      case ( 'velocity-inlet-parabolic', 'velocity_inlet_parabolic', &
             'parabolic-inlet', 'parabolic_inlet', 'inlet-parabolic' )
         ! Fully-developed (parabolic) velocity inlet:
         !   "<zone> velocity-inlet-parabolic <Umean> <span> <axis> <origin> [<T>]"
         ! The face velocity magnitude follows |u_f| = 6*Umean*xi*(1-xi) with
         ! xi = (x(axis) - origin)/span (clamped to [0,1]) and is applied
         ! normal to the patch, pointing INTO the domain (like mass-flow-inlet
         ! it uses the outward face normal: u_f = -|u_f| * n_outward).  Its
         ! span-average is exactly Umean, so the volume flow equals that of
         ! "velocity-inlet <Umean> 0 0 ..." on the same patch -- only the
         ! shape changes (this is the Betchen/paper fully-developed inlet).
         ! The optional last token is the inlet static temperature (same as
         ! velocity-inlet; defaults to 0, must be given when T is solved).
         spec%btype = BC_VINLET_PARAB
         read( val, *, iostat = ios ) spec%zone, bname, spec%umean, &
              spec%uspan, spec%uaxis, spec%uorigin, spec%tval
         if ( ios /= 0 ) then
            ! fewer than 7 records: retry without the temperature token
            ! (list-directed input assigns the leading items before it hits
            ! the end of record, but re-read for clarity/robustness)
            spec%tval = 0.0_dp
            read( val, *, iostat = ios ) spec%zone, bname, spec%umean, &
                 spec%uspan, spec%uaxis, spec%uorigin
            if ( ios /= 0 ) then
               write(*,'(a)') 'ERROR: velocity-inlet-parabolic needs ' // &
                  'Umean span axis origin: ' // trim(val)
               ier = 6
            end if
         end if
         if ( ier == 0 ) then
            if ( spec%uspan <= 0.0_dp ) then
               write(*,'(a)') 'ERROR: velocity-inlet-parabolic span must be > 0: ' &
                  // trim(val)
               ier = 6
            else if ( spec%uaxis < 1 .or. spec%uaxis > 3 ) then
               write(*,'(a)') 'ERROR: velocity-inlet-parabolic axis must be 1, 2 or 3: ' &
                  // trim(val)
               ier = 6
            end if
         end if
      case ( 'pressure-outlet', 'pressure_outlet', 'outlet' )
         spec%btype = BC_POUTLET
         read( val, *, iostat = ios ) spec%zone, bname, spec%pval
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: pressure-outlet needs p value: ' // trim(val)
            ier = 6
         end if
      case ( 'pressure-far-field', 'pressure_far_field', 'farfield', 'far-field' )
         ! pressure-far-field <p_far> <ux> <uy> <uz>
         ! Reuses uvel(3) for free-stream velocity and pval for free-stream pressure.
         ! At inflow (u.n < 0) velocity is fixed to u_far; at outflow velocity is
         ! extrapolated from the owner cell (zeroth-order). Pressure is always
         ! fixed to p_far (Dirichlet).
         spec%btype = BC_FARFIELD
         read( val, *, iostat = ios ) spec%zone, bname, spec%pval, spec%uvel
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: pressure-far-field needs p ux uy uz: ' // trim(val)
            ier = 6
         end if
      case ( 'slip-wall', 'slipwall', 'slip' )
         ! slip-wall: zero normal velocity (u.n=0), free tangential velocity.
         ! Used for inviscid flow simulations where the no-slip wall would
         ! create a singular boundary layer as mu -> 0.
         spec%btype = BC_SLIPWALL
         read( val, *, iostat = ios ) spec%zone, bname
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: slip-wall needs zone only: ' // trim(val)
            ier = 6
         end if
      case ( 'mass-flow-inlet', 'mass_flow_inlet', 'mass-inlet', 'mass_inlet' )
         ! mass-flux inlet: "<zone> mass-flow-inlet <mdot> [<Tin>]"
         ! mdot is mass flux per unit area [kg/m^2/s] INTO the domain; the
         ! face velocity is applied normal to each inlet face:
         ! u_f = -(mdot/rho) * n_outward (see bc_face_vel).  rho comes from
         ! the control file; for a planar inlet this gives a uniform inlet.
         ! The optional 4th token is the inlet static temperature (Dirichlet
         ! for the energy equation via bc_face_T); it defaults to 0 and must
         ! be given whenever the energy equation is solved.
         spec%btype = BC_MASSINLET
         read( val, *, iostat = ios ) spec%zone, bname, spec%mdot
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: mass-flow-inlet needs mdot value: ' // trim(val)
            ier = 6
         end if
         if ( spec%mdot <= 0.0_dp ) then
            write(*,'(a)') 'ERROR: mass-flow-inlet mdot must be positive (into domain)'
            ier = 6
         end if
         read( val, *, iostat = ios ) spec%zone, bname, spec%mdot, spec%tval
      case ( 'outflow', 'fully-developed-outlet', 'fully_developed_outlet' )
         ! Fully-developed outflow: zero normal gradient for velocity,
         ! pressure and temperature; no fixed pressure value.  The outflow
         ! face fluxes are rescaled every iteration (outflow_mass_sums /
         ! outflow_mass_scale in mod_uns_simple) so that the total outflow
         ! exactly matches the total inflow prescribed at fixed-flux
         ! boundaries -- this keeps the pure-Neumann PPE consistent.
         spec%btype = BC_OUTFLOW
         read( val, *, iostat = ios ) spec%zone, bname
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: outflow needs zone only: ' // trim(val)
            ier = 6
         end if
      case default
         write(*,'(a)') 'ERROR: unknown bc type: ' // trim(bname)
         ier = 7
      end select

   end subroutine parse_bc

   !----------------------------------------------------------------------------
   ! Parse "zone dir coord ux uy uz" of a lid line and attach it to the wall
   ! bc entry of that zone
   !----------------------------------------------------------------------------
   subroutine parse_lid( val, ctrl, ier )
      character(len=*), intent(in)    :: val
      type(ctrl_t),     intent(inout) :: ctrl
      integer,          intent(out)   :: ier

      integer :: ios, i, zid, ldir, ib
      real(dp) :: lcoord, lvel(3)

      ier = 0
      read( val, *, iostat = ios ) zid, ldir, lcoord, lvel
      if ( ios /= 0 .or. ldir < 1 .or. ldir > 3 ) then
         write(*,'(a)') 'ERROR: bad lid line (need zone dir coord ux uy uz): ' &
                        // trim(val)
         ier = 8
         return
      end if

      ! attach to the existing wall bc of this zone
      ib = 0
      do i = 1, ctrl%nbc
         if ( ctrl%bc(i)%zone == zid .and. ctrl%bc(i)%btype == BC_WALL ) ib = i
      end do
      if ( ib == 0 ) then
         write(*,'(a,i0)') 'ERROR: lid refers to a zone without a wall bc: ', &
                           zid
         ier = 8
         return
      end if

      ctrl%bc(ib)%has_lid   = .true.
      ctrl%bc(ib)%lid_dir   = ldir
      ctrl%bc(ib)%lid_coord = lcoord
      ctrl%bc(ib)%lid_vel   = lvel

   end subroutine parse_lid

   !----------------------------------------------------------------------------
   ! Parse "zone ttype [tval] [qval]" of a tbc line and attach it to the bc
   ! entry of that zone.  ttype: 0=adiabatic, 1=fixed-T (needs tval),
   ! 2=fixed-heat-flux (needs qval, q>0 heats the domain).
   !----------------------------------------------------------------------------
   subroutine parse_tbc( val, ctrl, ier )
      character(len=*), intent(in)    :: val
      type(ctrl_t),     intent(inout) :: ctrl
      integer,          intent(out)   :: ier

      integer :: ios, i, zid, tt, ib
      real(dp) :: tv, qv

      ier = 0
      read( val, *, iostat = ios ) zid, tt
      if ( ios /= 0 ) then
         write(*,'(a)') 'ERROR: bad tbc line (need zone ttype [tval/qval]): ' &
                        // trim(val)
         ier = 9; return
      end if
      if ( tt < 0 .or. tt > 2 ) then
         write(*,'(a,i0)') 'ERROR: tbc ttype must be 0/1/2, got ', tt
         ier = 9; return
      end if

      tv = 0.0_dp
      qv = 0.0_dp
      if ( tt == 1 ) then
         read( val, *, iostat = ios ) zid, tt, tv
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: tbc fixed-T needs tval: ' // trim(val)
            ier = 9; return
         end if
      else if ( tt == 2 ) then
         read( val, *, iostat = ios ) zid, tt, qv
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: tbc fixed-flux needs qval: ' // trim(val)
            ier = 9; return
         end if
      end if

      ib = 0
      do i = 1, ctrl%nbc
         if ( ctrl%bc(i)%zone == zid ) ib = i
      end do
      if ( ib == 0 ) then
         write(*,'(a,i0)') 'ERROR: tbc refers to a zone without a bc: ', zid
         ier = 9; return
      end if

      ctrl%bc(ib)%ttype = tt
      ctrl%bc(ib)%tval  = tv
      ctrl%bc(ib)%qval  = qv

   end subroutine parse_tbc

   !----------------------------------------------------------------------------
   ! Parse "zone dir coord ttype [tval] [qval]" of a tbc_plane line and attach
   ! it to the wall bc entry of that zone (plane-filtered thermal bc, mirrors
   ! the lid mechanism).
   !----------------------------------------------------------------------------
   subroutine parse_tbc_plane( val, ctrl, ier )
      character(len=*), intent(in)    :: val
      type(ctrl_t),     intent(inout) :: ctrl
      integer,          intent(out)   :: ier

      integer :: ios, i, zid, ldir, tt, ib, n
      real(dp) :: lcoord, tv, qv

      ier = 0
      read( val, *, iostat = ios ) zid, ldir, lcoord, tt
      if ( ios /= 0 .or. ldir < 1 .or. ldir > 3 .or. tt < 0 .or. tt > 2 ) then
         write(*,'(a)') 'ERROR: bad tbc_plane line (need zone dir coord ttype ' &
                        // '[tval/qval]): ' // trim(val)
         ier = 10; return
      end if

      tv = 0.0_dp
      qv = 0.0_dp
      if ( tt == 1 ) then
         read( val, *, iostat = ios ) zid, ldir, lcoord, tt, tv
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: tbc_plane fixed-T needs tval: ' // trim(val)
            ier = 10; return
         end if
      else if ( tt == 2 ) then
         read( val, *, iostat = ios ) zid, ldir, lcoord, tt, qv
         if ( ios /= 0 ) then
            write(*,'(a)') 'ERROR: tbc_plane fixed-flux needs qval: ' // trim(val)
            ier = 10; return
         end if
      end if

      ib = 0
      do i = 1, ctrl%nbc
         if ( ctrl%bc(i)%zone == zid .and. ctrl%bc(i)%btype == BC_WALL ) ib = i
      end do
      if ( ib == 0 ) then
         write(*,'(a,i0)') 'ERROR: tbc_plane refers to a zone without a wall bc: ', &
                           zid
         ier = 10; return
      end if

      n = ctrl%bc(ib)%ntbc_plane + 1
      if ( n > 4 ) then
         write(*,'(a,i0)') 'ERROR: too many tbc_plane for zone ', zid
         ier = 10; return
      end if
      ctrl%bc(ib)%ntbc_plane  = n
      ctrl%bc(ib)%tbc_pdir(n)   = ldir
      ctrl%bc(ib)%tbc_pcoord(n) = lcoord
      ctrl%bc(ib)%tbc_pttype(n) = tt
      ctrl%bc(ib)%tbc_ptval(n)  = tv
      ctrl%bc(ib)%tbc_pqval(n)  = qv

   end subroutine parse_tbc_plane

   !----------------------------------------------------------------------------
   ! Lower-case copy of a string
   !----------------------------------------------------------------------------
   function lowercase( s ) result( r )
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, ic
      r = s
      do i = 1, len_trim(s)
         ic = iachar( s(i:i) )
         if ( ic >= iachar('A') .and. ic <= iachar('Z') ) &
            r(i:i) = achar( ic + 32 )
      end do
   end function lowercase

   !----------------------------------------------------------------------------
   ! Name of a bc type (for printing)
   !----------------------------------------------------------------------------
   function bc_type_name( btype ) result( name )
      integer, intent(in) :: btype
      character(len=26) :: name
      select case ( btype )
      case ( BC_WALL );      name = 'wall'
      case ( BC_SYMMETRY );  name = 'symmetry'
      case ( BC_VINLET );    name = 'velocity-inlet'
      case ( BC_VINLET_PARAB ); name = 'velocity-inlet-parabolic'
      case ( BC_POUTLET );   name = 'pressure-outlet'
      case ( BC_FARFIELD );  name = 'pressure-far-field'
      case ( BC_SLIPWALL );  name = 'slip-wall'
      case ( BC_MASSINLET ); name = 'mass-flow-inlet'
      case ( BC_OUTFLOW );   name = 'outflow'
      case ( BC_INTERFACE ); name = 'coupling-interface'
      case default;          name = 'none'
      end select
   end function bc_type_name

end module mod_uns_control
