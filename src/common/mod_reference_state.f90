!===============================================================================
! mod_reference_state.f90 -- runtime-configurable reference state for the
! coupling interface (phase 5).
!
! The structured solver (OpenCFD-EC) works internally with non-dimensional
! variables; the unstructured solver works in SI.  At the coupling interface
! every quantity is converted to SI, exchanged, and converted back.  The
! scaling factors come from a *reference state* (rho_ref, T_ref, L_ref, ...).
!
! mod_constants.f90 holds compile-time defaults (sea-level standard air);
! this module provides a runtime-mutable copy that the coupling control file
! (mix.control, key=value style) can override per case.
!
! Non-dimensionalisation convention (OpenCFD-EC, as implemented in
! mod_struct_init/mod_struct_bc: u_inf* = 1, p_inf* = 1/(gamma*Ma^2),
! T* = gamma*Ma^2*p*/rho*):
!   rho* = rho / rho_ref
!   u*   = u   / U_inf        (U_inf = Ma * a_ref, a_ref = sqrt(gamma*R*T_ref))
!   T*   = T   / T_ref
!   p*   = p   / (rho_ref * U_inf**2)
! U_inf is the runtime key u_ref (0 -> a_ref, i.e. the Ma=1 special case).
!
! The reference state is process-global (shared by every rank); the coupling
! driver calls init_reference_state() once after reading mix.control.
!===============================================================================
module mod_reference_state
   use mod_precision, only: dp
   use mod_constants, only: GAMMA, R_GAS, RHO_REF, T_REF, L_REF
   implicit none
   private

   public :: reference_state_t
   public :: init_reference_state, read_mix_control
   public :: get_ref_state, get_coupling_params
   public :: scheduled_steps

   ! ---------------------------------------------------------------------------
   ! reference_state_t -- runtime copy of the scaling factors
   ! ---------------------------------------------------------------------------
   type :: reference_state_t
      real(dp) :: rho_ref = RHO_REF   ! [kg/m^3]
      real(dp) :: T_ref   = T_REF     ! [K]
      real(dp) :: L_ref   = L_REF     ! [m]
      real(dp) :: u_ref   = 0.0_dp    ! [m/s] free-stream U_inf (0 => a_ref)
      ! derived (set by compute_derived)
      real(dp) :: a_ref   = 0.0_dp    ! [m/s] speed of sound at T_ref
      real(dp) :: p_ref   = 0.0_dp    ! [Pa]   rho_ref * R * T_ref
      real(dp) :: p_scale = 0.0_dp    ! [Pa]   rho_ref * u_ref^2
   end type reference_state_t

   ! module-global instance
   type(reference_state_t), save :: g_ref

   ! ---------------------------------------------------------------------------
   ! Coupling-loop parameters (phase 7), also read from mix.control:
   !   n_couple     : number of weak-coupling iterations
   !   n_uns_steps  : SIMPLE outer steps per coupling iteration
   !   iface_relax  : under-relaxation factor for the interface BC applied on
   !                  the unstructured side:  BC = omega*struct + (1-omega)*uns,
   !                  omega = iface_relax * min(1, iter/iface_ramp)
   !   iface_ramp   : number of coupling iterations over which omega ramps
   !                  linearly from iface_relax/iface_ramp up to iface_relax
   ! ---------------------------------------------------------------------------
   integer,  save :: g_n_couple    = 3
   integer,  save :: g_n_uns_steps = 50
   integer,  save :: g_n_struct_steps = 1
   real(dp), save :: g_iface_relax = 1.0_dp
   integer,  save :: g_iface_ramp  = 1
   ! phase-13 start-up step schedule (see scheduled_steps below).  When
   ! n_*_steps_start > n_*_steps the corresponding side starts with the large
   ! count and halves it every g_step_decay_every coupling iterations until it
   ! reaches the steady count.  0 (or <= target) = schedule disabled.
   integer,  save :: g_n_struct_steps_start = 0
   integer,  save :: g_n_uns_steps_start    = 0
   integer,  save :: g_step_decay_every     = 1
   ! Interface velocity treatment used by the uns side when it imposes the
   ! peer (structured) interface state:
   !   0 = 'full'   : impose the whole peer velocity vector (legacy).
   !   1 = 'normal' : impose only the interface-normal component supplied by
   !                  the peer (mass-flux / Dirichlet-Neumann partition) and
   !                  take the tangential components from the local uns
   !                  interface state (zero-gradient).  This is the
   !                  physically consistent condition for a low-permeability
   !                  porous medium: the Darcian resistance cannot sustain the
   !                  free-stream tangential slip (~100 m/s), and forcing it
   !                  demands dp/dx = mu*u/K ~ 5e5 Pa/m -> divergence.
   !   2 = 'balance': ignore the peer velocity entirely and treat the interface
   !                  as a fully-developed outflow of the slab (own previous
   !                  interface state, zero gradient); the total interface
   !                  flux is then rescaled by the existing mass-balance
   !                  correction so that it exactly matches the coolant
   !                  injection.  Avoids the conflict between the peer normal
   !                  velocity and the flux correction (mode 1 imposes both on
   !                  the same quantity and the two fight each other).
   integer,  save :: g_iface_vel_mode = 0
   ! Interface pressure datum anchoring (opt-in, default off):
   ! The uns PPE has no pressure boundary (pure Neumann: the slab is closed
   ! except the coolant inlet and the coupling interface), so its gauge
   ! pressure level is an arbitrary free parameter that drifts during the
   ! run.  The struct, however, uses (uns gauge + p_ref) as an ABSOLUTE
   ! back pressure, so an un-anchored datum makes the two sides' interface
   ! pressures diverge (struct sucks / pressurises) and the loop blows up.
   ! When enabled, the whole uns gauge field is shifted every coupling
   ! iteration so that the mean interface absolute pressure matches the
   ! struct's -- a pure datum shift, the gradients (and hence the uns
   ! solution) are unchanged.
   integer,  save :: g_iface_p_anchor = 0
   ! Interface PRESSURE model (phase 14, opt-in; default 0 = legacy):
   !   0 = 'grad0'     : interface pressure is zero-gradient on the uns side
   !                     (bc_face_p returns pP).  The slab's PPE is then pure
   !                     Neumann (no pressure BC anywhere), so its gauge datum
   !                     is a free parameter -> iface_p_anchor / the uniform
   !                     interface-flux correction (du_n) exist to patch that.
   !   1 = 'dirichlet' : the peer (struct) interface pressure is imposed as a
   !                     per-face DIRICHLET pressure on the coupling interface
   !                     (converted to gauge).  The PPE becomes well posed
   !                     without pinning cell 1, the interface mass flux is
   !                     *solved* by continuity instead of being imposed and
   !                     patched, and iface_p_anchor/du_n are bypassed.  The
   !                     interface VELOCITY then must not be Dirichlet-ised on
   !                     the same faces (over-constrained): the normal component
   !                     switches to the zero-gradient/outflow treatment.
   !   2 = 'momentum'  : 'dirichlet' plus the Betchen 2006 Eq.43/44 normal
   !                     momentum balance across the flow-area change
   !                     (1-eps)/eps * rho*u_n^2, the two-sided inverse-distance
   !                     blend with the porous-side extrapolation
   !                     (weight g_iface_p_blend) and the deferred p-mdot
   !                     sub-iteration (g_iface_pm_subiter) of Betchen Sec.4.2.
   integer,  save :: g_iface_p_model = 0
   ! Interface TANGENTIAL (shear) treatment (phase 14, opt-in; 0 = off):
   !   1 = 'bj'     : Beavers-Joseph / Ochoa-Tapia-Whitaker stress jump at the
   !                  coupling interface through the series-resistance law
   !                  C = mu*A/(d_f + sqrt(K)/alpha + eps*d_p) of mod_iface_law
   !                  (the same law the solver applies to internal fluid/porous
   !                  faces via bj_alpha; d_f is the peer gap, d_p the slab
   !                  gap).
   !   2 = 'noslip' : tangential velocity CONTINUITY at the interface (Betchen
   !                  2006 Eq.13: u_fl = <u>_por) with the flush (no-slip)
   !                  porous-side conductance C = mu_e*A/d_p.  This is the
   !                  natural companion of iface_p_model (which supplies the
   !                  normal/PART of the interface constraint, so the
   !                  tangential part must be closed explicitly, Eq.16+Eq.13).
   integer,  save :: g_iface_slip = 0
   integer,  save :: g_iface_pm_subiter = 0      ! p-mdot sub-iterations
                                                 ! (0 = pure deferred: one
                                                 !  Eq.44 evaluation on the
                                                 !  current interface flux;
                                                 !  >=1 = damped local
                                                 !  sub-iteration with
                                                 !  omega = 1/(1+G))
   real(dp), save :: g_iface_p_blend    = 0.5_dp ! weight of the porous-side
                                                 ! interface-pressure estimate
   real(dp), save :: g_iface_slip_alpha = -1.0_dp ! BJ alpha for the coupling
                                                 ! interface; <0 => take the
                                                 ! mean bj_alpha of the adjacent
                                                 ! slab cells, else 1.0
   ! Zhang 2011 interface closure (Sec.3.5.1; phase 14, opt-in).
   !   0 = 'off' (default): legacy exchange (peer face state imposed directly).
   !   1 = 'zhang'        : the interface VELOCITY is evaluated from Eqs.22/25/26
   !                        -- the conductance-weighted normal balance and the
   !                        tangential stress-jump balance -- and imposed as a
   !                        full-vector Dirichlet on the porous side, i.e.
   !                        iface_velocity = zhang (mode 3).  This is the
   !                        literature closure for the tangential direction,
   !                        replacing the ad-hoc 'balance'/'normal' modes.
   ! The temperature side of the same section (Eqs.27-29) is switched
   ! independently by g_iface_t_model.
   integer,  save :: g_iface_zhang = 0
   !   beta, beta1 : the Eq.26 excess (viscous / inertial) stress-jump
   !   coefficients.  beta default 0 => derived as 1/bj_alpha (the standard
   !   Beavers-Joseph correspondence), so iface_slip_alpha keeps its meaning.
   real(dp), save :: g_iface_beta  = 0.0_dp
   real(dp), save :: g_iface_beta1 = 0.0_dp
   !   d_f/d_p used by Eqs.25/26 when the peer-side gap is not part of the
   !   exchange (default 1 = use the slab-side gap for both).
   real(dp), save :: g_iface_df_ratio = 1.0_dp
   ! Interface TEMPERATURE model (phase 14, opt-in; 0 = off):
   !   0 = 'off'   : legacy -- the peer temperature is imposed as a single
   !                 Dirichlet value on the porous side and the porous cell
   !                 temperature is handed back.
   !   1 = 'zhang' : Zhang 2011 Eqs.27-29 -- the clear-fluid side sees the
   !                 porosity-weighted VOLUME AVERAGE eps*<T_f>^f +
   !                 (1-eps)*<T_s>^s (Eq.27; identical to the plain value under
   !                 the LTE model), and the interface heat flux is split
   !                 between the porous phases by area ratio (porosity):
   !                 the fluid-phase equation receives eps*F (Eq.28), the
   !                 solid-phase equation (1-eps)*F (Eq.29).
   integer,  save :: g_iface_t_model = 0
   ! phase-10 auto-save / joint restart
   integer,  save :: g_save_interval = 0   ! 0 = off, >0 = save every N coupling iters
   integer,  save :: g_couple_restart = 0  ! 0 = cold start, 1 = joint restart

contains

   !---------------------------------------------------------------------------
   ! Initialise the global reference state to default sea-level values.
   ! Derived quantities (a_ref, p_ref, p_scale) are computed here.
   !---------------------------------------------------------------------------
   subroutine init_reference_state
      g_ref%rho_ref = RHO_REF
      g_ref%T_ref   = T_REF
      g_ref%L_ref   = L_REF
      g_ref%u_ref   = 0.0_dp
      call compute_derived( g_ref )
   end subroutine init_reference_state

   !---------------------------------------------------------------------------
   ! Read an optional mix.control file (key=value, '#' comments) and override
   ! the global reference state.  Missing keys keep their default value; an
   ! absent file is not an error (defaults are used).
   !
   ! Recognised keys (phase-5 minimal set; coupling iter params added phase 6;
   ! phase-13 start-up step schedule added at the end):
   !   rho_ref = <real>   [kg/m^3]
   !   T_ref   = <real>   [K]
   !   L_ref   = <real>   [m]
   !   u_ref   = <real>   [m/s]  (optional; 0 => sonic reference)
   !   n_couple            = <int>
   !   n_uns_steps         = <int>   steady SIMPLE steps per coupling iteration
   !   n_struct_steps      = <int>   steady struct substeps per coupling iteration
   !   n_struct_steps_start= <int>   cold-start struct substeps; halved every
   !                                 step_decay_every coupling iterations down
   !                                 to n_struct_steps (0 = constant cadence)
   !   n_uns_steps_start   = <int>   same, unstructured side (0 = off)
   !   step_decay_every    = <int>   halving period in coupling iterations [1]
   !   iface_relax, iface_ramp, save_interval, couple_restart
   !   iface_velocity  = full|normal|balance|zhang
   !   iface_p_anchor  = 0|1
   !   iface_p_model   = grad0|dirichlet|momentum   (phase 14)
   !   iface_slip      = off|bj|noslip              (phase 14, mod_iface_law)
   !   iface_t_model   = off|zhang                  (phase 14, Zhang Eq.27-29)
   !   iface_pm_subiter, iface_p_blend, iface_slip_alpha
   !   iface_beta, iface_beta1, iface_df_ratio      (phase 14, Zhang Eq.25/26)
   !---------------------------------------------------------------------------
   subroutine read_mix_control( filename, ier )
      character(len=*), intent(in)  :: filename
      integer,          intent(out) :: ier

      character(len=512) :: line, key, val
      integer :: u, ios, ie
      real(dp) :: rv

      ier = 0
      open( newunit = u, file = trim(filename), status = 'old', &
            action = 'read', iostat = ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'mix.control not found, using default reference state: ' &
                        // trim(filename)
         call init_reference_state()
         return
      end if

      call init_reference_state()

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

         select case ( trim(key) )
         case ( 'rho_ref' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_ref%rho_ref = rv
         case ( 'T_ref', 't_ref' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_ref%T_ref = rv
         case ( 'L_ref', 'l_ref' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_ref%L_ref = rv
         case ( 'u_ref' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_ref%u_ref = rv
         case ( 'n_couple' )
            read( val, *, iostat = ios ) g_n_couple
         case ( 'n_uns_steps' )
            read( val, *, iostat = ios ) g_n_uns_steps
         case ( 'n_struct_steps' )
            read( val, *, iostat = ios ) g_n_struct_steps
         case ( 'n_struct_steps_start' )
            read( val, *, iostat = ios ) g_n_struct_steps_start
         case ( 'n_uns_steps_start' )
            read( val, *, iostat = ios ) g_n_uns_steps_start
         case ( 'step_decay_every' )
            read( val, *, iostat = ios ) g_step_decay_every
            if ( g_step_decay_every < 1 ) g_step_decay_every = 1
         case ( 'iface_velocity', 'iface_vel_mode' )
            ! 'full' (default) | 'normal' | 'balance' | 'zhang'
            select case ( trim(lowercase(adjustl(val))) )
            case ( 'normal', 'normal-only', 'normal_only', '1' )
               g_iface_vel_mode = 1
            case ( 'balance', 'outflow', '2' )
               g_iface_vel_mode = 2
            case ( 'zhang', 'zhang2011', 'eq25', 'law25', 'zhang-l-25', '3' )
               g_iface_vel_mode = 3
            case ( 'full', '0' )
               g_iface_vel_mode = 0
            case default
               write(*,'(a)') 'WARNING: unknown iface_velocity value: ' // trim(val)
            end select
         case ( 'iface_p_anchor' )
            ! 0 = off (default), 1 = on (or on/off/yes/no)
            select case ( trim(lowercase(adjustl(val))) )
            case ( '1', 'on', 'yes', 'true', 't' )
               g_iface_p_anchor = 1
            case ( '0', 'off', 'no', 'false', 'f' )
               g_iface_p_anchor = 0
            case default
               read( val, *, iostat = ios ) g_iface_p_anchor
            end select
         case ( 'iface_p_model', 'iface_pmode' )
            ! 'grad0' (default) | 'dirichlet' | 'momentum'
            select case ( trim(lowercase(adjustl(val))) )
            case ( 'grad0', 'zero-gradient', 'zerograd', 'none', 'off', '0' )
               g_iface_p_model = 0
            case ( 'dirichlet', 'p-dirichlet', 'pdir', '1' )
               g_iface_p_model = 1
            case ( 'momentum', 'betchen', 'eq44', '2' )
               g_iface_p_model = 2
            case default
               read( val, *, iostat = ios ) g_iface_p_model
               if ( ios /= 0 .or. g_iface_p_model < 0 .or. &
                    g_iface_p_model > 2 ) then
                  write(*,'(a)') 'WARNING: unknown iface_p_model value: ' &
                                 // trim(val) // ' (using grad0)'
                  g_iface_p_model = 0
               end if
            end select
         case ( 'iface_slip', 'iface_slip_model' )
            ! 'off' (default) | 'bj' (stress jump, mod_iface_law)
            select case ( trim(lowercase(adjustl(val))) )
            case ( 'off', 'none', 'no', '0' )
               g_iface_slip = 0
            case ( 'bj', 'beavers-joseph', 'beavers_joseph', 'stressjump', &
                   'stress-jump', 'otw', '1' )
               g_iface_slip = 1
            case ( 'noslip', 'no-slip', 'cont', 'continuity', '2' )
               g_iface_slip = 2
            case default
               read( val, *, iostat = ios ) g_iface_slip
               if ( ios /= 0 .or. g_iface_slip < 1 .or. &
                    g_iface_slip > 2 ) then
                  write(*,'(a)') 'WARNING: unknown iface_slip value: ' &
                                 // trim(val) // ' (using off)'
                  g_iface_slip = 0
               end if
            end select
         case ( 'iface_t_model', 'iface_tmodel', 'iface_t_law' )
            ! 'off' (default) | 'zhang' (Eqs.27-29 volume average + flux split)
            select case ( trim(lowercase(adjustl(val))) )
            case ( 'off', 'none', 'no', 'legacy', '0' )
               g_iface_t_model = 0
            case ( 'zhang', 'zhang2011', 'eq27', 'volavg', 'volume-average', '1' )
               g_iface_t_model = 1
            case default
               read( val, *, iostat = ios ) g_iface_t_model
               if ( ios /= 0 .or. g_iface_t_model < 0 .or. &
                    g_iface_t_model > 1 ) then
                  write(*,'(a)') 'WARNING: unknown iface_t_model value: ' &
                                 // trim(val) // ' (using off)'
                  g_iface_t_model = 0
               end if
            end select
         case ( 'iface_beta' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_iface_beta = max( 0.0_dp, rv )
         case ( 'iface_beta1' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_iface_beta1 = max( 0.0_dp, rv )
         case ( 'iface_df_ratio', 'iface_d_f_ratio', 'iface_gap_ratio' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 .and. rv > 0.0_dp ) g_iface_df_ratio = rv
         case ( 'iface_pm_subiter' )
            read( val, *, iostat = ios ) g_iface_pm_subiter
            if ( ios /= 0 ) g_iface_pm_subiter = 0
            if ( g_iface_pm_subiter < 0 ) g_iface_pm_subiter = 0
         case ( 'iface_p_blend' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_iface_p_blend = max( 0.0_dp, min( 1.0_dp, rv ) )
         case ( 'iface_slip_alpha', 'iface_alpha' )
            read( val, *, iostat = ios ) rv
            if ( ios == 0 ) g_iface_slip_alpha = rv
         case ( 'iface_relax' )
            read( val, *, iostat = ios ) g_iface_relax
         case ( 'iface_ramp' )
            read( val, *, iostat = ios ) g_iface_ramp
         case ( 'save_interval' )
            read( val, *, iostat = ios ) g_save_interval
         case ( 'couple_restart' )
            read( val, *, iostat = ios ) g_couple_restart
         case default
            ! unknown keys are silently ignored (forward compatibility)
         end select
      end do
      close( u )

      call compute_derived( g_ref )

      write(*,'(a)') ''
      write(*,'(a)') '--- reference state (SI) ---'
      write(*,'(a,es12.4,a)') '  rho_ref = ', g_ref%rho_ref, ' [kg/m^3]'
      write(*,'(a,es12.4,a)') '  T_ref   = ', g_ref%T_ref,   ' [K]'
      write(*,'(a,es12.4,a)') '  L_ref   = ', g_ref%L_ref,   ' [m]'
      write(*,'(a,es12.4,a)') '  a_ref   = ', g_ref%a_ref,   ' [m/s]'
      write(*,'(a,es12.4,a)') '  u_ref   = ', g_ref%u_ref,   ' [m/s] (struct velocity scale U_inf)'
      write(*,'(a,es12.4,a)') '  p_ref   = ', g_ref%p_ref,   ' [Pa]'
      write(*,'(a,es12.4,a)') '  p_scale = ', g_ref%p_scale, ' [Pa] (rho_ref*u_ref^2)'
      write(*,'(a,i0)')       '  save_interval = ', g_save_interval
      write(*,'(a,i0)')       '  couple_restart = ', g_couple_restart
      write(*,'(a,i0)')       '  iface_vel_mode = ', g_iface_vel_mode
      write(*,'(a,i0)')       '  iface_p_anchor = ', g_iface_p_anchor
      write(*,'(a,i0,a)')     '  iface_p_model  = ', g_iface_p_model, &
           merge(' (dirichlet)', '            ', g_iface_p_model == 1)
      if ( g_iface_p_model == 2 ) &
         write(*,'(a,i0,a,f5.2)') '    momentum: pm_subiter = ', &
              g_iface_pm_subiter, '  p_blend = ', g_iface_p_blend
      write(*,'(a,i0)')       '  iface_slip     = ', g_iface_slip
      if ( g_iface_slip == 1 ) &
         write(*,'(a,f8.4)')  '    bj alpha      = ', g_iface_slip_alpha
      write(*,'(a,i0,a)')     '  iface_t_model  = ', g_iface_t_model, &
           merge(' (zhang Eq.27-29)', '                 ', g_iface_t_model == 1)
      if ( g_iface_vel_mode == 3 .or. g_iface_zhang == 1 ) &
         write(*,'(a,2f10.4,a,f8.4)') '    zhang beta/beta1 = ', &
              g_iface_beta, g_iface_beta1, '  d_f/d_p =', g_iface_df_ratio
      write(*,'(a)') '--- end reference state ---'
   end subroutine read_mix_control

   !---------------------------------------------------------------------------
   ! Accessor for the global reference state.
   !---------------------------------------------------------------------------
   function get_ref_state() result( st )
      type(reference_state_t) :: st
      st = g_ref
   end function get_ref_state

   !---------------------------------------------------------------------------
   ! Accessor for the coupling-loop parameters.
   !---------------------------------------------------------------------------
   subroutine get_coupling_params( n_couple, n_uns_steps, iface_relax, iface_ramp, &
                                   n_struct_steps, save_interval, couple_restart, &
                                   n_struct_steps_start, n_uns_steps_start, &
                                   step_decay_every, iface_vel_mode, iface_p_anchor, &
                                   iface_p_model, iface_slip, iface_pm_subiter, &
                                   iface_p_blend, iface_slip_alpha, &
                                   iface_t_model, iface_beta, iface_beta1, &
                                   iface_df_ratio )
      integer,  intent(out) :: n_couple, n_uns_steps, iface_ramp
      integer,  intent(out), optional :: n_struct_steps
      integer,  intent(out), optional :: save_interval, couple_restart
      integer,  intent(out), optional :: n_struct_steps_start, n_uns_steps_start
      integer,  intent(out), optional :: step_decay_every
      integer,  intent(out), optional :: iface_vel_mode, iface_p_anchor
      integer,  intent(out), optional :: iface_p_model, iface_slip
      integer,  intent(out), optional :: iface_pm_subiter
      real(dp), intent(out), optional :: iface_p_blend, iface_slip_alpha
      integer,  intent(out), optional :: iface_t_model
      real(dp), intent(out), optional :: iface_beta, iface_beta1, iface_df_ratio
      real(dp), intent(out) :: iface_relax
      n_couple    = g_n_couple
      n_uns_steps = g_n_uns_steps
      iface_relax = g_iface_relax
      iface_ramp  = g_iface_ramp
      if ( present(n_struct_steps) ) n_struct_steps = g_n_struct_steps
      if ( present(save_interval) )  save_interval  = g_save_interval
      if ( present(couple_restart) ) couple_restart = g_couple_restart
      if ( present(n_struct_steps_start) ) n_struct_steps_start = g_n_struct_steps_start
      if ( present(n_uns_steps_start) )    n_uns_steps_start    = g_n_uns_steps_start
      if ( present(step_decay_every) )     step_decay_every     = g_step_decay_every
      if ( present(iface_vel_mode) )       iface_vel_mode       = g_iface_vel_mode
      if ( present(iface_p_anchor) )       iface_p_anchor       = g_iface_p_anchor
      if ( present(iface_p_model) )        iface_p_model        = g_iface_p_model
      if ( present(iface_slip) )           iface_slip           = g_iface_slip
      if ( present(iface_pm_subiter) )     iface_pm_subiter     = g_iface_pm_subiter
      if ( present(iface_p_blend) )        iface_p_blend        = g_iface_p_blend
      if ( present(iface_slip_alpha) )     iface_slip_alpha     = g_iface_slip_alpha
      if ( present(iface_t_model) )        iface_t_model        = g_iface_t_model
      if ( present(iface_beta) )           iface_beta           = g_iface_beta
      if ( present(iface_beta1) )          iface_beta1          = g_iface_beta1
      if ( present(iface_df_ratio) )       iface_df_ratio       = g_iface_df_ratio
   end subroutine get_coupling_params

   !---------------------------------------------------------------------------
   ! scheduled_steps -- number of internal solver steps a side runs during
   ! coupling iteration `iter` under the optional start-up (halving) schedule.
   !
   !   n(iter) = max( n_target, n_start / 2**((iter-1)/decay_every) )
   !
   ! Rationale: a large step-per-exchange ratio at cold start lets each side
   ! settle onto its own manifold before the interface is updated, then the
   ! ratio is halved every `decay_every` coupling iterations until both sides
   ! exchange at the same cadence (n_start == n_target => constant cadence).
   ! The halving is done by repeated integer division so no power of two is
   ! ever formed (no overflow for a large iter).  n_start <= 0 (or
   ! n_start <= n_target) disables the schedule.
   !---------------------------------------------------------------------------
   pure function scheduled_steps( iter, n_start, n_target, decay_every ) result( n )
      integer, intent(in) :: iter, n_start, n_target
      integer, intent(in), optional :: decay_every
      integer :: n, k, every, nhalve

      n = max( n_target, 0 )
      if ( n_start <= 0 .or. n_start <= n_target ) return

      every = 1
      if ( present(decay_every) ) every = max( 1, decay_every )

      nhalve = ( max(iter, 1) - 1 ) / every
      n = n_start
      do k = 1, nhalve
         if ( n <= n_target ) exit
         n = max( n_target, n / 2 )
      end do
   end function scheduled_steps

   !---------------------------------------------------------------------------
   ! Compute derived reference quantities from the independent ones.
   !---------------------------------------------------------------------------
   subroutine compute_derived( st )
      type(reference_state_t), intent(inout) :: st
      st%a_ref   = sqrt( GAMMA * R_GAS * st%T_ref )
      st%p_ref   = st%rho_ref * R_GAS * st%T_ref
      if ( st%u_ref <= 0.0_dp ) st%u_ref = st%a_ref
      ! OpenCFD-EC scales pressure by the dynamic pressure rho_ref*U_inf^2,
      ! not by rho_ref*a_ref^2 (u* = u/U_inf, p_inf* = 1/(gamma*Ma^2)).
      st%p_scale = st%rho_ref * st%u_ref * st%u_ref
   end subroutine compute_derived

   !---------------------------------------------------------------------------
   ! local lowercase helper (avoids depending on solver-specific modules)
   !---------------------------------------------------------------------------
   function lowercase( s ) result( r )
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, c
      r = s
      do i = 1, len_trim(s)
         c = iachar( r(i:i) )
         if ( c >= iachar('A') .and. c <= iachar('Z') ) &
            r(i:i) = achar( c + 32 )
      end do
   end function lowercase

end module mod_reference_state
