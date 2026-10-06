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
   ! Recognised keys (phase-5 minimal set; coupling iter params added phase 6):
   !   rho_ref = <real>   [kg/m^3]
   !   T_ref   = <real>   [K]
   !   L_ref   = <real>   [m]
   !   u_ref   = <real>   [m/s]  (optional; 0 => sonic reference)
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
                                   n_struct_steps, save_interval, couple_restart )
      integer,  intent(out) :: n_couple, n_uns_steps, iface_ramp
      integer,  intent(out), optional :: n_struct_steps
      integer,  intent(out), optional :: save_interval, couple_restart
      real(dp), intent(out) :: iface_relax
      n_couple    = g_n_couple
      n_uns_steps = g_n_uns_steps
      iface_relax = g_iface_relax
      iface_ramp  = g_iface_ramp
      if ( present(n_struct_steps) ) n_struct_steps = g_n_struct_steps
      if ( present(save_interval) )  save_interval  = g_save_interval
      if ( present(couple_restart) ) couple_restart = g_couple_restart
   end subroutine get_coupling_params

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
