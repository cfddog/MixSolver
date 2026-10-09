!===============================================================================
! mod_iface_law.f90 -- fluid/porous interface closure laws (phase 14).
!
! Pure, dependency-free (mod_precision only) closures shared by the coupling
! driver and the unstructured solver:
!
!   1) tangential (shear) stress jump -- Beavers-Joseph / Ochoa-Tapia-Whitaker
!      (Ochoa-Tapia & Whitaker 1995; Zhang 2011 Eq.24-26; Betchen 2006 Eq.35).
!   2) interface pressure from the normal momentum balance across the abrupt
!      flow-area change (Betchen 2006 Eq.43/44) combined with the two-sided
!      inverse-distance estimate (Betchen 2006 Sec.4.2 / Eq.39-41).
!
! The stress-jump law is the SAME law the unstructured solver already applies
! to *internal* fluid/porous faces (the bj_alpha branch of
! mod_uns_simple:momentum_assembly, validated against the two-layer ODE
! reference in cases/beavers_joseph): this module owns the formula so that the
! internal face and the cross-solver coupling interface use one closure.
! With d_p = 0 (single-domain face) `iface_slip_conductance` reduces to the
! internal-face coefficient  C = mu*A/(d_Pf + lam/alpha).
!
! Everything here is pure (no state, no I/O) so it can be unit-tested by
! src/coupling/test_iface_law.f90 and called from either side of the partition.
!===============================================================================
module mod_iface_law
   use mod_precision, only: dp
   implicit none
   private
   public :: iface_slip_conductance, iface_slip_velocity, &
             iface_p_momentum, iface_p_blend, iface_p_flux_response, &
             iface_p_gain, &
             iface_zhang_vn, iface_zhang_vt, iface_zhang_velocity, &
             iface_zhang_T, iface_volavg_T

contains

   !----------------------------------------------------------------------------
   ! iface_slip_conductance -- series-resistance tangential conductance C [kg/s]
   ! of an interface patch of area `area`:
   !
   !    1/C = (d_f + lam/alpha)/mu + eps*d_p/mu_e ,     mu_e = mu/eps
   !
   ! so that the tangential momentum flux (force) across the interface is
   !
   !    F_t = C * ( u_t^fluid - u_t^porous )                [N]
   !
   ! Terms:
   !   d_f/mu        fluid-side one-sided viscous resistance (peer cell -> face)
   !   eps*d_p/mu    porous-side (Brinkman) resistance, eps-weighted exactly as
   !                 the harmonic-mean face coefficient of Betchen 2006 Eq.41
   !                 (derived from Eqs.37-40: D = mu*A/(d_f + eps*d_p))
   !   lam/alpha     Beavers-Joseph stress-jump (slip) resistance with
   !                 lam = sqrt(K), alpha = BJ slip-resistance coefficient
   !                 (alpha -> inf : no slip, i.e. plain stress continuity;
   !                  alpha > 0 small : slippery, C -> 0 = free slip)
   !
   ! Limits used by the unit test:
   !   d_p = 0            -> C = mu*A/(d_f + lam/alpha): the internal-face
   !                         bj_alpha branch of mod_uns_simple (whose continuous
   !                         limit 1/C = d_Pf/(mu*A) + lam/(mu*alpha*A) is the
   !                         two-layer reference of cases/beavers_joseph)
   !   alpha = 0 (or < 0) -> C = mu*A/(d_f + eps*d_p): stress continuity
   !   lam = 0 (K -> 0)   -> C = mu*A/(d_f + eps*d_p): impermeable limit
   !----------------------------------------------------------------------------
   pure function iface_slip_conductance( area, d_f, d_p, lam, mu, alpha, eps ) &
                                         result( C )
      real(dp), intent(in)  :: area, d_f, d_p, lam, mu, alpha, eps
      real(dp)              :: C, R
      real(dp), parameter   :: r_floor = 1.0e-30_dp

      R = d_f + eps * d_p
      if ( alpha > 0.0_dp .and. lam > 0.0_dp ) R = R + lam / alpha
      if ( R < r_floor ) R = r_floor        ! zero gap / zero slip: C -> huge
      if ( mu <= 0.0_dp .or. area <= 0.0_dp ) then
         C = 0.0_dp
      else
         C = mu * area / R
      end if
   end function iface_slip_conductance

   !----------------------------------------------------------------------------
   ! iface_slip_velocity -- the value a *Dirichlet* interface face must take to
   ! reproduce the stress-jump flux through a face whose flush ("no-slip")
   ! Dirichlet conductance is D_por = mu_e*A/d_p = mu*A/(eps*d_p):
   !
   !    D_por*( V_imp - u_t^por ) = C*( u_t^fl - u_t^por )
   ! => V_imp = u_t^por + (C/D_por) * ( u_t^fl - u_t^por )
   !
   ! C/D_por = eps*d_p/(d_f + lam/alpha + eps*d_p) is always in [0,1], so V_imp
   ! always lies between the porous cell value and the peer value; the limits are
   !   C -> 0       : V_imp = u_t^por  (zero-gradient = free slip / no grip)
   !   C -> D_por   : V_imp = u_t^fl   (no slip)
   ! For faces whose normal is aligned with a coordinate axis this Dirichlet
   ! value reproduces the Robin flux of the stress-jump condition exactly (the
   ! identity used by the coupling driver's iface_slip path).
   !----------------------------------------------------------------------------
   pure function iface_slip_velocity( C, D_por, u_t_fluid, u_t_porous ) result( V )
      real(dp), intent(in) :: C, D_por, u_t_fluid, u_t_porous
      real(dp) :: V, ratio

      if ( D_por <= 0.0_dp ) then
         V = u_t_porous
         return
      end if
      ratio = C / D_por                 ! <= 1 analytically; guard round-off
      ratio = max( 0.0_dp, min( 1.0_dp, ratio ) )
      V = u_t_porous + ratio * ( u_t_fluid - u_t_porous )
   end function iface_slip_velocity

   !----------------------------------------------------------------------------
   ! iface_p_momentum -- interface pressure estimate from the normal momentum
   ! balance across the flow-area change (Betchen 2006 Eq.43/44):
   !
   !    p_i = p_fluid - ( mdot * (u_i . n) / A ) * (1 - eps)/eps
   !
   ! n is the interface normal pointing OUT of the porous region, so that
   ! mdot = rho*A*(u.n) > 0 for coolant blowing from the porous slab into the
   ! mainstream; p_fluid is the *fluid-side* pressure extrapolated to the
   ! interface (the peer's face value).  Physical content (Betchen Sec.3.1/4.2):
   ! only a fraction eps of the interface normal stress is carried by the fluid
   ! constituent of the porous medium, so the remaining (1-eps)/eps fraction of
   ! the dynamic head rho*u_n^2 appears as a pressure drop on the porous side of
   ! the nominal interface (the pressure GRADIENT may jump even though p itself
   ! is continuous).
   !----------------------------------------------------------------------------
   pure function iface_p_momentum( p_fluid, mdot, u_n, area, eps ) result( p )
      real(dp), intent(in) :: p_fluid, mdot, u_n, area, eps
      real(dp) :: p, e

      e = eps
      if ( e < 1.0e-6_dp ) e = 1.0e-6_dp     ! guard eps -> 0 (solid)
      if ( e > 1.0_dp    ) e = 1.0_dp

      p = p_fluid
      if ( area > 0.0_dp ) p = p - ( mdot * u_n / area ) * ( 1.0_dp - e ) / e
   end function iface_p_momentum

   !----------------------------------------------------------------------------
   ! iface_p_blend -- two-sided interface pressure estimate.  Betchen 2006
   ! Sec.4.2: the fluid-side estimate of Eq.44 is *averaged* with an estimate
   ! extrapolated from the porous region; inverse-distance weighting is the
   ! "simpler estimate" implemented implicitly, the two-sided / momentum
   ! correction being deferred.
   !    p = (1-w)*p_fluid + w*p_porous,  w = weight of the porous-side estimate.
   ! w = 0.5 is the inverse-distance average for d_f = d_p.
   !----------------------------------------------------------------------------
   pure function iface_p_blend( p_fluid, p_porous, w ) result( p )
      real(dp), intent(in) :: p_fluid, p_porous, w
      real(dp) :: p, wc

      wc = max( 0.0_dp, min( 1.0_dp, w ) )
      p  = ( 1.0_dp - wc ) * p_fluid + wc * p_porous
   end function iface_p_blend

   !----------------------------------------------------------------------------
   ! iface_p_flux_response -- one step of the deferred p-mdot sub-iteration.
   ! The interface mass flux implied by a *changed* Dirichlet interface pressure
   ! follows from the PPE's own linear flux response
   !
   !    d(mdot)/d(p_face) = -af = -rho * (V/a_P) * A/dn
   !
   ! (same af as the Dirichlet-face treatment in
   ! mod_uns_simple:ppe_assembly / correct_fields).  Betchen 2006 Sec.4.2 notes
   ! that the interface pressure depends on the interface mass flow rate, which
   ! in turn depends on the interface pressure, so "a small number of iterations
   ! is required to ensure accuracy" -- this is that iteration, run in a deferred
   ! fashion (cell pressure and u_f frozen over the sub-loop).
   !----------------------------------------------------------------------------
   pure function iface_p_flux_response( mdot, af, p_new, p_old ) result( mdot_new )
      real(dp), intent(in) :: mdot, af, p_new, p_old
      real(dp) :: mdot_new

      mdot_new = mdot - af * ( p_new - p_old )
   end function iface_p_flux_response

   !----------------------------------------------------------------------------
   ! iface_p_gain -- dimensionless local loop gain of the p-mdot coupling
   !
   !     G = |d(p_face)/d(mdot)| * |d(mdot)/d(p_face)|
   !       = [ 2*(1-eps)*|mdot| / (eps*rho*A^2) ] * af
   !       = 2*(1-eps)*af*|mdot| / (eps*rho*A^2)
   !
   ! The first factor is the Eq.43/44 normal-momentum sensitivity (the pressure
   ! scales as rho*u_n^2), the second is the PPE's Dirichlet-face flux response
   ! (af = rho*(V/a_P)*A/dn, the same af as ppe_assembly / correct_fields).
   ! G is what decides how many sub-iterations the p-mdot pair tolerates:
   !   G <~ 1 : the deferred (explicit) estimate is already the fixed point and
   !            one evaluation suffices;
   !   G >> 1 : the undamped fixed point iteration DIVERGES (it oscillates with
   !            growth G per step), so the sub-iteration must be under-relaxed
   !            with omega <= 1/(1+G) -- or simply left deferred, in which case
   !            the outer coupling loop closes the p-mdot loop.
   ! Returned as the raw gain (not clamped) so callers can report it.
   !----------------------------------------------------------------------------
   pure function iface_p_gain( af, mdot, eps, rho, area ) result( G )
      real(dp), intent(in) :: af, mdot, eps, rho, area
      real(dp) :: G, e

      e = eps
      if ( e < 1.0e-6_dp ) e = 1.0e-6_dp
      if ( e > 1.0_dp    ) e = 1.0_dp
      G = 0.0_dp
      if ( rho > 0.0_dp .and. area > 0.0_dp ) &
         G = 2.0_dp * ( 1.0_dp - e ) * af * abs( mdot ) / ( e * rho * area**2 )
   end function iface_p_gain

   !===========================================================================
   ! Zhang 2011 (Sec.3.5.1) porous/clear-fluid interface closure.
   !
   ! Eq.22 : V|_fl = <V>|_p = V|_interface -- ONE interface velocity, shared
   !         by the two domains (mass balance).
   ! Eq.23/25 (normal stress):
   !         (mu_e/eps)(V_n - <V>_n^p)/d_p = mu_f (V_n^fl - V_n)/d_f
   ! Eq.24/26 (tangential stress jump, Ochoa-Tapia & Whitaker):
   !         (mu_e/eps)(V_t - <V>_t^p)/d_p - mu_f (V_t^fl - V_t)/d_f
   !           = beta (mu_e/sqrt(K)) V_t + beta1 rho_f |V_t| V_t
   !
   ! <V>^p and V^fl are the cell-centre velocity vectors of the porous and
   ! clear-fluid domains nearest the interface, at distances d_p and d_f from
   ! it.  With the two-sided viscous conductances
   !
   !    G_p = mu_e/(eps*d_p) ,   G_f = mu_f/d_f ,   Q = G_p U^p + G_f U^fl
   !
   ! Eq.25 is the harmonic (conductance-weighted) blend
   !
   !    V_n = Q_n / (G_p + G_f)
   !
   ! and Eq.26 is a Robin balance with an *excess* stress-jump resistance
   ! beta*mu_e/sqrt(K) plus an inertial (Forchheimer-like) term:
   !
   !    B = G_p + G_f + beta*mu_e/sqrt(K)
   !    beta1*rho*x^2 + B*x - |Q_t| = 0 ,  x = |V_t|
   !    V_t = (2*|Q_t| / (B + sqrt(B^2 + 4*beta1*rho*|Q_t|))) * Q_t/|Q_t|
   !
   ! (the numerically stable root; beta1 = 0 -> V_t = Q_t/B).  Parameter
   ! correspondence with the Beavers-Joseph numbering used elsewhere in this
   ! code (mod_uns_bc:bj_alpha, iface_slip_alpha): beta = 1/alpha, so
   ! alpha -> 0 (beta -> inf) drives V_t -> 0 = no slip, while alpha -> inf
   ! (beta -> 0) recovers the pure harmonic blend = stress continuity.
   !
   ! mu_e is the porous-domain effective (Brinkman) viscosity; passing
   ! mu_e = mu_f gives the "intrinsic viscosity" reading, passing
   ! mu_e = mu_f/eps reproduces the internal-face bj_alpha branch of
   ! mod_uns_simple (see iface_slip_conductance above).
   !===========================================================================

   !----------------------------------------------------------------------------
   ! iface_zhang_vn -- Eq.25: normal interface velocity [m/s].
   !----------------------------------------------------------------------------
   pure function iface_zhang_vn( un_p, un_f, d_p, d_f, mu, mu_e, eps ) &
                                          result( Vn )
      real(dp), intent(in) :: un_p, un_f, d_p, d_f, mu, mu_e, eps
      real(dp) :: Vn, Gp, Gf, e

      e  = max( eps, 1.0e-12_dp )
      Gp = 0.0_dp
      Gf = 0.0_dp
      if ( d_p > 0.0_dp .and. mu_e > 0.0_dp ) Gp = mu_e / ( e * d_p )
      if ( d_f > 0.0_dp .and. mu   > 0.0_dp ) Gf = mu   / d_f

      if ( Gp + Gf > 0.0_dp ) then
         Vn = ( Gp * un_p + Gf * un_f ) / ( Gp + Gf )
      else
         Vn = un_p                     ! no two-sided information
      end if
   end function iface_zhang_vn

   !----------------------------------------------------------------------------
   ! iface_zhang_vt -- Eq.26: *scalar* (one tangential direction) interface
   ! velocity [m/s].  Called by iface_zhang_velocity with the two cell values
   ! projected on the slip direction, which turns the vector law into exactly
   ! this scalar balance with ut_p/ut_f = the projections and the return value
   ! = the tangential speed.
   !   beta  = excess viscous stress-jump coefficient (beta = 1/alpha_BJ)
   !   beta1 = excess inertial coefficient (default 0 = linear Beavers-Joseph)
   !   lam   = sqrt(K) (interface-normal permeability length)
   !----------------------------------------------------------------------------
   pure function iface_zhang_vt( ut_p, ut_f, d_p, d_f, mu, mu_e, eps, lam, &
                                 beta, beta1, rho ) result( Vt )
      real(dp), intent(in) :: ut_p, ut_f, d_p, d_f, mu, mu_e, eps, lam, &
                              beta, beta1, rho
      real(dp) :: Vt, Gp, Gf, B, Q, disc, b2, e

      e  = max( eps, 1.0e-12_dp )
      Gp = 0.0_dp
      Gf = 0.0_dp
      if ( d_p > 0.0_dp .and. mu_e > 0.0_dp ) Gp = mu_e / ( e * d_p )
      if ( d_f > 0.0_dp .and. mu   > 0.0_dp ) Gf = mu   / d_f

      b2 = 0.0_dp
      if ( lam > 0.0_dp .and. beta > 0.0_dp ) b2 = beta * mu_e / lam
      B  = Gp + Gf + b2
      Q  = Gp * ut_p + Gf * ut_f

      if ( B <= 0.0_dp ) then
         Vt = ut_p
      else if ( beta1 * rho <= 0.0_dp ) then
         Vt = Q / B
      else
         disc = B * B + 4.0_dp * beta1 * rho * Q
         if ( disc > 0.0_dp ) then
            Vt = 2.0_dp * Q / ( B + sqrt( disc ) )
         else
            ! no real root (the inertial excess cannot be balanced by the
            ! available conductances): drop it and keep the linear root
            Vt = Q / B
         end if
      end if
   end function iface_zhang_vt

   !----------------------------------------------------------------------------
   ! iface_zhang_velocity -- Eqs.22 + 25 + 26: the single interface velocity
   ! vector [m/s] shared by both domains.
   !   u_p   : porous-side cell-centre velocity (point P', distance d_p)
   !   u_f   : clear-fluid-side cell-centre velocity (point F', distance d_f)
   !   nvec  : interface normal (any orientation; normalised internally),
   !           pointing OUT of the porous domain
   ! The tangential part is isotropic: the scalar Eq.26 balance is applied on
   ! the slip direction Q_t/|Q_t| (Q_t = the tangential part of G_p u_p +
   ! G_f u_f), which is exact whenever the friction is aligned with V_t.
   !----------------------------------------------------------------------------
   pure function iface_zhang_velocity( u_p, u_f, nvec, d_p, d_f, mu, mu_e, &
                                       eps, lam, beta, beta1, rho ) result( V )
      real(dp), intent(in) :: u_p(3), u_f(3), nvec(3)
      real(dp), intent(in) :: d_p, d_f, mu, mu_e, eps, lam, beta, beta1, rho
      real(dp) :: V(3), n(3), t_p(3), t_f(3), Qt(3), chat(3)
      real(dp) :: un_p, un_f, Vn, Vt, s, Gp, Gf, e

      e = max( eps, 1.0e-12_dp )
      n = nvec
      s = norm2( n )
      if ( s > 0.0_dp ) n = n / s

      un_p = dot_product( u_p, n )
      un_f = dot_product( u_f, n )
      Vn   = iface_zhang_vn( un_p, un_f, d_p, d_f, mu, mu_e, e )

      t_p = u_p - un_p * n
      t_f = u_f - un_f * n

      Gp = 0.0_dp
      Gf = 0.0_dp
      if ( d_p > 0.0_dp .and. mu_e > 0.0_dp ) Gp = mu_e / ( e * d_p )
      if ( d_f > 0.0_dp .and. mu   > 0.0_dp ) Gf = mu   / d_f

      Qt = Gp * t_p + Gf * t_f
      s  = norm2( Qt )
      if ( s > 1.0e-30_dp ) then
         chat = Qt / s
         Vt   = iface_zhang_vt( dot_product( t_p, chat ), &
                                dot_product( t_f, chat ), &
                                d_p, d_f, mu, mu_e, e, lam, beta, beta1, rho )
         V    = Vn * n + Vt * chat
      else
         V    = Vn * n
      end if
   end function iface_zhang_velocity

   !===========================================================================
   ! Zhang 2011 Eqs.27-29: interface temperature / heat-flux split.
   !
   !   Eq.27 : T_fl = <T>^p = eps <T_f>^f + (1-eps) <T_s>^s
   !   Eq.28 : eps     k_f dT_fl/dn|_fl = k_fe d<T_f>^f/dn
   !   Eq.29 : (1-eps) k_f dT_fl/dn|_fl = k_se d<T_s>^s/dn
   !
   ! The interface temperature seen by the CLEAR FLUID is therefore *unique*
   ! and equal to the porosity-weighted VOLUME AVERAGE of the two porous-phase
   ! temperatures -- it is neither <T_f>^f nor <T_s>^s alone.  The total heat
   ! flux arriving from the clear fluid is split between the phases by area
   ! ratio (porosity): eps*F to the fluid phase, (1-eps)*F to the solid phase.
   !===========================================================================

   !----------------------------------------------------------------------------
   ! iface_zhang_T -- closed-form solution of the discretised Eqs.27-29.
   !   in : T_flP = clear-fluid-side temperature next to the interface (peer)
   !        T_fP, T_sP = porous-side fluid/solid cell-centre temperatures
   !        eps, k_f, k_fe, k_se, d_f (clear-fluid gap), d_p (porous gap)
   !   out: T_fl = clear-fluid INTERFACE temperature = <T>^p  (Eq.27)
   !        T_fi, T_si = porous fluid/solid INTERFACE temperatures
   !        F    = interface heat-flux density [W/m^2], + = fluid -> porous
   !
   ! Energy conservation across the interface (flux leaving the clear fluid =
   ! flux entering the porous phases, split by area ratio eps/(1-eps)):
   !        F   = k_f (T_flP - T_fl)/d_f
   !        eps*F     = k_fe (T_fi - T_fP)/d_p            (Eq.28)
   !        (1-eps)*F = k_se (T_si - T_sP)/d_p            (Eq.29)
   !        T_fl = eps*T_fi + (1-eps)*T_si                (Eq.27)
   ! Substituting the two phase relations into Eq.27 and dividing by the total
   ! interface resistance gives the closed form
   !        F   = (T_flP - [eps*T_fP + (1-eps)*T_sP]) /
   !              ( d_f/k_f + d_p*( eps^2/k_fe + (1-eps)^2/k_se ) )
   ! (a SUM of the two one-sided resistances), and then the three interface
   ! temperatures follow by back-substitution.  The returned values satisfy
   ! Eqs.27-29 identically; the phase fluxes k_fe(T_fi-T_fP)/d_p and
   ! k_se(T_si-T_sP)/d_p are exactly eps*F and (1-eps)*F, so their sum is the
   ! clear-fluid flux k_f(T_flP-T_fl)/d_f.
   !----------------------------------------------------------------------------
   pure subroutine iface_zhang_T( T_flP, T_fP, T_sP, eps, k_f, k_fe, k_se, &
                                  d_f, d_p, T_fl, T_fi, T_si, F )
      real(dp), intent(in)  :: T_flP, T_fP, T_sP, eps, k_f, k_fe, k_se, &
                               d_f, d_p
      real(dp), intent(out) :: T_fl, T_fi, T_si, F
      real(dp) :: e, den, V

      e   = max( eps, 1.0e-12_dp )
      e   = min( e, 1.0_dp - 1.0e-12_dp )
      den = 0.0_dp
      if ( k_f > 0.0_dp .and. d_f > 0.0_dp ) den = d_f / k_f
      if ( k_fe > 0.0_dp ) den = den + d_p * e * e / k_fe
      if ( k_se > 0.0_dp ) den = den + d_p * ( 1.0_dp - e )**2 / k_se
      V   = e * T_fP + ( 1.0_dp - e ) * T_sP

      ! No clear-fluid conductance (k_f or d_f vanished) means the interface
      ! cannot transmit the split flux at all: return the one-sided Dirichlet
      ! with zero flux (finite, NaN-free) instead of a runaway interface value.
      if ( den <= 1.0e-30_dp .or. k_f <= 0.0_dp .or. d_f <= 0.0_dp ) then
         F    = 0.0_dp
         T_fl = T_flP
         T_fi = T_fP
         T_si = T_sP
      else
         F    = ( T_flP - V ) / den
         T_fl = T_flP
         T_fi = T_fP
         T_si = T_sP
         if ( k_f  > 0.0_dp ) T_fl = T_flP - F * d_f / k_f
         if ( k_fe > 0.0_dp ) T_fi = T_fP + e * F * d_p / k_fe
         if ( k_se > 0.0_dp ) T_si = T_sP + ( 1.0_dp - e ) * F * d_p / k_se
      end if
   end subroutine iface_zhang_T

   !----------------------------------------------------------------------------
   ! iface_volavg_T -- Eq.27 volume average  <T>^p = eps*<T_f>^f +
   ! (1-eps)*<T_s>^s : the temperature the clear-fluid side must see.  With
   ! the LTE model (T_f = T_s = T) it reduces to T unchanged, so the same
   ! call is valid for both thermal models.
   !----------------------------------------------------------------------------
   pure function iface_volavg_T( T_f, T_s, eps ) result( T )
      real(dp), intent(in) :: T_f, T_s, eps
      real(dp) :: T, e

      e = max( 0.0_dp, min( 1.0_dp, eps ) )
      T = e * T_f + ( 1.0_dp - e ) * T_s
   end function iface_volavg_T

end module mod_iface_law
