!===============================================================================
! test_iface_law.f90 -- phase-14 unit tests for the shared interface closures
! (mod_iface_law): the Beavers-Joseph / Ochoa-Tapia-Whitaker tangential
! stress-jump conductance, the Betchen 2006 Eq.43/44 interface pressure, the
! two-sided blend, the p-mdot flux response and its local loop gain.
!
! Run:  make iface_law_test && bin/iface_law_test
!
! The tests deliberately include the *discretizations of the two production
! cases* so that the same numbers can be followed in
! cases/beavers_joseph/README.md (internal fluid/porous faces, bj_alpha) and
! cases/../C_P_test (cross-solver coupling interface, iface_slip / iface_p_model).
!===============================================================================
program test_iface_law
   use mod_precision, only: dp
   use mod_iface_law, only: iface_slip_conductance, iface_slip_velocity, &
                            iface_p_momentum, iface_p_blend, &
                            iface_p_flux_response, iface_p_gain, &
                            iface_zhang_vn, iface_zhang_vt, &
                            iface_zhang_velocity, iface_zhang_T, &
                            iface_volavg_T
   implicit none

   integer :: nfail
   nfail = 0

   call test_bj_internal_equivalence( nfail )
   call test_bj_limits( nfail )
   call test_bj_case_discretization( nfail )
   call test_slip_velocity_identity( nfail )
   call test_p_momentum( nfail )
   call test_p_blend( nfail )
   call test_p_mdot_pair( nfail )
   call test_zhang_vn( nfail )
   call test_zhang_vt( nfail )
   call test_zhang_velocity_vector( nfail )
   call test_zhang_T( nfail )

   write(*,'(a,i0,a,i0,a)') '=== 11 tests, ', nfail, ' failures ==='
   if ( nfail > 0 ) then
      write(*,'(a)') 'FAIL'
      stop 1
   end if
   write(*,'(a)') 'PASS'

contains

   !----------------------------------------------------------------------------
   ! (1) The coupling-interface law with d_p = 0 must reproduce the *internal*
   ! face coefficient used by mod_uns_simple's bj_alpha branch,
   !     C_int = mu*A*(alpha/lam) / (1 + alpha*d_f/lam)
   ! (that branch's continuous limit is the two-layer reference solution of
   ! cases/beavers_joseph).  One law, two applications.
   !----------------------------------------------------------------------------
   subroutine test_bj_internal_equivalence( nfail )
      integer, intent(inout) :: nfail
      integer, parameter :: NC = 5
      integer :: k, bad
      real(dp) :: A, mu, alpha(NC), df(NC), lam(NC), C1, C2

      A  = 1.0e-4_dp
      mu = 1.846e-5_dp
      df   = (/ 1.25e-4_dp, 2.0e-4_dp, 1.0e-3_dp, 5.0e-5_dp, 2.5e-4_dp /)
      lam  = (/ 1.0e-3_dp,  7.071e-4_dp, 1.0e-3_dp, 3.0e-4_dp, 1.0e-3_dp /)
      alpha= (/ 1.0_dp, 0.5_dp, 2.0_dp, 10.0_dp, 0.1_dp /)

      bad = 0
      do k = 1, NC
         C1 = iface_slip_conductance( A, df(k), 0.0_dp, lam(k), mu, alpha(k), &
                                      1.0_dp )
         C2 = mu * A * ( alpha(k) / lam(k) ) / &
              ( 1.0_dp + alpha(k) * df(k) / lam(k) )
         if ( abs( C1 - C2 ) > 1.0e-14_dp * abs(C2) ) bad = bad + 1
      end do
      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0,a)') 'FAIL (1): coupling law != internal bj_alpha law for ', &
                             bad, ' cases'
      else
         write(*,'(a)') 'PASS (1): d_p=0 reduces to the internal bj_alpha law'
      end if
   end subroutine test_bj_internal_equivalence

   !----------------------------------------------------------------------------
   ! (2) Limits of the stress-jump conductance:
   !   alpha -> 0   : free slip, C -> 0
   !   alpha -> inf : no slip, C -> mu*A/(d_f + eps*d_p) (stress continuity)
   !   eps*d_p      : the porous-side Brinkman resistance, eps-weighted
   !----------------------------------------------------------------------------
   subroutine test_bj_limits( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp), parameter :: A = 1.0e-4_dp, mu = 1.846e-5_dp
      real(dp), parameter :: df = 2.5e-4_dp, dgap = 2.5e-4_dp
      real(dp), parameter :: lam = 1.0e-3_dp
      real(dp) :: C_free, C_noslip, C_ref, C_ref2, C_half

      bad = 0
      C_free   = iface_slip_conductance( A, df, dgap, lam, mu, 1.0e-13_dp, 1.0_dp )
      C_noslip = iface_slip_conductance( A, df, dgap, lam, mu, 1.0e13_dp,  1.0_dp )
      C_ref    = mu * A / ( df + dgap )
      if ( C_free > 1.0e-9_dp * C_ref ) bad = bad + 1            ! free slip
      if ( abs( C_noslip - C_ref ) > 1.0e-9_dp * C_ref ) bad = bad + 1

      ! eps weighting: the porous-side gap enters as eps*d_p (Betchen Eq.41)
      C_half = iface_slip_conductance( A, df, dgap, 0.0_dp, mu, 1.0_dp, 0.5_dp )
      C_ref2 = mu * A / ( df + 0.5_dp * dgap )
      if ( abs( C_half - C_ref2 ) > 1.0e-14_dp * C_ref2 ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0,a,3es12.4)') 'FAIL (2): BJ limits, bad=', bad, &
                                      C_free, C_noslip, C_ref
      else
         write(*,'(a)') 'PASS (2): BJ limits (free slip / stress continuity / eps)'
      end if
   end subroutine test_bj_limits

   !----------------------------------------------------------------------------
   ! (3) cases/beavers_joseph discretization: K = 1e-6 m^2 (lam = 1 mm), eps = 1,
   ! mu = 1.846e-5, dy = 0.25 mm -> the internal face sees the fluid cell at
   ! d_Pf = dy/2 = 0.125 mm.  The continuous limit of the discrete flux is the
   ! series-resistance interface of the README's two-layer reference:
   !   1/C = d_Pf/(mu*A) + lam/(mu*alpha*A)
   ! which must be exactly the resistance returned by the shared law with
   ! d_p = 0.  Both alpha = 1 and alpha = 2 (bj_a1 / bj_a2) are checked.
   !----------------------------------------------------------------------------
   subroutine test_bj_case_discretization( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp), parameter :: mu = 1.846e-5_dp, K = 1.0e-6_dp
      real(dp), parameter :: lam = 1.0e-3_dp, dPf = 0.125e-3_dp
      real(dp), parameter :: A = 1.0e-4_dp
      real(dp) :: R_law1, R_law2, R_ref1, R_ref2, C1, C2

      bad = 0
      C1 = iface_slip_conductance( A, dPf, 0.0_dp, lam, mu, 1.0_dp, 1.0_dp )
      C2 = iface_slip_conductance( A, dPf, 0.0_dp, lam, mu, 2.0_dp, 1.0_dp )
      R_law1 = 1.0_dp / C1
      R_law2 = 1.0_dp / C2
      ! README series form: 1/C = d_Pf/(mu*A) + lam/(mu*alpha*A)
      R_ref1 = dPf / ( mu * A ) + lam / ( mu * 1.0_dp * A )
      R_ref2 = dPf / ( mu * A ) + lam / ( mu * 2.0_dp * A )
      if ( abs( R_law1 - R_ref1 ) > 1.0e-14_dp * R_ref1 ) bad = bad + 1
      if ( abs( R_law2 - R_ref2 ) > 1.0e-14_dp * R_ref2 ) bad = bad + 1

      ! lambda = sqrt(K) must be what the caller passes
      if ( abs( lam - sqrt(K) ) > 1.0e-16_dp ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0,a,4es13.5)') 'FAIL (3): BJ case resistance, bad=', bad, &
                                      R_law1, R_ref1, R_law2, R_ref2
      else
         write(*,'(a,2es13.5,a,2es13.5)') 'PASS (3): BJ case series resistance al=1,2: ', &
               R_law1, R_law2, '  README: ', R_ref1, R_ref2
      end if
   end subroutine test_bj_case_discretization

   !----------------------------------------------------------------------------
   ! (4) `iface_slip_velocity` is the Dirichlet value that reproduces the
   ! stress-jump flux through a face whose flush conductance is D_por:
   !     D_por*(V - u_p) == C*(u_f - u_p)     (identity)
   !   plus the two limits (C = 0 -> V = u_p ; C = D_por -> V = u_f) and
   !   bracketing u_p <= V <= u_f for u_f > u_p.
   !----------------------------------------------------------------------------
   subroutine test_slip_velocity_identity( nfail )
      integer, intent(inout) :: nfail
      integer :: k, bad
      real(dp) :: C, Dpor, uf, up, V, lhs, rhs
      real(dp) :: Cs(4), Ds(4), ufs(4), ups(4)

      Cs  = (/ 1.0e-3_dp, 5.0e-4_dp, 9.0e-4_dp, 1.0e-5_dp /)
      Ds  = (/ 2.0e-3_dp, 5.0e-4_dp, 1.0e-3_dp, 1.0e-2_dp /)
      ufs = (/ 102.0_dp, -0.5_dp, 3.0_dp, 20.0_dp /)
      ups = (/ 0.30_dp,  0.10_dp, 1.5_dp, 1.0_dp  /)

      bad = 0
      do k = 1, 4
         C = Cs(k); Dpor = Ds(k); uf = ufs(k); up = ups(k)
         V   = iface_slip_velocity( C, Dpor, uf, up )
         lhs = Dpor * ( V - up )
         rhs = C    * ( uf - up )
         if ( abs( lhs - rhs ) > 1.0e-13_dp * abs(rhs) ) bad = bad + 1
         if ( V < min(uf,up) - 1.0e-12_dp .or. V > max(uf,up) + 1.0e-12_dp ) &
            bad = bad + 1
      end do
      ! limits
      if ( abs( iface_slip_velocity( 0.0_dp, 1.0e-3_dp, 10.0_dp, 2.0_dp ) &
                - 2.0_dp ) > 1.0e-14_dp ) bad = bad + 1
      if ( abs( iface_slip_velocity( 1.0e-3_dp, 1.0e-3_dp, 10.0_dp, 2.0_dp ) &
                - 10.0_dp ) > 1.0e-14_dp ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0)') 'FAIL (4): slip-velocity identity, bad=', bad
      else
         write(*,'(a)') 'PASS (4): slip-velocity identity + limits'
      end if
   end subroutine test_slip_velocity_identity

   !----------------------------------------------------------------------------
   ! (5) Betchen 2006 Eq.43/44 interface pressure:
   !     p_i = p_fluid - (mdot*(u.n)/A)*(1-eps)/eps
   ! With mdot = rho*A*u_n this is exactly p_fluid - rho*u_n^2*(1-eps)/eps.
   ! Numbers are the C_P_test coupling interface (one of its 60 faces,
   ! A = 0.6/60 = 0.01 m^2, rho = 1.177, eps = 0.5, u_n = 0.298 m/s as reported
   ! in VALIDATION_NOTES) plus a high-Re blowing check (u_n = 20 m/s).
   !----------------------------------------------------------------------------
   subroutine test_p_momentum( nfail )
      integer, intent(inout) :: nfail
      integer  :: bad
      real(dp), parameter :: rho = 1.177_dp, A = 1.0e-2_dp, eps = 0.5_dp
      real(dp), parameter :: pfl = 120.0_dp        ! gauge Pa
      real(dp) :: un, mdot, p, p_ref, p2, p3, p4

      bad = 0
      un   = 0.298_dp
      mdot = rho * A * un
      p    = iface_p_momentum( pfl, mdot, un, A, eps )
      p_ref = pfl - rho * un * un * ( 1.0_dp - eps ) / eps
      if ( abs( p - p_ref ) > 1.0e-12_dp ) bad = bad + 1
      ! magnitude is the dynamic head rho*u_n^2 = 0.1058 Pa for the transpiration
      ! velocity above, i.e. negligible against the ~100 Pa interface scale
      if ( abs( (p - pfl) + 0.104522_dp ) > 1.0e-5_dp ) bad = bad + 1

      ! high-Re blowing: rho*u_n^2*(1-eps)/eps = 1.177*400 = 470.8 Pa
      un   = 20.0_dp
      mdot = rho * A * un
      p2   = iface_p_momentum( pfl, mdot, un, A, eps )
      if ( abs( (p2 - pfl) + 470.8_dp ) > 1.0e-2_dp ) bad = bad + 1

      ! eps = 1 (no area change / pure fluid): pressure is continuous
      p3 = iface_p_momentum( pfl, mdot, un, A, 1.0_dp )
      if ( abs( p3 - pfl ) > 1.0e-14_dp ) bad = bad + 1
      ! degenerate inputs: zero area -> pass-through; eps -> 0 stays finite
      p4 = iface_p_momentum( pfl, mdot, un, 0.0_dp, eps )
      if ( abs( p4 - pfl ) > 1.0e-14_dp ) bad = bad + 1
      p4 = iface_p_momentum( pfl, mdot, un, A, 0.0_dp )
      if ( .not. ( abs(p4) < 1.0e30_dp ) ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0,a,3es13.5)') 'FAIL (5): Eq.44 momentum pressure, bad=', bad, &
              ' p/p_ref/p2=', p, p_ref, p2
      else
         write(*,'(a,2es13.5,a,es13.5,a)') 'PASS (5): Eq.44 dp = ', p-pfl, p2-pfl, &
              ' (rho*un^2*(1-eps)/eps =', &
              rho*20.0_dp**2 * (1.0_dp-eps)/eps, ')'
      end if
   end subroutine test_p_momentum

   !----------------------------------------------------------------------------
   ! (6) Two-sided blend (Betchen Sec.4.2): weights, limits and clamping
   !----------------------------------------------------------------------------
   subroutine test_p_blend( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp), parameter :: pf = 100.0_dp, pp = 104.0_dp
      real(dp) :: p

      bad = 0
      if ( abs( iface_p_blend( pf, pp, 0.0_dp ) - pf ) > 1.0e-14_dp ) bad = bad + 1
      if ( abs( iface_p_blend( pf, pp, 1.0_dp ) - pp ) > 1.0e-14_dp ) bad = bad + 1
      if ( abs( iface_p_blend( pf, pp, 0.5_dp ) - 0.5_dp*(pf+pp) ) > 1.0e-14_dp ) &
         bad = bad + 1
      ! out-of-range weights are clamped (never extrapolating into a new
      ! spurious pressure extremum)
      if ( abs( iface_p_blend( pf, pp, -3.0_dp ) - pf ) > 1.0e-14_dp ) bad = bad + 1
      if ( abs( iface_p_blend( pf, pp,  9.0_dp ) - pp ) > 1.0e-14_dp ) bad = bad + 1
      p = iface_p_blend( pf, pp, 0.25_dp )
      if ( abs( p - ( 0.75_dp*pf + 0.25_dp*pp ) ) > 1.0e-14_dp ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0)') 'FAIL (6): p blend, bad=', bad
      else
         write(*,'(a)') 'PASS (6): p blend weights/limits'
      end if
   end subroutine test_p_blend

   !----------------------------------------------------------------------------
   ! (7) p-mdot pair (Betchen Sec.4.2 "a small number of iterations is
   ! required"): local loop gain, stiffness of the undamped iteration and the
   ! quality of the deferred estimate.
   !
   ! Model of the pair, exactly as the coupling driver sees it:
   !   PPE   : mdot(p) = mdot0 - af*(p - p_prev)      (iface_p_flux_response)
   !   Eq.44 : p(mdot) = p_fl - (mdot*u_n/area)*(1-eps)/eps
   ! with the C_P_test interface discretization: rho = 1.177, eps = 0.5,
   ! A = 0.01 m^2 (one of 60 faces of a 0.6 m^2 patch), mu = 2.452e-3,
   ! K = 5e-7 (lam = 7.071e-4), slab cell V = 2.5e-5 m^3, d_p = 1.25e-3 m,
   ! a_P ~ 0.175 kg/s -> af ~ rho*(V/a_P)*A/d_p = 1.35e-3 kg/(s*Pa).
   !----------------------------------------------------------------------------
   subroutine test_p_mdot_pair( nfail )
      integer, intent(inout) :: nfail
      integer :: k, bad
      real(dp), parameter :: rho = 1.177_dp, A = 1.0e-2_dp, eps = 0.5_dp
      real(dp), parameter :: af = 1.35e-3_dp, pfl = 120.0_dp
      real(dp), parameter :: p_prev = 115.0_dp
      real(dp) :: un, mdot0, G, G_hi, p_def, p_num, p_k, om
      real(dp) :: d1, d3

      bad = 0

      ! ---- local loop gain, hand values -------------------------------------
      G = iface_p_gain( af, rho*A*0.298_dp, eps, rho, A )
      if ( abs( G - 0.08054_dp ) > 1.0e-4_dp ) bad = bad + 1
      ! high-Re blowing (u_n = 20 m/s) is stiff: G = 5.405 >> 1
      G_hi = iface_p_gain( af, rho*A*20.0_dp, eps, rho, A )
      if ( abs( G_hi - 5.405_dp ) > 1.0e-2_dp ) bad = bad + 1

      ! ---- C_P_test regime (G << 1): the undamped iteration CONTRACTS -------
      un    = 0.298_dp
      mdot0 = rho * A * un
      p_k   = p_prev
      do k = 1, 3
         p_num = iface_p_momentum( pfl, &
                    iface_p_flux_response( mdot0, af, p_k, p_prev ), &
                    iface_p_flux_response( mdot0, af, p_k, p_prev )/(rho*A), &
                    A, eps )
         if ( k == 1 ) d1 = abs( p_num - p_k )
         if ( k == 3 ) d3 = abs( p_num - p_k )
         p_k = p_num
      end do
      if ( .not. ( d3 < d1 ) ) bad = bad + 1

      ! ---- stiff regime (G >> 1): the pair has a REPELLING fixed point -----
      ! The composed map is monotone increasing with slope G > 1, so neither
      ! the raw nor the under-relaxed (omega in (0,1)) iteration can converge
      ! -- this is exactly why the driver only runs the sub-iteration when
      ! G < 1 and otherwise keeps the deferred (explicit) estimate, leaving the
      ! p-mdot loop to the outer coupling iteration.
      un    = 20.0_dp
      mdot0 = rho * A * un
      p_k   = p_prev
      do k = 1, 3
         p_num = iface_p_momentum( pfl, &
                    iface_p_flux_response( mdot0, af, p_k, p_prev ), &
                    iface_p_flux_response( mdot0, af, p_k, p_prev )/(rho*A), &
                    A, eps )
         if ( k == 1 ) d1 = abs( p_num - p_k )
         if ( k == 3 ) d3 = abs( p_num - p_k )
         p_k = p_num
      end do
      if ( .not. ( d3 > d1 ) ) bad = bad + 1             ! undamped grows
      om = 1.0_dp / ( 1.0_dp + G_hi )
      p_k = p_prev
      do k = 1, 3
         p_num = p_k + om * ( iface_p_momentum( pfl, &
                    iface_p_flux_response( mdot0, af, p_k, p_prev ), &
                    iface_p_flux_response( mdot0, af, p_k, p_prev )/(rho*A), &
                    A, eps ) - p_k )
         if ( k == 1 ) d1 = abs( p_num - p_k )
         if ( k == 3 ) d3 = abs( p_num - p_k )
         p_k = p_num
      end do
      if ( .not. ( d3 > d1 ) ) bad = bad + 1             ! damped still grows

      ! ---- deferred estimate vs the true fixed point when G << 1 ----------
      un    = 0.298_dp
      mdot0 = rho * A * un
      p_k   = p_prev
      do k = 1, 200
         p_k = iface_p_momentum( pfl, &
                  iface_p_flux_response( mdot0, af, p_k, p_prev ), &
                  iface_p_flux_response( mdot0, af, p_k, p_prev )/(rho*A), &
                  A, eps )
      end do
      p_def = iface_p_momentum( pfl, mdot0, un, A, eps )
      if ( abs( p_def - p_k ) > 0.02_dp * abs( p_def - p_prev ) ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0,a,4es12.4)') 'FAIL (7): p-mdot pair, bad=', bad, G, G_hi
      else
         write(*,'(a,es13.5,a,es13.5,a,es13.5)') 'PASS (7): p-mdot pair; G(C_P_test)=', G, &
              ' G(hi-Re)=', G_hi, '  deferred-fixedpoint dp=', p_def - p_k
      end if
   end subroutine test_p_mdot_pair

   !----------------------------------------------------------------------------
   ! (8) Zhang 2011 Eq.25 -- normal interface velocity is the two-sided
   ! conductance blend of the near-interface cell values:
   !     V_n = (G_p*<V>_n^p + G_f*V_n^fl)/(G_p + G_f),
   !     G_p = mu_e/(eps*d_p), G_f = mu/d_f
   ! Limits: d_f -> inf (no peer information) => V_n = <V>_n^p;
   !         d_p -> 0 (porous resistance vanishes) => V_n = <V>_n^p;
   !         d_f -> 0 (fluid resistance vanishes) => V_n = V_n^fl;
   ! plus the bracketing of the two inputs at finite gaps.
   !----------------------------------------------------------------------------
   subroutine test_zhang_vn( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp) :: mu, mu_e, eps, d_p, df, up, uf, Vn, Gp, Gf, exact

      mu   = 1.846e-5_dp
      eps  = 0.5_dp
      mu_e = mu / eps            ! code convention (see iface_slip_conductance)
      d_p   = 5.0e-4_dp
      df   = 2.0e-3_dp
      up   = 0.30_dp             ! porous-side normal velocity
      uf   = 0.05_dp             ! clear-fluid-side normal velocity

      bad = 0
      Gp    = mu_e / ( eps * d_p )
      Gf    = mu / df
      exact = ( Gp*up + Gf*uf ) / ( Gp + Gf )
      Vn    = iface_zhang_vn( up, uf, d_p, df, mu, mu_e, eps )
      if ( abs( Vn - exact ) > 1.0e-14_dp * abs(exact) ) bad = bad + 1
      if ( Vn < min(up,uf) - 1.0e-12_dp .or. Vn > max(up,uf) + 1.0e-12_dp ) &
         bad = bad + 1

      ! limits
      if ( abs( iface_zhang_vn( up, uf, d_p, 1.0e30_dp, mu, mu_e, eps ) &
                - up ) > 1.0e-12_dp ) bad = bad + 1
      if ( abs( iface_zhang_vn( up, uf, 1.0e-30_dp, df, mu, mu_e, eps ) &
                - up ) > 1.0e-12_dp ) bad = bad + 1
      if ( abs( iface_zhang_vn( up, uf, d_p, 1.0e-30_dp, mu, mu_e, eps ) &
                - uf ) > 1.0e-12_dp ) bad = bad + 1
      ! equal gaps: the porous side is stiffer by 1/eps^2 with mu_e = mu/eps
      exact = ( up/eps**2 + uf ) / ( 1.0_dp/eps**2 + 1.0_dp )
      if ( abs( iface_zhang_vn( up, uf, d_p, d_p, mu, mu_e, eps ) - exact ) &
           > 1.0e-13_dp ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0)') 'FAIL (8): Zhang Eq.25 normal blend, bad=', bad
      else
         write(*,'(a)') 'PASS (8): Zhang Eq.25 normal velocity blend + limits'
      end if
   end subroutine test_zhang_vn

   !----------------------------------------------------------------------------
   ! (9) Zhang 2011 Eq.26 -- tangential stress-jump balance.
   !   beta1 = 0             : V_t = Q/(G_p+G_f+beta*mu_e/lam)
   !   beta -> inf (alpha->0): V_t -> 0                     (no slip)
   !   beta  = 0             : the same harmonic blend as Eq.25 (stress
   !                           continuity)
   !   beta1 > 0             : the excess inertia always REDUCES |V_t|, and the
   !                           returned value solves the quadratic exactly
   !   no real root          : falls back to the linear root (no NaN)
   !----------------------------------------------------------------------------
   subroutine test_zhang_vt( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp) :: mu, mu_e, eps, d_p, df, lam, rho, up, uf, Gp, Gf
      real(dp) :: V0, V1, Vb, Q, B, disc, root, big

      mu   = 1.846e-5_dp
      eps  = 0.5_dp
      mu_e = mu / eps
      d_p   = 5.0e-4_dp
      df   = 2.0e-3_dp
      lam  = 1.0e-3_dp
      rho  = 1.177_dp
      up   = 0.10_dp
      uf   = 12.0_dp

      Gp = mu_e / ( eps * d_p )
      Gf = mu / df
      Q  = Gp*up + Gf*uf

      bad = 0
      ! beta = 0, beta1 = 0 : pure harmonic blend (identical to Eq.25)
      V0 = iface_zhang_vt( up, uf, d_p, df, mu, mu_e, eps, lam, &
                           0.0_dp, 0.0_dp, rho )
      if ( abs( V0 - Q/(Gp+Gf) ) > 1.0e-14_dp*abs(V0) .or. &
           abs( V0 - iface_zhang_vn( up, uf, d_p, df, mu, mu_e, eps ) ) &
             > 1.0e-14_dp*abs(V0) ) bad = bad + 1
      if ( V0 < min(up,uf) - 1.0e-12_dp .or. V0 > max(up,uf) + 1.0e-12_dp ) &
         bad = bad + 1

      ! finite beta: Robin (excess slip resistance) -> strictly below the blend
      V1 = iface_zhang_vt( up, uf, d_p, df, mu, mu_e, eps, lam, &
                           1.0_dp, 0.0_dp, rho )
      B  = Gp + Gf + mu_e/lam
      if ( abs( V1 - Q/B ) > 1.0e-14_dp*abs(V1) ) bad = bad + 1
      if ( V1 >= V0 ) bad = bad + 1

      ! no slip: beta -> inf drives V_t -> 0
      big = iface_zhang_vt( up, uf, d_p, df, mu, mu_e, eps, lam, &
                            1.0e12_dp, 0.0_dp, rho )
      if ( abs( big ) > 1.0e-6_dp ) bad = bad + 1

      ! inertial term: exact stable quadratic root, and a smaller slip
      Vb   = iface_zhang_vt( up, uf, d_p, df, mu, mu_e, eps, lam, &
                             0.0_dp, 2.0_dp, rho )
      B    = Gp + Gf
      disc = B*B + 4.0_dp*2.0_dp*rho*Q
      root = 2.0_dp*Q / ( B + sqrt(disc) )
      if ( abs( Vb - root ) > 1.0e-14_dp*abs(root) ) bad = bad + 1
      if ( abs( Vb ) >= abs( V0 ) ) bad = bad + 1

      ! unreachable inertial term (disc < 0): linear fallback, finite
      B   = Gp + Gf
      V1  = iface_zhang_vt( -up, -uf, d_p, df, mu, mu_e, eps, lam, &
                            0.0_dp, 1.0e12_dp, rho )
      if ( V1 /= V1 ) bad = bad + 1                         ! NaN guard
      if ( abs( V1 - ( -Q/B ) ) > 1.0e-12_dp*abs(Q/B) ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0)') 'FAIL (9): Zhang Eq.26 tangential jump, bad=', bad
      else
         write(*,'(a)') 'PASS (9): Zhang Eq.26 tangential stress jump + limits'
      end if
   end subroutine test_zhang_vt

   !----------------------------------------------------------------------------
   ! (10) Zhang 2011 Eqs.22+25+26 -- the vector interface velocity.
   !    - both cell velocities normal to the interface => V purely normal and
   !      equal to the Eq.25 blend;
   !    - pure tangential slip with beta = 0 => the direction of Q_t is kept
   !      and the interface normal component carries the Eq.25 blend;
   !    - beta1 = 0 reconstitutes the vector balance
   !      (G_p+G_f+beta*mu_e/lam) V_t = Q_t on the tangential projection.
   !----------------------------------------------------------------------------
   subroutine test_zhang_velocity_vector( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp) :: mu, mu_e, eps, d_p, df, lam, rho, beta
      real(dp) :: nv(3), up(3), uf(3), V(3), Vt(3), Qt(3), Gp, Gf, B
      real(dp) :: nrm

      mu   = 1.846e-5_dp
      eps  = 0.5_dp
      mu_e = mu / eps
      d_p   = 5.0e-4_dp
      df   = 2.0e-3_dp
      lam  = 1.0e-3_dp
      rho  = 1.177_dp
      beta = 1.0_dp

      Gp = mu_e / ( eps * d_p )
      Gf = mu / df
      B  = Gp + Gf + beta*mu_e/lam

      bad = 0
      ! (a) pure normal
      nv  = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
      up  = 0.30_dp * nv
      uf  = 0.05_dp * nv
      V   = iface_zhang_velocity( up, uf, nv, d_p, df, mu, mu_e, eps, &
                                  lam, beta, 0.0_dp, rho )
      if ( abs( V(1) ) > 1.0e-14_dp .or. abs( V(3) ) > 1.0e-14_dp ) bad = bad + 1
      if ( abs( V(2) - iface_zhang_vn( 0.30_dp, 0.05_dp, d_p, df, mu, mu_e, &
                                       eps ) ) > 1.0e-14_dp ) bad = bad + 1

      ! (b) pure tangential, beta = 0: direction of Q_t preserved, n-component
      !     still the Eq.25 blend (here zero)
      nv = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
      up = (/ 2.0_dp, -3.0_dp, 0.0_dp /)
      uf = (/ 40.0_dp, -1.0_dp, 0.0_dp /)
      V  = iface_zhang_velocity( up, uf, nv, d_p, df, mu, mu_e, eps, &
                                 lam, 0.0_dp, 0.0_dp, rho )
      if ( abs( V(3) ) > 1.0e-14_dp ) bad = bad + 1
      Qt = Gp*up + Gf*uf
      if ( norm2( V ) <= 0.0_dp .or. &
           abs( V(1)*Qt(2) - V(2)*Qt(1) ) > 1.0e-12_dp * norm2(V)*norm2(Qt) ) &
         bad = bad + 1
      if ( abs( V(1) - Qt(1)/(Gp+Gf) ) > 1.0e-14_dp*abs(V(1)) .or. &
           abs( V(2) - Qt(2)/(Gp+Gf) ) > 1.0e-14_dp*abs(V(2)) ) bad = bad + 1

      ! (c) reconstitute the vector balance with the slip resistance
      V  = iface_zhang_velocity( up, uf, nv, d_p, df, mu, mu_e, eps, &
                                 lam, beta, 0.0_dp, rho )
      Vt = V - dot_product( V, nv )*nv              ! tangential part (n = z)
      if ( norm2( B*Vt - Qt ) > 1.0e-14_dp * norm2(Qt) ) bad = bad + 1
      nrm = iface_zhang_vn( dot_product(up,nv), dot_product(uf,nv), d_p, df, &
                            mu, mu_e, eps )
      if ( abs( dot_product( V, nv ) - nrm ) > 1.0e-14_dp ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0)') 'FAIL (10): Zhang vector interface velocity, bad=', bad
      else
         write(*,'(a)') 'PASS (10): Zhang Eq.22/25/26 vector interface velocity'
      end if
   end subroutine test_zhang_velocity_vector

   !----------------------------------------------------------------------------
   ! (11) Zhang 2011 Eqs.27-29 -- interface temperature / flux split.
   !   (a) Eq.27 : eps*T_fi + (1-eps)*T_si = T_fl = <T>^p
   !   (b) Eqs.28/29 : the phase fluxes are exactly eps*F and (1-eps)*F, and
   !       their sum equals the clear-fluid flux k_f(T_fl-T_flP)/d_f
   !       (energy conservation across the interface);
   !   (c) LTE reduction (T_f = T_s, k_fe = eps*k_e, k_se = (1-eps)*k_e):
   !       all three interface values collapse onto the series-resistance
   !       average (k_f/d_f*T_flP + k_e/d_p*T_p)/(k_f/d_f + k_e/d_p);
   !   (d) the volume-average helper returns Eq.27's value;
   !   (e) a degenerate (non-positive) denominator falls back to the one-sided
   !       Dirichlet with zero flux (no NaN).
   !----------------------------------------------------------------------------
   subroutine test_zhang_T( nfail )
      integer, intent(inout) :: nfail
      integer :: bad
      real(dp) :: eps, kf, kfe, kse, d_p, df, ke
      real(dp) :: TflP, TfP, TsP, Tfl, Tfi, Tsi, F, Ff, Fs, Fl, exact

      eps  = 0.5_dp
      kf   = 0.026_dp                 ! air
      kfe  = eps * kf
      kse  = ( 1.0_dp - eps ) * 0.5_dp
      d_p   = 5.0e-4_dp
      df   = 1.0e-3_dp
      TflP = 900.0_dp
      TfP  = 300.0_dp
      TsP  = 380.0_dp

      bad = 0
      call iface_zhang_T( TflP, TfP, TsP, eps, kf, kfe, kse, df, d_p, &
                          Tfl, Tfi, Tsi, F )
      ! (a) Eq.27
      if ( abs( eps*Tfi + (1.0_dp-eps)*Tsi - Tfl ) > 1.0e-10_dp ) bad = bad + 1
      ! (b) flux split + energy conservation
      Ff = kfe * ( Tfi - TfP ) / d_p
      Fs = kse * ( Tsi - TsP ) / d_p
      Fl = kf  * ( TflP - Tfl ) / df
      if ( abs( Ff - eps*F ) > 1.0e-10_dp*max(abs(Ff),1.0_dp) ) bad = bad + 1
      if ( abs( Fs - (1.0_dp-eps)*F ) > 1.0e-10_dp*max(abs(Fs),1.0_dp) ) &
         bad = bad + 1
      if ( abs( Ff + Fs - Fl ) > 1.0e-9_dp*max(abs(Fl),1.0_dp) ) bad = bad + 1
      ! heat flows from the hot clear fluid into the colder porous region
      if ( F <= 0.0_dp ) bad = bad + 1
      ! and the interface values are bracketed by the two sides
      if ( Tfl > TflP .or. Tfl < TfP .or. Tfi < TfP .or. Tsi < TsP ) &
         bad = bad + 1

      ! (c) LTE reduction: all three values = the series-resistance average
      ke = 0.8_dp
      call iface_zhang_T( 700.0_dp, 400.0_dp, 400.0_dp, eps, kf, &
                          eps*ke, (1.0_dp-eps)*ke, df, d_p, &
                          Tfl, Tfi, Tsi, F )
      exact = ( kf/df*700.0_dp + ke/d_p*400.0_dp ) / ( kf/df + ke/d_p )
      if ( abs( Tfl - exact ) > 1.0e-9_dp .or. &
           abs( Tfi - exact ) > 1.0e-9_dp .or. &
           abs( Tsi - exact ) > 1.0e-9_dp ) bad = bad + 1

      ! (d) volume-average helper
      if ( abs( iface_volavg_T( 400.0_dp, 200.0_dp, 0.5_dp ) - 300.0_dp ) &
           > 1.0e-14_dp ) bad = bad + 1
      if ( abs( iface_volavg_T( 400.0_dp, 400.0_dp, 0.3_dp ) - 400.0_dp ) &
           > 1.0e-14_dp ) bad = bad + 1
      call iface_zhang_T( TflP, TfP, TsP, eps, kf, kfe, kse, df, d_p, &
                          Tfl, Tfi, Tsi, F )
      if ( abs( Tfl - iface_volavg_T( Tfi, Tsi, eps ) ) > 1.0e-12_dp ) &
         bad = bad + 1

      ! (e) a vanishing fluid-side conductance (k_f -> 0 is represented by the
      !     degenerate guard) leaves the one-sided Dirichlet with zero flux,
      !     finite and NaN-free; a *large* porous-side resistance (huge d_p)
      !     still gives the physical, bracketed solution with a tiny flux.
      call iface_zhang_T( 900.0_dp, 300.0_dp, 380.0_dp, eps, 0.0_dp, kfe, kse, &
                          df, d_p, Tfl, Tfi, Tsi, F )
      if ( Tfl /= Tfl ) bad = bad + 1                      ! NaN guard
      if ( abs( Tfl - 900.0_dp ) > 1.0e-12_dp ) bad = bad + 1
      if ( abs( Tfi - 300.0_dp ) > 1.0e-12_dp ) bad = bad + 1
      if ( abs( Tsi - 380.0_dp ) > 1.0e-12_dp ) bad = bad + 1
      if ( abs( F ) > 1.0e-30_dp ) bad = bad + 1

      call iface_zhang_T( 900.0_dp, 300.0_dp, 380.0_dp, eps, kf, kfe, kse, &
                          df, 1.0e6_dp, Tfl, Tfi, Tsi, F )
      if ( Tfl /= Tfl .or. F /= F ) bad = bad + 1          ! NaN guard
      if ( F > 2.0_dp ) bad = bad + 1                      ! nearly insulating
      if ( Tfl > 900.0_dp .or. Tfl < 340.0_dp ) bad = bad + 1

      if ( bad > 0 ) then
         nfail = nfail + 1
         write(*,'(a,i0)') 'FAIL (11): Zhang Eq.27-29 interface T, bad=', bad
      else
         write(*,'(a)') 'PASS (11): Zhang Eq.27-29 interface temperature + split'
      end if
   end subroutine test_zhang_T

end program test_iface_law



