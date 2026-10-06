!===============================================================================
! mod_forces.f90 -- Pressure/viscous force integration on wall zones
!
! Provides:
!   compute_forces       -- integrate pressure (and viscous if mu>0) force on
!                           each wall zone; report total Fx/Fy/Fz and, when a
!                           pressure-far-field BC is present, the drag and lift
!                           coefficients CD/CL using q_inf = 0.5*rho*|U_inf|^2
!                           and A_ref = (x_max-x_min)*(z_max-z_min) (frontal
!                           projected area of the body).
!   write_cp_distribution -- write "theta(deg)  Cp  Cp_exact  x  y  z" for every
!                           wall face, with the analytical potential-flow Cp
!                           Cp = 1 - 4 sin^2(theta) (cylinder, no circulation).
!
! Force sign convention:
!   g%sf(:,i) is the face area vector pointing from the owner cell c0 outward
!   (toward the wall, i.e. into the body). Pressure p in the fluid pushes the
!   wall in the +Sf direction, so the pressure force on the body is
!       F = sum_i (p_i - p_inf) * Sf_i
!   and the drag (along +x, the free-stream direction) is Fx = sum (p-p_inf)*Sf_x.
!   For potential flow around a cylinder, Fx = 0 (d'Alembert paradox).
!
! All routines are rank-0 only (operate on the global mesh m_g and global
! fld_g after gathering); they make no MPI calls.
!===============================================================================
module mod_uns_forces
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_geometry
   use mod_uns_control
   use mod_uns_bc
   use mod_uns_fields
   implicit none
   private
   public :: compute_forces, write_cp_distribution, compute_separation

contains

   !----------------------------------------------------------------------------
   ! Integrate pressure (and optional viscous) force on each wall zone.
   ! Returns total Fx/Fy/Fz, plus CD/CL when a far-field BC exists.
   !----------------------------------------------------------------------------
   subroutine compute_forces( m, g, ctrl, bcs, fld, &
                              Fx, Fy, Fz, CD, CL, u_inf, A_ref, p_inf )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(ctrl_t),   intent(in)  :: ctrl
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld
      real(dp),       intent(out) :: Fx, Fy, Fz
      real(dp),       intent(out) :: CD, CL
      real(dp),       intent(out) :: u_inf, A_ref, p_inf

      integer  :: i, gi, c0
      integer  :: nwall
      real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
      real(dp) :: D_char, L_char, q_dyn, p_gauge
      real(dp) :: nx, ny, nz, tx, ty, tz, dn, ut, tau
      real(dp) :: Fx_p, Fy_p, Fz_p, Fx_v, Fy_v, Fz_v
      real(dp), allocatable :: fx_z(:), fy_z(:), fz_z(:)
      real(dp), allocatable :: vfx_z(:), vfy_z(:), vfz_z(:)

      Fx = 0.0_dp; Fy = 0.0_dp; Fz = 0.0_dp
      CD = 0.0_dp; CL = 0.0_dp
      u_inf = 0.0_dp; A_ref = 0.0_dp; p_inf = 0.0_dp
      Fx_p = 0.0_dp; Fy_p = 0.0_dp; Fz_p = 0.0_dp
      Fx_v = 0.0_dp; Fy_v = 0.0_dp; Fz_v = 0.0_dp

      ! ---- free-stream reference from far-field BC, if any ----------------
      do gi = 1, bcs%nb
         if ( bcs%gb(gi)%btype == BC_FARFIELD ) then
            u_inf = norm2( bcs%gb(gi)%uvel )
            p_inf = bcs%gb(gi)%pval
            exit
         end if
      end do

      ! ---- count wall zones; bail out if none ------------------------------
      nwall = 0
      do gi = 1, bcs%nb
         if ( bcs%gb(gi)%btype == BC_WALL ) nwall = nwall + 1
      end do
      if ( nwall == 0 ) return

      ! ---- bounding box of wall-face centroids (for D) and mesh nodes (L) ---
      ! For a 2D-extruded single-layer mesh the wall face centroids all share
      ! the same z (mid-layer), so zmax-zmin=0. Use the mesh NODE z-range for
      ! the span L instead; this is the physical extrusion length.
      xmin = huge(1.0_dp); xmax = -huge(1.0_dp)
      ymin = huge(1.0_dp); ymax = -huge(1.0_dp)
      zmin = huge(1.0_dp); zmax = -huge(1.0_dp)
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         if ( bcs%gb(gi)%btype /= BC_WALL ) cycle
         xmin = min( xmin, g%xf(1,i) ); xmax = max( xmax, g%xf(1,i) )
         ymin = min( ymin, g%xf(2,i) ); ymax = max( ymax, g%xf(2,i) )
         zmin = min( zmin, g%xf(3,i) ); zmax = max( zmax, g%xf(3,i) )
      end do
      ! Frontal projected area for a cylinder in x-flow with z-axis along span:
      !   A_ref = D * L = (x_max - x_min) * (z_node_max - z_node_min)
      D_char  = xmax - xmin
      ! Use mesh node z-range for L (handles single-layer extrusion correctly)
      L_char  = maxval( m%x(3,:) ) - minval( m%x(3,:) )
      A_ref   = D_char * L_char

      ! ---- pressure + viscous force on each wall zone -----------------------
      ! Pressure force: F_p = sum (p - p_inf) * Sf
      ! Viscous force (no-slip wall): tau_w = mu * u_cell_t / dn,
      !   F_v = sum tau_w * t_hat * area, where t_hat is the wall tangent.
      allocate( fx_z(bcs%nb), fy_z(bcs%nb), fz_z(bcs%nb) )
      allocate( vfx_z(bcs%nb), vfy_z(bcs%nb), vfz_z(bcs%nb) )
      fx_z = 0.0_dp; fy_z = 0.0_dp; fz_z = 0.0_dp
      vfx_z = 0.0_dp; vfy_z = 0.0_dp; vfz_z = 0.0_dp

      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         if ( bcs%gb(gi)%btype /= BC_WALL ) cycle
         c0 = m%f(i)%c0
         ! pressure force
         p_gauge = fld%p(c0) - p_inf
         fx_z(gi) = fx_z(gi) + p_gauge * g%sf(1,i)
         fy_z(gi) = fy_z(gi) + p_gauge * g%sf(2,i)
         fz_z(gi) = fz_z(gi) + p_gauge * g%sf(3,i)
         ! viscous force (only if mu > 0)
         if ( ctrl%mu > 0.0_dp ) then
            nx = g%sf(1,i) / norm2( g%sf(:,i) )
            ny = g%sf(2,i) / norm2( g%sf(:,i) )
            nz = g%sf(3,i) / norm2( g%sf(:,i) )
            ! tangent direction = n x z_hat (dominant in-plane tangent)
            tx = ny; ty = -nx; tz = 0.0_dp
            dn = norm2( g%xf(:,i) - g%xc(:,c0) )
            ut = fld%u(1,c0)*tx + fld%u(2,c0)*ty + fld%u(3,c0)*tz
            tau = ctrl%mu * ut / dn
            ! force on the wall: fluid drags wall in +t direction
            vfx_z(gi) = vfx_z(gi) + tau * tx * g%area(i)
            vfy_z(gi) = vfy_z(gi) + tau * ty * g%area(i)
            vfz_z(gi) = vfz_z(gi) + tau * tz * g%area(i)
         end if
      end do

      Fx_p = sum( fx_z ); Fy_p = sum( fy_z ); Fz_p = sum( fz_z )
      Fx_v = sum( vfx_z ); Fy_v = sum( vfy_z ); Fz_v = sum( vfz_z )
      Fx = Fx_p + Fx_v
      Fy = Fy_p + Fy_v
      Fz = Fz_p + Fz_v

      ! ---- coefficients -----------------------------------------------------
      if ( u_inf > 0.0_dp .and. A_ref > 0.0_dp ) then
         q_dyn = 0.5_dp * ctrl%rho * u_inf**2
         CD = Fx / ( q_dyn * A_ref )
         CL = Fy / ( q_dyn * A_ref )
      end if

      ! ---- report -----------------------------------------------------------
      write(*,'(a)') ''
      write(*,'(a)') '--- Force summary (pressure + viscous on wall zones) ---'
      write(*,'(a,es14.6)') '  Free-stream U_inf  : ', u_inf
      write(*,'(a,es14.6)') '  Free-stream p_inf  : ', p_inf
      write(*,'(a,es14.6)') '  Reference area     : ', A_ref
      write(*,'(a,es14.6)') '  Pressure  Fx (drag): ', Fx_p
      write(*,'(a,es14.6)') '  Viscous   Fx (drag): ', Fx_v
      write(*,'(a,es14.6)') '  Total     Fx (drag): ', Fx
      write(*,'(a,es14.6)') '  Total     Fy (lift): ', Fy
      write(*,'(a,es14.6)') '  Total     Fz       : ', Fz
      if ( u_inf > 0.0_dp .and. A_ref > 0.0_dp ) then
         write(*,'(a,es14.6)') '  CD_p (pressure)    : ', Fx_p / ( q_dyn * A_ref )
         write(*,'(a,es14.6)') '  CD_v (viscous)     : ', Fx_v / ( q_dyn * A_ref )
         write(*,'(a,es14.6)') '  CD  = total Fx/(q*A): ', CD
         write(*,'(a,es14.6)') '  CL  = total Fy/(q*A): ', CL
         write(*,'(a,es14.6)') '  CD lit. (Re=40)    : ', 1.5_dp
         write(*,'(a,es14.6)') '  |CD - 1.5|         : ', abs( CD - 1.5_dp )
      end if

      deallocate( fx_z, fy_z, fz_z, vfx_z, vfy_z, vfz_z )

   end subroutine compute_forces

   !----------------------------------------------------------------------------
   ! Write Cp distribution on wall faces for comparison with the analytical
   ! potential-flow solution Cp(theta) = 1 - 4 sin^2(theta) (cylinder, no
   ! circulation). theta is the polar angle measured from +x axis (CCW), so the
   ! front stagnation point at (-R,0) corresponds to theta=pi, where
   ! sin^2(theta)=0 and Cp=1; the top at (0,+R) corresponds to theta=pi/2,
   ! sin^2(theta)=1, Cp=-3 -- which matches the textbook formula 1-4sin^2(theta)
   ! measured from the front stagnation point.
   !----------------------------------------------------------------------------
   subroutine write_cp_distribution( m, g, ctrl, bcs, fld, filename )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(ctrl_t),   intent(in)  :: ctrl
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld
      character(len=*), intent(in) :: filename

      integer  :: u, ios, i, gi, c0, n
      real(dp) :: u_inf, p_inf, q_dyn, Cp, Cp_exact, theta
      real(dp) :: xcyl, ycyl, x, y, z

      ! ---- free-stream reference from far-field BC --------------------------
      u_inf = 0.0_dp; p_inf = 0.0_dp
      do gi = 1, bcs%nb
         if ( bcs%gb(gi)%btype == BC_FARFIELD ) then
            u_inf = norm2( bcs%gb(gi)%uvel )
            p_inf = bcs%gb(gi)%pval
            exit
         end if
      end do
      if ( u_inf <= 0.0_dp ) then
         write(*,'(a)') 'WARNING: no far-field BC; skip Cp distribution write'
         return
      end if
      q_dyn = 0.5_dp * ctrl%rho * u_inf**2

      ! ---- cylinder center = centroid of wall face centroids in x-y --------
      xcyl = 0.0_dp; ycyl = 0.0_dp; n = 0
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         if ( bcs%gb(gi)%btype /= BC_WALL ) cycle
         xcyl = xcyl + g%xf(1,i)
         ycyl = ycyl + g%xf(2,i)
         n    = n + 1
      end do
      if ( n == 0 ) return
      xcyl = xcyl / real( n, dp )
      ycyl = ycyl / real( n, dp )

      ! ---- write Cp table --------------------------------------------------
      open( newunit = u, file = trim(filename), status = 'replace', &
            action = 'write', iostat = ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'WARNING: cannot open Cp output file: ' // trim(filename)
         return
      end if

      write(u,'(a)')     '# theta(deg)   Cp           Cp_exact     x            y            z'
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         if ( bcs%gb(gi)%btype /= BC_WALL ) cycle
         c0 = m%f(i)%c0
         x = g%xf(1,i); y = g%xf(2,i); z = g%xf(3,i)
         theta    = atan2( y - ycyl, x - xcyl )
         Cp       = ( fld%p(c0) - p_inf ) / q_dyn
         Cp_exact = 1.0_dp - 4.0_dp * sin( theta )**2
         write(u,'(f10.3,2(1x,es12.4),3(1x,es11.3))') &
            theta * 180.0_dp / pi, Cp, Cp_exact, x, y, z
      end do
      close( u )

      write(*,'(a,a)')     '  Cp distribution written to : ', trim(filename)
      write(*,'(a,i0)')    '  wall faces                : ', n
      write(*,'(a)')       '  Cp_exact = 1 - 4 sin^2(theta)  (potential flow)'

   end subroutine write_cp_distribution

   !----------------------------------------------------------------------------
   ! Compute separation angle and wake length for a cylinder in cross-flow.
   !
   ! Separation angle: where the wall shear stress changes sign. For a no-slip
   ! wall tau_w ~ mu * u_cell_tangent / delta_n, so we track the tangential
   ! velocity of the owner cell along the wall and find its zero crossings.
   ! The tangential direction is t = n x z_hat = (n_y, -n_x, 0), where n is
   ! the outward face normal; positive t points counter-clockwise (CCW) when
   ! viewed from +z.
   !
   ! theta is measured CCW from +x axis:
   !   front stagnation  (-R, 0) -> theta = 180 deg
   !   top               (0, +R) -> theta =  90 deg
   !   rear stagnation   (+R, 0) -> theta =   0 deg
   ! The "separation angle from front stagnation" = 180 - theta_sep (deg).
   !
   ! Wake length: maximum x-coordinate of cells with u_x < 0 (reverse flow)
   ! behind the cylinder, minus the cylinder center x.
   !----------------------------------------------------------------------------
   subroutine compute_separation( m, g, ctrl, bcs, fld )
      type(mesh_t),   intent(in)  :: m
      type(geom_t),   intent(in)  :: g
      type(ctrl_t),   intent(in)  :: ctrl
      type(bc_t),     intent(in)  :: bcs
      type(fields_t), intent(in)  :: fld

      integer  :: i, gi, c0, n, k, j
      real(dp) :: xcyl, ycyl, x, y, nx, ny, tx, ty, ut
      real(dp) :: theta_sep_top, theta_sep_bot, x_wake_end
      real(dp), parameter :: PI = 3.141592653589793238462643_dp
      logical  :: found_top, found_bot
      ! per-wall-face arrays
      integer,  allocatable :: idx(:)
      real(dp), allocatable :: th(:), ut_arr(:)
      real(dp) :: ut_prev, th_prev, tmp
      integer  :: itmp

      theta_sep_top = 0.0_dp; theta_sep_bot = 0.0_dp
      found_top = .false.; found_bot = .false.

      ! ---- cylinder center from wall face centroids -------------------------
      xcyl = 0.0_dp; ycyl = 0.0_dp; n = 0
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         if ( bcs%gb(gi)%btype /= BC_WALL ) cycle
         xcyl = xcyl + g%xf(1,i)
         ycyl = ycyl + g%xf(2,i)
         n    = n + 1
      end do
      if ( n == 0 ) return
      xcyl = xcyl / real( n, dp )
      ycyl = ycyl / real( n, dp )

      ! ---- collect wall face angle and tangential velocity ------------------
      allocate( idx(n), th(n), ut_arr(n) )
      j = 0
      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         if ( bcs%gb(gi)%btype /= BC_WALL ) cycle
         c0 = m%f(i)%c0
         x = g%xf(1,i); y = g%xf(2,i)
         nx = g%sf(1,i) / norm2( g%sf(:,i) )
         ny = g%sf(2,i) / norm2( g%sf(:,i) )
         tx = ny; ty = -nx                 ! tangent = n x z_hat (CCW positive)
         ut = fld%u(1,c0)*tx + fld%u(2,c0)*ty
         j = j + 1
         idx(j) = j
         th(j)  = atan2( y - ycyl, x - xcyl )   ! -pi .. pi
         ut_arr(j) = ut
      end do

      ! ---- sort by theta (insertion sort, n is small ~80) -------------------
      do i = 2, n
         itmp = idx(i); tmp = th(i); ut = ut_arr(i)
         j = i - 1
         do while ( j >= 1 .and. th(j) > tmp )
            idx(j+1) = idx(j); th(j+1) = th(j); ut_arr(j+1) = ut_arr(j)
            j = j - 1
         end do
         idx(j+1) = itmp; th(j+1) = tmp; ut_arr(j+1) = ut
      end do

      ! ---- find sign change of ut along the sorted perimeter ----------------
      ! Sorted theta goes -pi -> 0 -> +pi (bottom-front -> rear -> top-front).
      ! Tangent t = n x z_hat with n pointing toward the body.
      !   bottom (theta<0): attached flow (+x along bottom) gives ut>0;
      !                     separated wake gives ut<0. Look for + -> -.
      !   top    (theta>0): sorted from rear (theta~0, wake ut>0) to front
      !                     (theta~pi, attached ut<0). Look for + -> -.
      ! Both sides look for a positive-to-negative crossing. Faces with
      ! |ut| < 1e-4 (near stagnation/separation) are skipped entirely and do
      ! NOT update ut_prev, so the crossing is detected between the last
      ! robust positive and the first robust negative.
      ut_prev = 0.0_dp; th_prev = 0.0_dp
      do i = 1, n
         if ( abs( ut_arr(i) ) < 1.0e-4_dp ) cycle   ! skip noisy faces
         if ( ut_arr(i) < 0.0_dp .and. ut_prev > 0.0_dp ) then
            if ( th(i) < 0.0_dp .and. .not. found_bot ) then
               theta_sep_bot = th_prev + (0.0_dp - ut_prev) &
                  * ( th(i) - th_prev ) / ( ut_arr(i) - ut_prev )
               found_bot = .true.
            else if ( th(i) > 0.0_dp .and. .not. found_top ) then
               theta_sep_top = th_prev + (0.0_dp - ut_prev) &
                  * ( th(i) - th_prev ) / ( ut_arr(i) - ut_prev )
               found_top = .true.
            end if
         end if
         ut_prev = ut_arr(i); th_prev = th(i)
      end do

      ! ---- wake length: max x of reverse-flow cells behind cylinder --------
      x_wake_end = xcyl
      do k = 1, m%ncells
         if ( fld%u(1,k) < 0.0_dp .and. g%xc(1,k) > xcyl ) then
            x_wake_end = max( x_wake_end, g%xc(1,k) )
         end if
      end do

      ! ---- report -----------------------------------------------------------
      write(*,'(a)') ''
      write(*,'(a)') '--- Separation and wake (Re = rho*U*D/mu) ---'
      write(*,'(a,es12.4)') '  Re (based on D)        : ', &
         ctrl%rho * 1.0_dp * 1.0_dp / max( ctrl%mu, tiny(1.0_dp) )
      if ( found_top ) then
         write(*,'(a,f8.2,a)') '  Top separation theta   : ', &
            theta_sep_top * 180.0_dp / PI, '  deg (from +x, CCW)'
         write(*,'(a,f8.2,a)') '  Top sep. from front stg: ', &
            180.0_dp - theta_sep_top * 180.0_dp / PI, &
            '  deg (lit. ~125-130 at Re=40)'
      else
         write(*,'(a)') '  Top separation         : not detected'
      end if
      if ( found_bot ) then
         write(*,'(a,f8.2,a)') '  Bot separation theta   : ', &
            theta_sep_bot * 180.0_dp / PI, '  deg (from +x, CCW)'
         write(*,'(a,f8.2,a)') '  Bot sep. from front stg: ', &
            180.0_dp - abs( theta_sep_bot ) * 180.0_dp / PI, &
            '  deg (lit. ~125-130 at Re=40)'
      else
         write(*,'(a)') '  Bot separation         : not detected'
      end if
      write(*,'(a,f8.3,a)') '  Wake length (x_wake - xcyl)/D : ', &
         ( x_wake_end - xcyl ) / 1.0_dp, '  (lit. ~1.5-2.0 at Re=40)'

      deallocate( idx, th, ut_arr )

   end subroutine compute_separation

end module mod_uns_forces
