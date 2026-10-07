!===============================================================================
! mod_interface_match.f90 -- phase 4: geometric pairing of structured and
! unstructured coupling interface faces stored in mod_interface's
! Interface_List.
!
! Strategy (conservative, works for non-conformal overlaps):
!   1. Partition Interface_List into structured (PEER_STRUCT) and unstructured
!      (PEER_UNS) faces, keeping their global indices.
!   2. For each unstructured face centroid P_uns:
!        - project P_uns onto the plane of each structured quad face
!        - invert the bilinear map of that quad to get parametric (u,v)
!        - if (u,v) lies in [0,1]^2 the projection is *inside* the quad;
!          pick the inside quad with the smallest normal distance to P_uns
!        - store bilinear weights w = [(1-u)(1-v), u(1-v), uv, (1-u)v]
!          on the uns face (w sums to 1 by construction) and set peer_id
!          on both faces.
!   3. Print a match report: coverage, max/mean projection distance, worst
!      weight-sum deviation.
!
! This is a *geometric* matcher only: it does NOT exchange flow data and does
! NOT convert units (that is phase 5).  It assumes both sides registered into
! the same Interface_List using the same length unit (mm for grid_BC).
!
! Bilinear inversion uses Newton iteration (1 step is exact for a planar
! quad; a few iterations handle mild warping).  Fallback: nearest structured
! face centroid if no quad contains the projection.
!===============================================================================
module mod_interface_match
   use mod_precision, only: dp
   use mod_interface, only: Interface_FACE_TYPE, Interface_List, Num_Interface, &
                            PEER_STRUCT, PEER_UNS, MATCH_MATCHED, MATCH_UNMATCHED, &
                            classify_interface, iface_type_name, iface_recipe_string, &
                            IFACE_UNKNOWN, IFACE_COMP_FLUID_FLUID, &
                            IFACE_COMP_FLUID_POROUS, IFACE_UNS_FLUID_POROUS
   implicit none
   private
   public :: match_interfaces, dispatch_interfaces

   ! module-level state used to pass the selected component indices between
   ! jacobian22 and the residual extraction helpers (avoided passing through
   ! the call chain for brevity).
   integer, save :: ic1 = 1, ic2 = 2

contains

   !--------------------------------------------------------------------------
   ! Main entry point.  tol = maximum allowed normal distance (in the same
   ! length unit as the registered coordinates) for a match to be accepted.
   ! Default 1.0e-3 is ~1 um for the mm-scale grid_BC mesh, comfortably below
   ! the smallest face size while tolerating round-off.
   !--------------------------------------------------------------------------
   subroutine match_interfaces( tol )
      real(dp), optional, intent(in) :: tol

      integer, allocatable :: is_idx(:), iu_idx(:)
      integer :: ns, nu, i, j, best_is
      real(dp) :: dist_tol, P(3), dist, best_d, u, v, w(4)
      real(dp) :: max_d, sum_d, max_werr
      integer  :: n_matched, ns_covered
      logical  :: inside

      dist_tol = 1.0e-3_dp
      if ( present(tol) ) dist_tol = tol

      ! ---- split into structured / unstructured index lists ----
      ns = count( Interface_List(1:Num_Interface)%solver == PEER_STRUCT )
      nu = count( Interface_List(1:Num_Interface)%solver == PEER_UNS   )
      allocate( is_idx(ns), iu_idx(nu) )
      ns = 0; nu = 0
      do i = 1, Num_Interface
         if ( Interface_List(i)%solver == PEER_STRUCT ) then
            ns = ns + 1; is_idx(ns) = i
         else if ( Interface_List(i)%solver == PEER_UNS ) then
            nu = nu + 1; iu_idx(nu) = i
         end if
      end do

      write(*,'(a)') ''
      write(*,'(a)') '--- Interface geometric matching (phase 4) ---'
      write(*,'(a,i0,a,i0,a)') '  structured faces: ', ns, &
                               '   unstructured faces: ', nu
      if ( ns == 0 .or. nu == 0 ) then
         write(*,'(a)') '  WARNING: one side has no interface faces; nothing to match.'
         return
      end if

      ! ---- reset pairing state ----
      do i = 1, Num_Interface
         Interface_List(i)%match_state = MATCH_UNMATCHED
         Interface_List(i)%peer_id     = 0
         if ( allocated(Interface_List(i)%peer_w) ) &
            deallocate( Interface_List(i)%peer_w )
      end do

      ! ---- match each unstructured face to a structured quad ----
      n_matched = 0
      max_d    = 0.0_dp
      sum_d    = 0.0_dp
      max_werr = 0.0_dp

      do j = 1, nu
         P = Interface_List(iu_idx(j))%centroid

         best_is = 0
         best_d  = huge(1.0_dp)
         u = 0.0_dp; v = 0.0_dp

         do i = 1, ns
            call quad_project_and_invert( Interface_List(is_idx(i)), P, &
                                          dist, u, v, inside )
            if ( inside .and. abs(dist) < best_d ) then
               best_d  = abs(dist)
               best_is = i
            end if
         end do

         ! fallback: nearest structured centroid if no quad strictly contains
         ! the projection (e.g. centroid lands on a shared edge and round-off
         ! flips it just outside)
         if ( best_is == 0 ) then
            do i = 1, ns
               dist = norm2( P - Interface_List(is_idx(i))%centroid )
               if ( dist < best_d ) then
                  best_d = dist; best_is = i
               end if
            end do
            call quad_project_and_invert( Interface_List(is_idx(best_is)), P, &
                                          dist, u, v, inside )
         end if

         if ( best_d > dist_tol ) cycle   ! no acceptable match

         ! bilinear weights for the 4 struct quad vertices
         w(1) = (1.0_dp - u) * (1.0_dp - v)
         w(2) = u * (1.0_dp - v)
         w(3) = u * v
         w(4) = (1.0_dp - u) * v

         Interface_List(iu_idx(j))%match_state = MATCH_MATCHED
         Interface_List(iu_idx(j))%peer_id     = is_idx(best_is)
         allocate( Interface_List(iu_idx(j))%peer_w(4) )
         Interface_List(iu_idx(j))%peer_w = w

         ! reverse link: struct face -> uns face.  peer_w on the struct side
         ! is left unallocated; struct->uns interpolation is done by
         ! area-weighted averaging over the uns faces pointing here (phase 5).
         Interface_List(is_idx(best_is))%match_state = MATCH_MATCHED
         Interface_List(is_idx(best_is))%peer_id     = iu_idx(j)

         n_matched = n_matched + 1
         max_d  = max( max_d, best_d )
         sum_d  = sum_d + best_d
         max_werr = max( max_werr, abs( sum(w) - 1.0_dp ) )
      end do

      ! ---- report ----
      write(*,'(a,i0,a,i0,a,f0.1,a)') '  matched: ', n_matched, ' / ', nu, &
            '  uns faces  (', 100.0_dp*real(n_matched,dp)/max(nu,1), '%)'
      ! struct-side coverage: how many struct faces received at least one peer
      ns_covered = count( Interface_List(is_idx)%match_state == MATCH_MATCHED )
      write(*,'(a,i0,a,i0,a,f0.1,a)') '  struct faces referenced: ', ns_covered, &
            ' / ', ns, '  (', 100.0_dp*real(ns_covered,dp)/max(ns,1), '%)'
      if ( n_matched > 0 ) then
         write(*,'(a,es12.4)') '  max normal projection distance : ', max_d
         write(*,'(a,es12.4)') '  mean normal projection distance: ', sum_d / n_matched
         write(*,'(a,es12.4)') '  max |sum(w)-1| deviation        : ', max_werr
      else
         write(*,'(a)') '  WARNING: no faces matched; check overlap / unit consistency.'
      end if
      write(*,'(a)') '--- end matching ---'

   end subroutine match_interfaces

   !---------------------------------------------------------------------------
   ! Phase-11 dispatch: run after match_interfaces.  Classify each matched pair
   ! from the solver + cell-zone class of both sides (classify_interface) and
   ! print the interface dispatch table (type + exchange-quantity recipe per
   ! type).  Unmatched structured faces stay IFACE_UNKNOWN and are counted.
   !---------------------------------------------------------------------------
   subroutine dispatch_interfaces
      integer :: i, j, it, cnt_ff, cnt_fp, cnt_uu, cnt_unk, ns_side

      do i = 1, Num_Interface
         Interface_List(i)%iface_type = IFACE_UNKNOWN
      end do

      cnt_ff = 0; cnt_fp = 0; cnt_uu = 0; cnt_unk = 0; ns_side = 0
      do i = 1, Num_Interface
         if ( Interface_List(i)%solver /= PEER_STRUCT ) cycle
         ns_side = ns_side + 1
         j = Interface_List(i)%peer_id
         if ( Interface_List(i)%match_state /= MATCH_MATCHED .or. &
              j < 1 .or. j > Num_Interface ) then
            cnt_unk = cnt_unk + 1
            cycle
         end if
         it = classify_interface( Interface_List(i)%solver, &
                                  Interface_List(i)%cz_type, &
                                  Interface_List(j)%solver, &
                                  Interface_List(j)%cz_type )
         Interface_List(i)%iface_type = it
         Interface_List(j)%iface_type = it
         select case ( it )
         case ( IFACE_COMP_FLUID_FLUID );  cnt_ff = cnt_ff + 1
         case ( IFACE_COMP_FLUID_POROUS ); cnt_fp = cnt_fp + 1
         case ( IFACE_UNS_FLUID_POROUS );  cnt_uu = cnt_uu + 1
         case default;                     cnt_unk = cnt_unk + 1
         end select
      end do

      write(*,'(a)') ''
      write(*,'(a)') '--- Interface dispatch table (phase 11) ---'
      write(*,'(a,i0)') '  struct faces           : ', ns_side
      write(*,'(a,i0)') '  unclassified (no match): ', cnt_unk
      if ( cnt_ff > 0 ) then
         write(*,'(a,i0,a,a)') '  faces ', cnt_ff, '  ', &
            trim(iface_type_name(IFACE_COMP_FLUID_FLUID))
         write(*,'(a,a)')      '        recipe:', &
            trim(iface_recipe_string(IFACE_COMP_FLUID_FLUID))
      end if
      if ( cnt_fp > 0 ) then
         write(*,'(a,i0,a,a)') '  faces ', cnt_fp, '  ', &
            trim(iface_type_name(IFACE_COMP_FLUID_POROUS))
         write(*,'(a,a)')      '        recipe:', &
            trim(iface_recipe_string(IFACE_COMP_FLUID_POROUS))
      end if
      if ( cnt_uu > 0 ) then
         write(*,'(a,i0,a,a)') '  faces ', cnt_uu, '  ', &
            trim(iface_type_name(IFACE_UNS_FLUID_POROUS))
         write(*,'(a,a)')      '        recipe:', &
            trim(iface_recipe_string(IFACE_UNS_FLUID_POROUS))
      end if
      write(*,'(a)') '--- end dispatch ---'
   end subroutine dispatch_interfaces

   !--------------------------------------------------------------------------
   ! Project point P onto the plane of a quad face, invert the bilinear map
   ! to get parametric (u,v), and report whether the projection lands inside
   ! the quad (u,v in [0,1]).  dist = signed normal distance from P to plane.
   !
   ! Quad vertices (verts(:,1:4)) are A,B,C,D ordered CCW around the face
   ! normal.  Bilinear surface:
   !   S(u,v) = (1-u)(1-v) A + u(1-v) B + u v C + (1-u) v D
   ! Newton solve S(u,v) = Pproj using the two in-plane directions of largest
   ! span to keep the 2x2 Jacobian well-conditioned.
   !--------------------------------------------------------------------------
   subroutine quad_project_and_invert( F, P, dist, u, v, inside )
      type(Interface_FACE_TYPE), intent(in) :: F
      real(dp), intent(in)  :: P(3)
      real(dp), intent(out) :: dist, u, v
      logical,  intent(out) :: inside

      real(dp) :: VA(3), VB(3), VC(3), VD(3), nrm(3), Pproj(3)
      real(dp) :: Su(3), Sv(3), R(3), J(2,2), Jinv(2,2), det, du(2)
      real(dp) :: f1, f2, tol_local
      integer :: iter

      inside = .false.
      u = 0.5_dp
      v = 0.5_dp

      if ( F%nv < 3 ) return
      VA = F%verts(:,1)
      VB = F%verts(:,2)
      VC = F%verts(:,3)
      if ( F%nv >= 4 ) then
         VD = F%verts(:,4)
      else
         VD = VC   ! triangle: degenerate quad, weight on D is zero
      end if

      ! plane normal (from the first triangle)
      nrm = cross( VB - VA, VC - VA )
      if ( norm2(nrm) < 1.0e-30_dp ) return
      nrm = nrm / norm2(nrm)

      ! signed distance and projection onto the plane
      dist = dot_product( P - VA, nrm )
      Pproj = P - dist * nrm

      ! Newton iteration on the two in-plane components of R = S(u,v) - Pproj
      tol_local = 1.0e-12_dp
      do iter = 1, 20
         Su = (1.0_dp - v)*(VB - VA) + v*(VC - VD)
         Sv = (1.0_dp - u)*(VD - VA) + u*(VC - VB)
         R  = bilinear_eval(VA,VB,VC,VD,u,v) - Pproj

         call jacobian22( Su, Sv, J )
         det = J(1,1)*J(2,2) - J(1,2)*J(2,1)
         if ( abs(det) < 1.0e-30_dp ) exit
         Jinv(1,1) =  J(2,2)/det
         Jinv(1,2) = -J(1,2)/det
         Jinv(2,1) = -J(2,1)/det
         Jinv(2,2) =  J(1,1)/det

         f1 = residual_comp1( R )
         f2 = residual_comp2( R )
         du(1) = -( Jinv(1,1)*f1 + Jinv(1,2)*f2 )
         du(2) = -( Jinv(2,1)*f1 + Jinv(2,2)*f2 )

         u = u + du(1)
         v = v + du(2)

         if ( norm2(du) < tol_local ) exit
      end do

      ! inside if (u,v) within a small epsilon of [0,1]^2
      inside = ( u >= -1.0e-9_dp .and. u <= 1.0_dp + 1.0e-9_dp .and. &
                 v >= -1.0e-9_dp .and. v <= 1.0_dp + 1.0e-9_dp )
      ! clamp to [0,1] for weight computation
      u = max(0.0_dp, min(1.0_dp, u))
      v = max(0.0_dp, min(1.0_dp, v))

   end subroutine quad_project_and_invert

   pure function bilinear_eval( A, B, C, D, u, v ) result( S )
      real(dp), intent(in) :: A(3), B(3), C(3), D(3), u, v
      real(dp) :: S(3)
      S = (1-u)*(1-v)*A + u*(1-v)*B + u*v*C + (1-u)*v*D
   end function bilinear_eval

   pure function cross( a, b ) result( c )
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: c(3)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross

   ! select the two components of (Su,Sv) with the largest combined magnitude
   ! and build the 2x2 Jacobian from them.  Stores the chosen indices in the
   ! module variables ic1/ic2 for residual_comp1/2.
   subroutine jacobian22( Su, Sv, J )
      real(dp), intent(in)  :: Su(3), Sv(3)
      real(dp), intent(out) :: J(2,2)
      real(dp) :: mag(3)
      mag = abs(Su) + abs(Sv)
      ic1 = maxloc(mag, 1)
      mag(ic1) = -1.0_dp
      ic2 = maxloc(mag, 1)
      J(1,1) = Su(ic1); J(1,2) = Sv(ic1)
      J(2,1) = Su(ic2); J(2,2) = Sv(ic2)
   end subroutine jacobian22

   real(dp) function residual_comp1( R )
      real(dp), intent(in) :: R(3)
      residual_comp1 = R(ic1)
   end function residual_comp1

   real(dp) function residual_comp2( R )
      real(dp), intent(in) :: R(3)
      residual_comp2 = R(ic2)
   end function residual_comp2

end module mod_interface_match
