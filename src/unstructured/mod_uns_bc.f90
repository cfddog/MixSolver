!===============================================================================
! mod_bc.f90 -- Boundary-condition groups on face zones (phase 3)
!
! Each control bc line becomes one boundary group holding the faces of its
! zone. A wall group may carry a lid plane filter: faces whose centroid
! coordinate along lid_dir matches lid_coord move with lid_vel, the rest are
! stationary.
!===============================================================================
module mod_uns_bc
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_control
   implicit none
   private
   public :: bc_t, build_bc, bc_face_vel, bc_face_p, bc_face_T, bc_type_name, &
             set_interface_vel, set_interface_p, set_interface_T, &
             set_inlet_ramp_factor, set_pval_ramp_factor, &
             bc_face_p_dirichlet, set_iface_p_dirichlet, iface_p_is_dirichlet, &
             set_interface_slip, set_iface_slip_mode, iface_slip_is_on, &
             set_interface_T_ltne, set_iface_t_ltne_mode, iface_t_ltne_is_on, &
             set_iface_zhang_vel, iface_zhang_vel_is_on

   ! Module-level ramp factor for mass-flow-inlet faces.
   ! simple_run/simple_run_mpi set it each outer iteration to
   !   min(1, iter / max(inlet_ramp,1))
   ! bc_face_vel multiplies the MASSINLET uspeed by this factor.
   real(dp), save :: g_inlet_ramp_factor = 1.0_dp

   ! Module-level ramp factor for pressure-Dirichlet (POUTLET / FARFIELD)
   ! faces.  During a cold start the mass-flow inlet is ramped to a tiny
   ! velocity, so the convective scale rho*u^2 is tiny; a fixed outlet
   ! pressure applied at full strength is then enormous by comparison and
   ! blows up the first pressure correction.  simple_run/simple_run_mpi ramp
   ! it in lock-step with the inlet (and it is 1 once the ramp is complete,
   ! so the converged solution is unaffected).
   real(dp), save :: g_pval_ramp_factor = 1.0_dp

   ! ---- coupling-interface closure switches (phase 14) ---------------------
   ! g_iface_p_dir : the interface pressure is a per-face DIRICHLET value
   !                 (bcs%iface_p) instead of zero-gradient.  Set by the
   !                 coupling driver when mix.control has iface_p_model =
   !                 dirichlet|momentum.  When .true. the PPE treats those
   !                 faces exactly like POUTLET/FARFIELD (af on the owner
   !                 diagonal, p'_f = 0, flux correction) -> the otherwise
   !                 pure-Neumann slab PPE becomes well posed, cell 1 is not
   !                 pinned and the uniform interface-flux correction is not
   !                 needed.  Default .false. = legacy zero-gradient.
   ! g_iface_slip  : tangential stress-jump (Beavers-Joseph) friction at the
   !                 interface faces; bcs%iface_slip_C holds the per-face
   !                 conductance [kg/s] and bcs%iface_peer_vel the peer face
   !                 velocity used as the free-fluid tangential target.
   logical,  save :: g_iface_p_dir = .false.
   logical,  save :: g_iface_slip  = .false.
   ! g_iface_zhang_vel : the interface velocity comes from the Zhang 2011
   !                 Eqs.22/25/26 law (opt-in, iface_velocity = zhang), i.e.
   !                 it is a VISCOUS closure value -- not an advective inflow.
   !                 The momentum assembly then treats the interface face as a
   !                 pure viscous Dirichlet: the prescribed face velocity
   !                 enters only through the diffusive coefficient D, while the
   !                 mass-flux (upwind) term carries the INTERIOR value on
   !                 outflow, i.e. rhs += D*uf - min(F,0)*uf instead of
   !                 D*uf - F*uf.  Without this, the (physically O(10 m/s))
   !                 tangential slip of Eq.26 is convected into the first slab
   !                 cell with the full mass-flux weight F >> D and drags the
   !                 porous slab toward the free stream (measured |u|max ->
   !                 100 m/s on C_P_test).  Default .false. keeps the legacy
   !                 path bit-identical.
   logical,  save :: g_iface_zhang_vel = .false.
   ! g_iface_T_ltne : with the LTNE thermal model the interface heat flux must
   !                 be split between the porous phases by area ratio
   !                 (porosity): the fluid-phase equation sees bcs%iface_T and
   !                 the solid-phase equation bcs%iface_Ts, both supplied per
   !                 face by the coupling driver from the Zhang 2011 Eq.27-29
   !                 closure.  Default .false. = one shared Dirichlet value
   !                 (legacy).
   logical,  save :: g_iface_T_ltne = .false.

   ! one boundary group (a set of boundary faces with common condition)
   type :: bcgroup_t
      integer              :: zone  = 0
      integer              :: btype = BC_NONE
      real(dp)             :: uvel(3)  = 0.0_dp
      real(dp)             :: pval     = 0.0_dp
      real(dp)             :: mdot     = 0.0_dp   ! mass-flux inlet kg/m^2/s
      real(dp)             :: uspeed   = 0.0_dp   ! inlet speed mdot/rho (m/s)
      ! --- parabolic (fully-developed) velocity inlet (BC_VINLET_PARAB) ---
      real(dp)             :: prfl_mean   = 0.0_dp ! profile mean speed Umean
      real(dp)             :: prfl_span   = 1.0_dp ! span over xi = 0..1
      real(dp)             :: prfl_origin = 0.0_dp ! coordinate where |u| = 0
      integer              :: prfl_axis   = 2      ! 1=x, 2=y, 3=z
      logical              :: has_lid  = .false.
      integer              :: lid_dir  = 0
      real(dp)             :: lid_coord = 0.0_dp
      real(dp)             :: lid_tol   = 0.0_dp
      real(dp)             :: lid_vel(3) = 0.0_dp
      ! --- thermal boundary conditions (phase 10) ---
      integer              :: ttype    = 0          ! 0=adiabatic,1=fixed-T,2=fixed-flux
      real(dp)             :: tval     = 0.0_dp     ! fixed-temperature value
      real(dp)             :: qval     = 0.0_dp     ! fixed-heat-flux value
      integer              :: ntbc_plane = 0        ! number of plane filters
      integer              :: tbc_pdir(4)   = 0      ! 1=x, 2=y, 3=z
      real(dp)             :: tbc_pcoord(4) = 0.0_dp
      real(dp)             :: tbc_tol   = 0.0_dp    ! plane match tolerance
      integer              :: tbc_pttype(4) = 0
      real(dp)             :: tbc_ptval(4) = 0.0_dp
      real(dp)             :: tbc_pqval(4) = 0.0_dp
      integer              :: nf = 0
      integer, allocatable :: faces(:)     ! boundary face indices
   end type bcgroup_t

   ! boundary-condition container
   type :: bc_t
      integer              :: nb = 0
      type(bcgroup_t), allocatable :: gb(:)
      integer, allocatable :: fgrp(:)      ! face -> group (0 = interior)
      ! ---- coupling-interface face data (phase 6) ----
      ! Per-face velocity / pressure supplied by the peer solver via the
      ! exchange layer.  Only faces whose group has btype == BC_INTERFACE
      ! read these; others are ignored.  Sized to (3, nfaces) / (nfaces).
      real(dp), allocatable :: iface_vel(:,:)
      real(dp), allocatable :: iface_p(:)
      ! per-face temperature supplied by the peer solver (Dirichlet at the
      ! velocity-specified coupling interface)
      real(dp), allocatable :: iface_T(:)
      ! per-face SOLID-PHASE interface temperature (phase 14, iface_t_model =
      ! zhang + LTNE): the complement of the porosity-weighted split, Eq.27-29.
      real(dp), allocatable :: iface_Ts(:)
      ! ---- interface tangential stress jump (phase 14, iface_slip = bj) ----
      real(dp), allocatable :: iface_slip_C(:)     ! (nfaces) conductance [kg/s]
      real(dp), allocatable :: iface_peer_vel(:,:) ! (3,nfaces) peer face vel
   end type bc_t

contains

   !----------------------------------------------------------------------------
   ! Build boundary groups from the control specs and the mesh zones.
   ! allow_empty = .true. (MPI local meshes): a bc zone with no boundary face
   ! on THIS rank is legal (its faces live on other ranks) -- build an empty
   ! group instead of erroring.  Serial / global meshes keep the strict check.
   !----------------------------------------------------------------------------
   subroutine build_bc( m, ctrl, bcs, ier, allow_empty )
      type(mesh_t),  intent(in)  :: m
      type(ctrl_t),  intent(in)  :: ctrl
      type(bc_t),    intent(out) :: bcs
      integer,       intent(out) :: ier
      logical, intent(in), optional :: allow_empty

      integer :: i, k, nz, iz, nbnd
      real(dp) :: coord_range(3)
      logical  :: allow

      allow = .false.
      if ( present(allow_empty) ) allow = allow_empty

      ier = 0
      allocate( bcs%gb(max(ctrl%nbc,1)) )
      bcs%nb = ctrl%nbc
      allocate( bcs%fgrp(m%nfaces) )
      bcs%fgrp = 0

      ! coordinate range for the lid matching tolerance
      coord_range = maxval( m%x, dim=2 ) - minval( m%x, dim=2 )

      ! ---- collect faces per group ------------------------------------------
      do i = 1, ctrl%nbc
         associate( gb => bcs%gb(i), sp => ctrl%bc(i) )
            gb%zone   = sp%zone
            gb%btype  = sp%btype
            gb%uvel   = sp%uvel
            gb%pval   = sp%pval
            gb%mdot   = sp%mdot
            ! pure-polynomial (fully-developed) inlet profile parameters
            gb%prfl_mean   = sp%umean
            gb%prfl_span   = sp%uspan
            gb%prfl_origin = sp%uorigin
            gb%prfl_axis   = sp%uaxis
            if ( sp%btype == BC_VINLET_PARAB ) then
               if ( gb%prfl_span <= 0.0_dp .or. gb%prfl_axis < 1 .or. &
                    gb%prfl_axis > 3 ) then
                  write(*,'(a)') 'ERROR: velocity-inlet-parabolic needs ' // &
                     'span > 0 and axis 1|2|3'
                  ier = 14
                  return
               end if
            end if
            ! mass-flux inlet: convert to an inlet speed magnitude; the
            ! per-face direction is taken from the outward face normal in
            ! bc_face_vel (u_f = -uspeed * n_outward).
            if ( sp%btype == BC_MASSINLET ) then
               if ( ctrl%rho <= 0.0_dp ) then
                  write(*,'(a)') 'ERROR: mass-flow-inlet needs rho > 0'
                  ier = 13
                  return
               end if
               gb%uspeed = sp%mdot / ctrl%rho
            end if
            gb%has_lid = sp%has_lid
            gb%lid_dir = sp%lid_dir
            gb%lid_coord = sp%lid_coord
            gb%lid_vel  = sp%lid_vel
            if ( gb%has_lid ) &
               gb%lid_tol = 1.0e-4_dp * maxval( coord_range )
            ! --- thermal boundary (phase 10) ---
            gb%ttype    = sp%ttype
            gb%tval     = sp%tval
            gb%qval     = sp%qval
            gb%ntbc_plane = sp%ntbc_plane
            gb%tbc_pdir  = sp%tbc_pdir
            gb%tbc_pcoord = sp%tbc_pcoord
            gb%tbc_pttype = sp%tbc_pttype
            gb%tbc_ptval = sp%tbc_ptval
            gb%tbc_pqval = sp%tbc_pqval
            if ( gb%ntbc_plane > 0 ) &
               gb%tbc_tol = 1.0e-4_dp * maxval( coord_range )

            ! zone id must exist
            nz = 0
            do iz = 1, m%nzone
               if ( m%zone(iz)%id == sp%zone ) nz = iz
            end do
            if ( nz == 0 ) then
               write(*,'(a,i0)') 'ERROR: bc refers to unknown zone ', sp%zone
               ier = 10
               return
            end if

            ! count boundary faces of this zone
            nbnd = 0
            do k = 1, m%nfaces
               if ( m%f(k)%c1 == 0 .and. m%f(k)%zone == sp%zone ) nbnd = nbnd+1
            end do
            if ( nbnd == 0 ) then
               if ( .not. allow ) then
                  write(*,'(a,i0)') 'ERROR: bc zone has no boundary faces: ', &
                                    sp%zone
                  ier = 11
                  return
               end if
            end if

            allocate( gb%faces(nbnd) )
            gb%nf = 0
            do k = 1, m%nfaces
               if ( m%f(k)%c1 == 0 .and. m%f(k)%zone == sp%zone ) then
                  gb%nf = gb%nf + 1
                  gb%faces(gb%nf) = k
                  bcs%fgrp(k) = i
               end if
            end do
         end associate
      end do

      ! ---- every boundary face must belong to a group ------------------------
      ! Exception: faces of a CAS "interface" zone are coupling surfaces
      ! handled by the exchange layer (phase 5), not physical BCs.  They are
      ! registered in mod_interface by register_interface_zones; here they are
      ! simply left without a group and reported once per zone.
      do k = 1, m%nfaces
         if ( m%f(k)%c1 == 0 .and. bcs%fgrp(k) == 0 ) then
            if ( zone_is_interface( m, m%f(k)%zone ) ) cycle
            write(*,'(a,i0,a,i0)') 'ERROR: boundary face ', k, ' (zone ', &
                                   m%f(k)%zone, ') has no bc specification'
            ier = 12
            return
         end if
      end do

      ! ---- coupling-interface faces get a BC_INTERFACE group (phase 6) ----
      ! Faces of CAS "interface" zones are coupling surfaces.  They used to be
      ! left without a group (fgrp=0) and skipped by the flux assembly.  Now
      ! they are collected into a single BC_INTERFACE group so that the face
      ! velocity/pressure can be supplied by the peer solver via iface_vel /
      ! iface_p (set by set_interface_vel / set_interface_p before each step).
      nbnd = 0
      do k = 1, m%nfaces
         if ( m%f(k)%c1 == 0 .and. bcs%fgrp(k) == 0 .and. &
              zone_is_interface( m, m%f(k)%zone ) ) nbnd = nbnd + 1
      end do
      if ( nbnd > 0 ) then
         bcs%nb = bcs%nb + 1
         if ( bcs%nb > size(bcs%gb) ) then
            ! gb was allocated to max(ctrl%nbc,1); grow if interface group
            ! pushes past that (rare: ctrl has no bc at all).
            block
               type(bcgroup_t), allocatable :: tmp(:)
               allocate( tmp(bcs%nb) )
               tmp(1:bcs%nb-1) = bcs%gb(1:bcs%nb-1)
               call move_alloc( tmp, bcs%gb )
            end block
         end if
         associate( gb => bcs%gb(bcs%nb) )
            gb%zone  = 0
            gb%btype = BC_INTERFACE
            gb%nf    = nbnd
            allocate( gb%faces(nbnd) )
            nbnd = 0
            do k = 1, m%nfaces
               if ( m%f(k)%c1 == 0 .and. bcs%fgrp(k) == 0 .and. &
                    zone_is_interface( m, m%f(k)%zone ) ) then
                  nbnd = nbnd + 1
                  gb%faces(nbnd) = k
                  bcs%fgrp(k) = bcs%nb
               end if
            end do
         end associate
         write(*,'(a,i0,a)') '  coupling interface: ', nbnd, &
                              ' faces -> BC_INTERFACE group'
      end if

      ! allocate per-face interface data arrays (zero-initialised; the driver
      ! fills them via set_interface_vel / set_interface_p each coupling step)
      allocate( bcs%iface_vel(3, m%nfaces), bcs%iface_p(m%nfaces), &
                bcs%iface_T(m%nfaces), bcs%iface_slip_C(m%nfaces), &
                bcs%iface_peer_vel(3, m%nfaces), bcs%iface_Ts(m%nfaces) )
      bcs%iface_vel = 0.0_dp
      bcs%iface_p   = 0.0_dp
      bcs%iface_T   = 0.0_dp
      bcs%iface_slip_C   = 0.0_dp
      bcs%iface_peer_vel = 0.0_dp
      bcs%iface_Ts       = 0.0_dp

      do i = 1, m%nzone
         if ( zone_is_interface( m, m%zone(i)%id ) ) &
            write(*,'(a,i0,a,a,a)') '  zone', m%zone(i)%id, '  ', &
               trim(m%zone(i)%cond_name), &
               '  -> coupling interface'
      end do

      ! ---- report -------------------------------------------------------------
      write(*,'(a)') ''
      write(*,'(a)') '--- Boundary conditions ---'
      do i = 1, bcs%nb
         associate( gb => bcs%gb(i) )
            write(*,'(a,i6,a,a,a,i8)') '  zone', gb%zone, '  ', &
               bc_type_name( gb%btype ), '  faces: ', gb%nf
            if ( gb%btype == BC_VINLET ) &
               write(*,'(a,3(es10.3,1x))') '    velocity        : ', gb%uvel
            if ( gb%btype == BC_VINLET_PARAB ) then
               write(*,'(a,es10.3,a,es10.3)') '    mean speed Umean: ', &
                  gb%prfl_mean, '   span H = ', gb%prfl_span
               write(*,'(a,i0,a,es10.3,a)') '    profile axis    : ', &
                  gb%prfl_axis, '   vanishing at ', gb%prfl_origin, &
                  '  (u = 6*Umean*xi*(1-xi), normal inflow)'
            end if
            if ( gb%btype == BC_MASSINLET ) then
               write(*,'(a,es12.4,a)') '    mass flux m''   : ', gb%mdot, ' kg/m^2/s'
               write(*,'(a,es10.3,a)') '    inlet speed     : ', gb%uspeed, ' m/s (normal)'
            end if
            if ( gb%btype == BC_POUTLET ) &
               write(*,'(a,es10.3)') '    pressure        : ', gb%pval
            if ( gb%btype == BC_FARFIELD ) then
               write(*,'(a,es10.3)') '    pressure (far)  : ', gb%pval
               write(*,'(a,3(es10.3,1x))') '    velocity (far) : ', gb%uvel
            end if
            if ( gb%has_lid ) &
               write(*,'(a,i0,a,es10.3,a,3(es10.3,1x))') &
                  '    lid on plane ', gb%lid_dir, ' at', gb%lid_coord, &
                  '  u = ', gb%lid_vel
            if ( gb%ttype == 1 ) &
               write(*,'(a,es10.3)') '    fixed-T wall   : ', gb%tval
            if ( gb%ttype == 2 ) &
               write(*,'(a,es10.3)') '    fixed-flux wall: ', gb%qval
            do k = 1, gb%ntbc_plane
               if ( gb%tbc_pttype(k) == 1 ) &
                  write(*,'(a,i0,a,es10.3,a,es10.3)') &
                     '    tbc_plane dir ', gb%tbc_pdir(k), ' at', &
                     gb%tbc_pcoord(k), '  fixed-T = ', gb%tbc_ptval(k)
               if ( gb%tbc_pttype(k) == 2 ) &
                  write(*,'(a,i0,a,es10.3,a,es10.3)') &
                     '    tbc_plane dir ', gb%tbc_pdir(k), ' at', &
                     gb%tbc_pcoord(k), '  fixed-flux = ', gb%tbc_pqval(k)
            end do
         end associate
      end do

   end subroutine build_bc

   !----------------------------------------------------------------------------
   ! Face velocity at boundary face i (for convection/diffusion evaluation).
   ! uP is the velocity of the owner cell (used for symmetry mirroring and
   ! outlet extrapolation).
   !----------------------------------------------------------------------------
   subroutine bc_face_vel( bcs, i, xf, sf, uP, uf )
      type(bc_t),   intent(in)  :: bcs
      integer,      intent(in)  :: i
      real(dp),     intent(in)  :: xf(3), sf(3), uP(3)
      real(dp),     intent(out) :: uf(3)

      integer  :: g
      real(dp) :: un, nvec(3), xi, uspc

      g = bcs%fgrp(i)

      select case ( bcs%gb(g)%btype )
      case ( BC_WALL )
         if ( bcs%gb(g)%has_lid .and. &
              abs( xf(bcs%gb(g)%lid_dir) - bcs%gb(g)%lid_coord ) &
                  <= bcs%gb(g)%lid_tol ) then
            uf = bcs%gb(g)%lid_vel
         else
            uf = 0.0_dp
         end if

      case ( BC_SYMMETRY )
         ! specular reflection: normal component zeroed, tangential kept
         nvec = sf / norm2( sf )
         un   = dot_product( uP, nvec )
         uf   = uP - 2.0_dp * un * nvec

      case ( BC_VINLET )
         uf = bcs%gb(g)%uvel

      case ( BC_VINLET_PARAB )
         ! Fully-developed (parabolic) inlet: profile along the INWARD normal.
         !   |u_f| = 6 * Umean * xi * (1 - xi),  xi = (x_axis - origin)/span
         ! xi is clamped to [0,1] (a patch lying slightly outside the nominal
         ! channel keeps |u_f| >= 0 instead of going negative); over the span
         ! the profile averages exactly Umean, so the volume flow equals a
         ! uniform 'velocity-inlet Umean' on the same patch.
         nvec = sf / norm2( sf )
         xi   = ( xf(bcs%gb(g)%prfl_axis) - bcs%gb(g)%prfl_origin ) &
                / bcs%gb(g)%prfl_span
         xi   = max( 0.0_dp, min( 1.0_dp, xi ) )
         uspc = 6.0_dp * bcs%gb(g)%prfl_mean * xi * ( 1.0_dp - xi )
         uf   = -uspc * nvec

      case ( BC_MASSINLET )
         ! Uniform-normal mass-flux inlet.  sf points outward from the fluid
         ! domain, so the inflow velocity points inward:
         !    u_f = -(mdot/rho) * n_outward * ramp_factor
         ! Tangential components are zero (normal-only inflow).
         nvec = sf / norm2( sf )
         uf   = -bcs%gb(g)%uspeed * g_inlet_ramp_factor * nvec

      case ( BC_POUTLET )
         uf = uP                        ! zeroth-order extrapolation

      case ( BC_OUTFLOW )
         uf = uP                        ! zeroth-order extrapolation (zero
                                        ! normal gradient; flux rescaled
                                        ! globally in flux_rhiechow path)

      case ( BC_INTERFACE )
         ! Coupling interface.  Two regimes (phase 14):
         !   iface_p_model = grad0 (legacy): the peer supplies the whole face
         !     VELOCITY vector through the exchange layer (bcs%iface_vel), i.e.
         !     this is a velocity-Dirichlet patch.
         !   iface_p_model = dirichlet|momentum: the peer supplies the face
         !     PRESSURE (bcs%iface_p, see bc_face_p), so the velocity here must
         !     be EXTRAPOLATED from the owner cell exactly like pressure-outlet
         !     -- deriving it from the stale exchanged vector instead makes the
         !     PPE face flux and the momentum upwind term inconsistent with the
         !     field (measured: divergent interface flux in C_P_test).
         if ( g_iface_p_dir ) then
            uf = uP
         else
            uf = bcs%iface_vel(:,i)
         end if

      case ( BC_FARFIELD )
         ! Characteristic far-field convention (sf points outward from the
         ! fluid domain).  The in/out-flow decision is made with the
         ! prescribed EXTERIOR state u_far, not the interior state uP:
         !   inflow  (u_far.n < 0): impose the free-stream velocity
         !   outflow (u_far.n > 0): zeroth-order extrapolation from owner cell
         ! NB: testing uP.n here fails when the field is started from rest:
         ! uP=0 gives uP.n=0 on every far-field face, the upstream boundary
         ! is misclassified as outflow (uf = uP = 0), no mass enters the
         ! domain, and the zero solution is stationary ("converges" at it=1).
         nvec = sf / norm2( sf )
         un   = dot_product( bcs%gb(g)%uvel, nvec )
         if ( un < 0.0_dp ) then
            uf = bcs%gb(g)%uvel
         else
            uf = uP
         end if

      case ( BC_SLIPWALL )
         ! Slip wall: zero normal velocity, tangential preserved.
         !   uf = uP - (uP.n) * n
         ! This gives uf.n = 0 (no mass flux through wall) and tangential
         ! velocity equal to the cell-center tangential component. Correct
         ! for inviscid flow where no-slip would create a singular BL.
         nvec = sf / norm2( sf )
         un   = dot_product( uP, nvec )
         uf   = uP - un * nvec

      case default
         uf = uP
      end select

   end subroutine bc_face_vel

   !----------------------------------------------------------------------------
   ! Face pressure at boundary face i; pP is the owner cell pressure
   !----------------------------------------------------------------------------
   pure subroutine bc_face_p( bcs, i, pP, pf )
      type(bc_t),   intent(in)  :: bcs
      integer,      intent(in)  :: i
      real(dp),     intent(in)  :: pP
      real(dp),     intent(out) :: pf

      integer :: g

      g = bcs%fgrp(i)
      if ( bcs%gb(g)%btype == BC_POUTLET .or. &
           bcs%gb(g)%btype == BC_FARFIELD ) then
         pf = bcs%gb(g)%pval * g_pval_ramp_factor   ! Dirichlet free-stream p
                                                    ! (ramped on cold start)
      else if ( bcs%gb(g)%btype == BC_INTERFACE ) then
         if ( g_iface_p_dir ) then
            ! phase 14 (iface_p_model = dirichlet|momentum): the peer solver
            ! supplies a per-face DIRICHLET interface pressure (bcs%iface_p,
            ! gauge Pa).  The interface VELOCITY is no longer Dirichlet on
            ! these faces -- the momentum assembly switches to the
            ! zero-gradient/outflow treatment and the interface mass flux is
            ! solved from continuity (closing the partition the literature way:
            ! Betchen 2006 Eq.16 + Eq.43/44, pressure is the interface
            ! constraint on the porous side).
            pf = bcs%iface_p(i)
         else
            ! Legacy: the peer supplies the face VELOCITY on this patch, so
            ! pressure must be zero-gradient here (imposing both over-constrains
            ! SIMPLE: p pinned at the inlet and the interface, the boundary
            ! momentum cannot balance and the inlet velocity explodes).
            pf = pP
         end if
      else
         pf = pP                        ! zero normal gradient
      end if

   end subroutine bc_face_p

   !----------------------------------------------------------------------------
   ! Face temperature boundary condition (phase 10).
   ! T_cell is the owner-cell temperature.  For a Dirichlet face (fixed-T wall
   ! or inlet) is_neumann=.false. and T_face holds the prescribed value.
   ! For a Neumann face (adiabatic / fixed-flux / symmetry / outlet)
   ! is_neumann=.true., q_face holds the prescribed heat flux (positive =
   ! heating into the domain), and T_face is set to T_cell so that the
   ! Green-Gauss gradient gives a zero normal gradient contribution.
   !----------------------------------------------------------------------------
   subroutine bc_face_T( bcs, i, xf, T_cell, T_face, q_face, is_neumann, phase )
      type(bc_t),   intent(in)  :: bcs
      integer,      intent(in)  :: i
      real(dp),     intent(in)  :: xf(3), T_cell
      real(dp),     intent(out) :: T_face, q_face
      logical,      intent(out) :: is_neumann
      ! phase: which porous phase the caller is assembling (0/absent = the
      ! legacy single-temperature treatment; 2 = the solid-phase equation).
      ! Only the coupling interface differentiates: with iface_t_model = zhang
      ! under LTNE the fluid phase sees bcs%iface_T and the solid phase
      ! bcs%iface_Ts (Zhang 2011 Eq.27-29, porosity-weighted flux split).
      integer,      intent(in), optional :: phase

      integer  :: g, tt, ipl
      real(dp) :: tv, qv

      g = bcs%fgrp(i)
      T_face  = T_cell
      q_face  = 0.0_dp
      is_neumann = .true.

      select case ( bcs%gb(g)%btype )
      case ( BC_WALL )
         ! zone-level default; plane filters override on matching faces
         tt = bcs%gb(g)%ttype
         tv = bcs%gb(g)%tval
         qv = bcs%gb(g)%qval
         do ipl = 1, bcs%gb(g)%ntbc_plane
            if ( abs( xf(bcs%gb(g)%tbc_pdir(ipl)) - bcs%gb(g)%tbc_pcoord(ipl) ) &
                 <= bcs%gb(g)%tbc_tol ) then
               tt = bcs%gb(g)%tbc_pttype(ipl)
               tv = bcs%gb(g)%tbc_ptval(ipl)
               qv = bcs%gb(g)%tbc_pqval(ipl)
               exit
            end if
         end do

         select case ( tt )
         case ( 0 )                           ! adiabatic
            is_neumann = .true.
            q_face     = 0.0_dp
            T_face     = T_cell
         case ( 1 )                           ! fixed-temperature (Dirichlet)
            is_neumann = .false.
            T_face     = tv
            q_face     = 0.0_dp
         case ( 2 )                           ! fixed-heat-flux (Neumann)
            is_neumann = .true.
            q_face     = qv
            T_face     = T_cell
         end select

      case ( BC_SYMMETRY )
         is_neumann = .true.
         q_face     = 0.0_dp
         T_face     = T_cell

      case ( BC_VINLET, BC_VINLET_PARAB )
         ! inlet fixed temperature (Dirichlet); tval defaults to 0
         is_neumann = .false.
         T_face     = bcs%gb(g)%tval
         q_face     = 0.0_dp

      case ( BC_MASSINLET )
         ! same thermal treatment as a velocity inlet (Dirichlet tval;
         ! adiabatic by default when no tbc is given)
         is_neumann = .false.
         T_face     = bcs%gb(g)%tval
         q_face     = 0.0_dp

      case ( BC_POUTLET, BC_OUTFLOW )
         ! zero normal gradient (convective outflow) by default; a zone-level
         ! tbc (or plane filter) with ttype = 2 injects a fixed heat flux
         ! through the outlet face (e.g. heated end face of a porous sample)
         tt = bcs%gb(g)%ttype
         tv = bcs%gb(g)%tval
         qv = bcs%gb(g)%qval
         do ipl = 1, bcs%gb(g)%ntbc_plane
            if ( abs( xf(bcs%gb(g)%tbc_pdir(ipl)) - bcs%gb(g)%tbc_pcoord(ipl) ) &
                 <= bcs%gb(g)%tbc_tol ) then
               tt = bcs%gb(g)%tbc_pttype(ipl)
               tv = bcs%gb(g)%tbc_ptval(ipl)
               qv = bcs%gb(g)%tbc_pqval(ipl)
               exit
            end if
         end do
         select case ( tt )
         case ( 2 )                           ! fixed-heat-flux (Neumann)
            is_neumann = .true.
            q_face     = qv
            T_face     = T_cell
         case default                         ! convective outflow
            is_neumann = .true.
            q_face     = 0.0_dp
            T_face     = T_cell
         end select

      case ( BC_FARFIELD )
         ! Simplified: always Neumann (zero gradient). For full characteristic
         ! treatment one would Dirichlet-fix T to a free-stream value at inflow
         ! and extrapolate at outflow; the cylinder inviscid case has no energy
         ! equation so this suffices.
         is_neumann = .true.
         q_face     = 0.0_dp
         T_face     = T_cell

      case ( BC_INTERFACE )
         ! The peer solver supplies the inflow temperature (velocity-specified
         ! coupling interface).  Dirichlet is mandatory: with a zero-gradient
         ! inflow the scalar has no exterior source and advection drains it.
         ! Phase 14 (iface_t_model = zhang, LTNE): the porosity-weighted split
         ! of Zhang 2011 Eq.27-29 gives each phase its own interface value --
         ! the fluid phase bcs%iface_T, the solid phase bcs%iface_Ts -- so the
         ! individual fluxes are eps*F and (1-eps)*F and their sum is exactly
         ! the clear-fluid flux (energy conserving).
         is_neumann = .false.
         q_face     = 0.0_dp
         T_face     = bcs%iface_T(i)
         if ( g_iface_T_ltne ) then
            if ( present(phase) ) then
               if ( phase == 2 ) T_face = bcs%iface_Ts(i)
            end if
         end if
      end select

   end subroutine bc_face_T

   !----------------------------------------------------------------------------
   ! Set face velocity / pressure on coupling-interface faces (phase 6).
   ! The coupling driver calls these after the exchange layer has interpolated
   ! the peer solver's state onto each interface face.  Only faces belonging
   ! to a BC_INTERFACE group use these values; other faces ignore them.
   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   ! set_inlet_ramp_factor -- update the module-level mass-flow-inlet ramp
   ! factor (called once per outer SIMPLE iteration by the solver drivers).
   !----------------------------------------------------------------------------
   subroutine set_inlet_ramp_factor( f )
      real(dp), intent(in) :: f
      g_inlet_ramp_factor = max( 0.0_dp, min( 1.0_dp, f ) )
   end subroutine set_inlet_ramp_factor

   !----------------------------------------------------------------------------
   ! set_pval_ramp_factor -- update the module-level pressure-Dirichlet ramp
   ! factor (POUTLET/FARFIELD).  Ramped in lock-step with the inlet for a
   ! stable cold start with a non-zero outlet pressure; = 1 when converged.
   !----------------------------------------------------------------------------
   subroutine set_pval_ramp_factor( f )
      real(dp), intent(in) :: f
      g_pval_ramp_factor = max( 0.0_dp, min( 1.0_dp, f ) )
   end subroutine set_pval_ramp_factor

   subroutine set_interface_vel( bcs, faces, u )
      type(bc_t), intent(inout) :: bcs
      integer,    intent(in)    :: faces(:)
      real(dp),   intent(in)    :: u(3, size(faces))
      integer :: k
      do k = 1, size(faces)
         bcs%iface_vel(:, faces(k)) = u(:,k)
      end do
   end subroutine set_interface_vel

   subroutine set_interface_p( bcs, faces, p )
      type(bc_t), intent(inout) :: bcs
      integer,    intent(in)    :: faces(:)
      real(dp),   intent(in)    :: p(size(faces))
      integer :: k
      do k = 1, size(faces)
         bcs%iface_p(faces(k)) = p(k)
      end do
   end subroutine set_interface_p

   subroutine set_interface_T( bcs, faces, T )
      type(bc_t), intent(inout) :: bcs
      integer,    intent(in)    :: faces(:)
      real(dp),   intent(in)    :: T(size(faces))
      integer :: k
      do k = 1, size(faces)
         bcs%iface_T(faces(k)) = T(k)
      end do
   end subroutine set_interface_T

   !----------------------------------------------------------------------------
   ! set_interface_T_ltne -- per-face porous-phase interface temperatures from
   ! the Zhang 2011 Eq.27-29 closure: Tf (fluid phase, Eq.28 share eps*F) and
   ! Ts (solid phase, Eq.29 share (1-eps)*F).
   !----------------------------------------------------------------------------
   subroutine set_interface_T_ltne( bcs, faces, Tf, Ts )
      type(bc_t), intent(inout) :: bcs
      integer,    intent(in)    :: faces(:)
      real(dp),   intent(in)    :: Tf(size(faces)), Ts(size(faces))
      integer :: k
      do k = 1, size(faces)
         bcs%iface_T(faces(k))  = Tf(k)
         bcs%iface_Ts(faces(k)) = Ts(k)
      end do
   end subroutine set_interface_T_ltne

   !----------------------------------------------------------------------------
   ! Interface closure switches (phase 14).  Set once by the coupling driver
   ! after mix.control has been read; the standalone uns solver never touches
   ! them, so its behaviour is bit-identical to the pre-phase-14 code.
   !----------------------------------------------------------------------------
   subroutine set_iface_p_dirichlet( flag )
      logical, intent(in) :: flag
      g_iface_p_dir = flag
   end subroutine set_iface_p_dirichlet

   logical function iface_p_is_dirichlet()
      iface_p_is_dirichlet = g_iface_p_dir
   end function iface_p_is_dirichlet

   subroutine set_iface_slip_mode( flag )
      logical, intent(in) :: flag
      g_iface_slip = flag
   end subroutine set_iface_slip_mode

   logical function iface_slip_is_on()
      iface_slip_is_on = g_iface_slip
   end function iface_slip_is_on

   subroutine set_iface_t_ltne_mode( flag )
      logical, intent(in) :: flag
      g_iface_T_ltne = flag
   end subroutine set_iface_t_ltne_mode

   logical function iface_t_ltne_is_on()
      iface_t_ltne_is_on = g_iface_T_ltne
   end function iface_t_ltne_is_on

   subroutine set_iface_zhang_vel( flag )
      logical, intent(in) :: flag
      g_iface_zhang_vel = flag
   end subroutine set_iface_zhang_vel

   logical function iface_zhang_vel_is_on()
      iface_zhang_vel_is_on = g_iface_zhang_vel
   end function iface_zhang_vel_is_on

   !----------------------------------------------------------------------------
   ! set_interface_slip -- per-face tangential stress-jump conductance and the
   ! peer (free-fluid) face velocity used as the tangential target.
   !   C      : conductance [kg/s], C = mu*A/(d_f + sqrt(K)/alpha + eps*d_p)
   !            (mod_iface_law:iface_slip_conductance)
   !   u_peer : peer face velocity [m/s] (only its tangential part is used)
   !----------------------------------------------------------------------------
   subroutine set_interface_slip( bcs, faces, C, u_peer )
      type(bc_t), intent(inout) :: bcs
      integer,    intent(in)    :: faces(:)
      real(dp),   intent(in)    :: C(size(faces))
      real(dp),   intent(in)    :: u_peer(3, size(faces))
      integer :: k
      do k = 1, size(faces)
         bcs%iface_slip_C(faces(k))     = C(k)
         bcs%iface_peer_vel(:,faces(k)) = u_peer(:,k)
      end do
   end subroutine set_interface_slip

   !----------------------------------------------------------------------------
   ! bc_face_p_dirichlet -- .true. when face i is a pressure-DIRICHLET boundary
   ! face, i.e. POUTLET / FARFIELD (always) or a coupling-interface face while
   ! the phase-14 interface pressure model is active.  Used by ppe_assembly /
   ! correct_fields so that the interface Dirichlet pressure gets the same
   ! matrix, pinning and flux-correction treatment as the ordinary pressure
   ! boundaries.
   !----------------------------------------------------------------------------
   logical function bc_face_p_dirichlet( bcs, i )
      type(bc_t), intent(in) :: bcs
      integer,    intent(in) :: i
      integer :: g

      bc_face_p_dirichlet = .false.
      g = bcs%fgrp(i)
      if ( g <= 0 ) return
      select case ( bcs%gb(g)%btype )
      case ( BC_POUTLET, BC_FARFIELD )
         bc_face_p_dirichlet = .true.
      case ( BC_INTERFACE )
         bc_face_p_dirichlet = g_iface_p_dir
      end select
   end function bc_face_p_dirichlet

   !----------------------------------------------------------------------------
   ! True when the given face zone id is a coupling "interface" zone, i.e.
   ! its condition or user name contains the substring "interface"
   ! (case-insensitive).  Must stay consistent with the matching rule used in
   ! mod_uns_geometry:register_interface_zones.
   !----------------------------------------------------------------------------
   logical function zone_is_interface( m, zid )
      type(mesh_t), intent(in) :: m
      integer,      intent(in) :: zid
      integer :: i
      character(len=64) :: nm

      zone_is_interface = .false.
      do i = 1, m%nzone
         if ( m%zone(i)%id /= zid ) cycle
         nm = trim(adjustl(lower(m%zone(i)%cond_name))) // ' ' // &
              trim(adjustl(lower(m%zone(i)%user_name)))
         if ( index( nm, 'interface' ) > 0 ) zone_is_interface = .true.
         return
      end do
   end function zone_is_interface

   pure function lower( s ) result( r )
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, c
      r = s
      do i = 1, len(s)
         c = iachar( s(i:i) )
         if ( c >= iachar('A') .and. c <= iachar('Z') ) &
            r(i:i) = achar( c + 32 )
      end do
   end function lower

end module mod_uns_bc
