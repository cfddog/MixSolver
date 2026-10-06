!===============================================================================
! mod_fields.f90 -- Flow-field storage and Green-Gauss gradients (phase 3)
!
! Collocated (cell-centred) arrangement:
!   u     : cell velocity (3,ncells)
!   p     : cell pressure (ncells)
!   gp    : pressure gradient (3,ncells)
!   gu    : velocity gradient (component, direction, ncells)
!   flux  : face mass flux rho*u.S, positive from c0 to c1
!   lf    : linear interpolation weight of the c0 value (interior faces);
!           lf = d1/(d0+d1), so the closer cell gets the larger weight
!   apc   : relaxed momentum-matrix diagonal (stored for Rhie-Chow)
!   u_old    : velocity at the previous time level (transient PISO, phase 8)
!   u_old_old: velocity two time levels back (BDF2, phase 9)
!   ts_order : time scheme flag, 0=steady, 1=implicit Euler, 2=BDF2 (phase 9)
!===============================================================================
module mod_uns_fields
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_geometry
   use mod_uns_bc
   use mod_uns_control, only: ctrl_t, cell_zone_t, CZ_AUTO, CZ_POROUS
   implicit none
   private
   public :: fields_t, init_fields, compute_gradients, grad_scalar, &
             grad_vector, setup_porous_fields

   type :: fields_t
      real(dp), allocatable :: u(:,:)      ! (3,ncells) velocity
      real(dp), allocatable :: p(:)        ! (ncells)   pressure
      real(dp), allocatable :: gp(:,:)     ! (3,ncells) pressure gradient
      real(dp), allocatable :: gu(:,:,:)   ! (3,3,ncells) velocity gradient
      real(dp), allocatable :: flux(:)     ! (nfaces) mass flux, c0 -> c1 > 0
      real(dp), allocatable :: lf(:)       ! (nfaces) interp weight of c0 value
      real(dp), allocatable :: apc(:)      ! (ncells) relaxed momentum diag
      real(dp), allocatable :: u_old(:,:)     ! (3,ncells) previous time level
      real(dp), allocatable :: u_old_old(:,:)! (3,ncells) two levels back (BDF2)
      integer  :: ts_order = 0               ! 0=steady, 1=Euler, 2=BDF2
      ! --- temperature field (phase 10, passive scalar / Boussinesq) ---
      real(dp), allocatable :: T(:)          ! (ncells) temperature
      real(dp), allocatable :: T_old(:)      ! (ncells) previous time level
      real(dp), allocatable :: T_old_old(:)  ! (ncells) two levels back (BDF2)
      real(dp), allocatable :: gt(:,:)       ! (3,ncells) temperature gradient
      ! --- solid-phase temperature (LTNE, phase 12) ---
      real(dp), allocatable :: T_s(:)        ! (ncells) solid temperature
      real(dp), allocatable :: T_s_old(:)    ! (ncells) previous time level
      real(dp), allocatable :: T_s_old_old(:)! (ncells) two levels back (BDF2)
      real(dp), allocatable :: gts(:,:)      ! (3,ncells) solid temperature gradient
      ! --- per-cell porous-medium coefficients (phase 12) ---
      ! For fluid cells: porosity=1, perm=0, inertial=0, solid props=0.
      ! Populated by setup_porous_fields from ctrl%cz (matched via m%czone).
      real(dp), allocatable :: perm(:)       ! (ncells) isotropic permeability (m^2)
      real(dp), allocatable :: perm_dir(:,:) ! (3,ncells) diagonal tensor K_xx..K_zz
                                             !  (axis-aligned anisotropy; unset
                                             !   components fall back to perm)
      real(dp), allocatable :: inertial(:)   ! (ncells) Forchheimer coeff (1/m)
      real(dp), allocatable :: porosity(:)   ! (ncells) void fraction
      real(dp), allocatable :: disp_l(:)     ! (ncells) longitudinal thermal
                                             !  dispersivity (m), Bear model
      real(dp), allocatable :: disp_t(:)     ! (ncells) transverse thermal
                                             !  dispersivity (m), Bear model
      real(dp), allocatable :: bj_alpha(:)   ! (ncells) Beavers-Joseph slip
                                             !  coefficient (porous cells only)
      real(dp), allocatable :: k_s(:)        ! (ncells) solid conductivity
      real(dp), allocatable :: cp_s(:)       ! (ncells) solid specific heat
      real(dp), allocatable :: rho_s(:)      ! (ncells) solid density
      real(dp), allocatable :: h_sf(:)       ! (ncells) fluid-solid HTC
      real(dp), allocatable :: a_sf(:)       ! (ncells) specific surface area
   end type fields_t

contains

   !----------------------------------------------------------------------------
   ! Allocate and initialize all fields to a uniform state at rest
   !----------------------------------------------------------------------------
   subroutine init_fields( m, g, fld )
      type(mesh_t),  intent(in)  :: m
      type(geom_t),  intent(in)  :: g
      type(fields_t), intent(out) :: fld

      integer :: i, c0, c1, icd
      real(dp) :: d0, d1, xm

      allocate( fld%u(3,m%ncells), fld%p(m%ncells) )
      allocate( fld%gp(3,m%ncells), fld%gu(3,3,m%ncells) )
      allocate( fld%flux(m%nfaces), fld%lf(m%nfaces), fld%apc(m%ncells) )
      allocate( fld%u_old(3,m%ncells) )
      allocate( fld%u_old_old(3,m%ncells) )
      allocate( fld%T(m%ncells), fld%T_old(m%ncells), fld%T_old_old(m%ncells) )
      allocate( fld%gt(3,m%ncells) )
      allocate( fld%T_s(m%ncells), fld%T_s_old(m%ncells), fld%T_s_old_old(m%ncells) )
      allocate( fld%gts(3,m%ncells) )
      allocate( fld%perm(m%ncells), fld%inertial(m%ncells), fld%porosity(m%ncells) )
      allocate( fld%perm_dir(3,m%ncells) )
      allocate( fld%disp_l(m%ncells), fld%disp_t(m%ncells) )
      allocate( fld%bj_alpha(m%ncells) )
      allocate( fld%k_s(m%ncells), fld%cp_s(m%ncells), fld%rho_s(m%ncells) )
      allocate( fld%h_sf(m%ncells), fld%a_sf(m%ncells) )

      fld%u        = 0.0_dp
      fld%u(1,:)   = 0.3472_dp  ! TEMP DEBUG: small 1D init eigenmode test
      fld%T        = 300.0_dp   ! TEMP DEBUG
      fld%p        = 0.0_dp
      ! TEMP DEBUG: analytic pressure: 302.1 Pa in fluid region, linear in bed
      do icd = 1, m%ncells
         xm = g%xc(1,icd) * 1.0e3_dp
         if ( xm >= 100.0_dp ) fld%p(icd) = 0.03021_dp * (200.0_dp - xm)
         if ( xm <  100.0_dp ) fld%p(icd) = 3.021_dp
      end do
      fld%gp       = 0.0_dp
      fld%gu       = 0.0_dp
      fld%flux     = 0.0_dp
      fld%apc      = 1.0_dp
      fld%u_old    = 0.0_dp
      fld%u_old_old= 0.0_dp
      fld%ts_order = 0
      fld%T        = 0.0_dp
      fld%T_old    = 0.0_dp
      fld%T_old_old= 0.0_dp
      fld%gt       = 0.0_dp
      fld%T_s      = 0.0_dp
      fld%T_s_old  = 0.0_dp
      fld%T_s_old_old = 0.0_dp
      fld%gts      = 0.0_dp
      ! default: every cell is pure fluid (no porous resistance, porosity=1)
      fld%perm     = 0.0_dp
      fld%perm_dir = 0.0_dp
      fld%disp_l   = 0.0_dp
      fld%disp_t   = 0.0_dp
      fld%bj_alpha = 0.0_dp
      fld%inertial = 0.0_dp
      fld%porosity = 1.0_dp
      fld%k_s      = 0.0_dp
      fld%cp_s     = 0.0_dp
      fld%rho_s    = 0.0_dp
      fld%h_sf     = 0.0_dp
      fld%a_sf     = 0.0_dp

      do i = 1, m%nfaces
         c1 = m%f(i)%c1
         if ( c1 == 0 ) then
            fld%lf(i) = 1.0_dp
         else
            c0 = m%f(i)%c0
            d0 = norm2( g%xf(:,i) - g%xc(:,c0) )
            d1 = norm2( g%xc(:,c1) - g%xf(:,i) )
            fld%lf(i) = d1 / ( d0 + d1 )
         end if
      end do

   end subroutine init_fields

   !----------------------------------------------------------------------------
   ! Compute pressure and velocity gradients (Green-Gauss, cell based)
   !----------------------------------------------------------------------------
   subroutine compute_gradients( m, g, bcs, fld )
      type(mesh_t),   intent(in)    :: m
      type(geom_t),   intent(in)    :: g
      type(bc_t),     intent(in)    :: bcs
      type(fields_t), intent(inout) :: fld

      integer  :: i, k, comp
      real(dp), allocatable :: phif(:)

      allocate( phif(m%nfaces) )

      ! ---- pressure ------------------------------------------------------------
      do i = 1, m%nfaces
         if ( m%f(i)%c1 > 0 ) then
            phif(i) = fld%lf(i)          * fld%p(m%f(i)%c0) &
                    + (1.0_dp-fld%lf(i)) * fld%p(m%f(i)%c1)
         else
            call bc_face_p( bcs, i, fld%p(m%f(i)%c0), phif(i) )
         end if
      end do
      call grad_scalar( m, g, phif, fld%gp )

      ! ---- velocity components ---------------------------------------------------
      do comp = 1, 3
         do i = 1, m%nfaces
            if ( m%f(i)%c1 > 0 ) then
               phif(i) = fld%lf(i)          * fld%u(comp,m%f(i)%c0) &
                       + (1.0_dp-fld%lf(i)) * fld%u(comp,m%f(i)%c1)
            else
               block
                  real(dp) :: uf(3)
                  call bc_face_vel( bcs, i, g%xf(:,i), g%sf(:,i), &
                                    fld%u(:,m%f(i)%c0), uf )
                  phif(i) = uf(comp)
               end block
            end if
         end do
         call grad_scalar( m, g, phif, fld%gu(comp,:,:) )
      end do

      ! ---- temperature (Green-Gauss; Neumann walls mirror cell value) ----------
      do i = 1, m%nfaces
         if ( m%f(i)%c1 > 0 ) then
            phif(i) = fld%lf(i)          * fld%T(m%f(i)%c0) &
                    + (1.0_dp-fld%lf(i)) * fld%T(m%f(i)%c1)
         else
            block
               real(dp) :: tf, qf
               logical  :: neum
               call bc_face_T( bcs, i, g%xf(:,i), fld%T(m%f(i)%c0), tf, qf, neum )
               phif(i) = tf
            end block
         end if
      end do
      call grad_scalar( m, g, phif, fld%gt )

      ! ---- solid temperature (LTNE; same boundary treatment as fluid T) -------
      do i = 1, m%nfaces
         if ( m%f(i)%c1 > 0 ) then
            phif(i) = fld%lf(i)          * fld%T_s(m%f(i)%c0) &
                    + (1.0_dp-fld%lf(i)) * fld%T_s(m%f(i)%c1)
         else
            block
               real(dp) :: tf, qf
               logical  :: neum
               ! solid phase shares the wall thermal BC (T or q) with the fluid
               call bc_face_T( bcs, i, g%xf(:,i), fld%T_s(m%f(i)%c0), tf, qf, neum )
               phif(i) = tf
            end block
         end if
      end do
      call grad_scalar( m, g, phif, fld%gts )

      deallocate( phif )

   end subroutine compute_gradients

   !----------------------------------------------------------------------------
   ! Green-Gauss gradient of a cell-centred scalar given its face values:
   ! grad_P = (1/V) sum_f phi_f S_f (S outward from the owning cell)
   !----------------------------------------------------------------------------
   subroutine grad_scalar( m, g, phif, grad )
      type(mesh_t), intent(in)  :: m
      type(geom_t), intent(in)  :: g
      real(dp),     intent(in)  :: phif(:)
      real(dp),     intent(out) :: grad(:,:)    ! (3,ncells)

      integer :: i, c0, c1

      grad = 0.0_dp
      do i = 1, m%nfaces
         c0 = m%f(i)%c0
         grad(:,c0) = grad(:,c0) + phif(i) * g%sf(:,i)
         c1 = m%f(i)%c1
         if ( c1 > 0 ) grad(:,c1) = grad(:,c1) - phif(i) * g%sf(:,i)
      end do
      do i = 1, m%ncells
         grad(:,i) = grad(:,i) / g%vol(i)
      end do

   end subroutine grad_scalar

   !----------------------------------------------------------------------------
   ! Green-Gauss gradient of a cell-centred vector field given face values
   !----------------------------------------------------------------------------
   subroutine grad_vector( m, g, vphif, grad )
      type(mesh_t), intent(in)  :: m
      type(geom_t), intent(in)  :: g
      real(dp),     intent(in)  :: vphif(:,:)   ! (3,nfaces)
      real(dp),     intent(out) :: grad(:,:,:)  ! (3,3,ncells) (comp,dir)

      integer :: i, c0, c1, comp

      grad = 0.0_dp
      do i = 1, m%nfaces
         c0 = m%f(i)%c0
         c1 = m%f(i)%c1
         do comp = 1, 3
            grad(comp,:,c0) = grad(comp,:,c0) + vphif(comp,i) * g%sf(:,i)
            if ( c1 > 0 ) &
               grad(comp,:,c1) = grad(comp,:,c1) - vphif(comp,i) * g%sf(:,i)
         end do
      end do
      do i = 1, m%ncells
         grad(:,:,i) = grad(:,:,i) / g%vol(i)
      end do

   end subroutine grad_vector

   !----------------------------------------------------------------------------
   ! Populate per-cell porous-medium coefficients from ctrl%cz.
   !
   ! For every cell, the owning cell-zone id (m%czone) is matched against the
   ! control-file cell_zone entries (ctrl%cz) by zone id, then by zone name
   ! (cond_name / user_name, case-insensitive).  Coefficients are copied to a
   ! cell when (a) the cell's resolved block type is CZ_POROUS (explicit
   ! cell_zone type or a "VC: porous" zone-name tag) and (b) the matched entry
   ! carries porous coefficients (ztype = CZ_POROUS or CZ_AUTO).  Otherwise the
   ! arrays keep the fluid default (porosity=1, perm=0, inertial=0,
   ! solid props=0) set in init_fields.
   !
   ! This routine must be called AFTER resolve_cell_zones (which sets m%cztype)
   ! and AFTER init_fields.
   !----------------------------------------------------------------------------
   subroutine setup_porous_fields( m, ctrl, fld )
      type(mesh_t),    intent(in)    :: m
      type(ctrl_t),    intent(in)    :: ctrl
      type(fields_t),  intent(inout) :: fld

      integer :: i, iz, j, czid
      logical :: found
      character(len=32) :: want

      if ( .not. allocated(m%czone) ) return

      do i = 1, m%ncells
         czid = m%czone(i)
         ! find the mesh cell-zone table entry with this id
         iz = 0
         do j = 1, m%nczone
            if ( m%czt(j)%id == czid ) then
               iz = j
               exit
            end if
         end do
         if ( iz == 0 ) cycle   ! unknown zone -> keep fluid default

         ! find the matching control-file cell_zone entry
         found = .false.
         do j = 1, ctrl%ncz
            if ( ctrl%cz(j)%id /= 0 ) then
               found = ( ctrl%cz(j)%id == m%czt(iz)%id )
            else
               want = trim(adjustl(ctrl%cz(j)%name))
               found = ( trim(adjustl(m%czt(iz)%cond_name)) == want .or. &
                         trim(adjustl(m%czt(iz)%user_name)) == want )
            end if
            if ( found ) exit
         end do

         if ( found .and. m%cztype(i) == CZ_POROUS .and. &
              ( ctrl%cz(j)%ztype == CZ_POROUS .or. ctrl%cz(j)%ztype == CZ_AUTO ) ) then
            fld%perm(i)     = ctrl%cz(j)%perm
            ! diagonal permeability tensor: an explicitly given component wins,
            ! unset components (0) fall back to the scalar perm
            fld%perm_dir(1,i) = merge( ctrl%cz(j)%perm_xx, ctrl%cz(j)%perm, &
                                       ctrl%cz(j)%perm_xx > 0.0_dp )
            fld%perm_dir(2,i) = merge( ctrl%cz(j)%perm_yy, ctrl%cz(j)%perm, &
                                       ctrl%cz(j)%perm_yy > 0.0_dp )
            fld%perm_dir(3,i) = merge( ctrl%cz(j)%perm_zz, ctrl%cz(j)%perm, &
                                       ctrl%cz(j)%perm_zz > 0.0_dp )
            fld%inertial(i) = ctrl%cz(j)%inertial
            fld%porosity(i) = ctrl%cz(j)%porosity
            fld%disp_l(i)   = ctrl%cz(j)%disp_l
            fld%disp_t(i)   = ctrl%cz(j)%disp_t
            fld%bj_alpha(i) = ctrl%cz(j)%bj_alpha
            fld%k_s(i)      = ctrl%cz(j)%k_s
            fld%cp_s(i)     = ctrl%cz(j)%cp_s
            fld%rho_s(i)    = ctrl%cz(j)%rho_s
            fld%h_sf(i)     = ctrl%cz(j)%h_sf
            fld%a_sf(i)     = ctrl%cz(j)%a_sf
         end if
         ! else: fluid defaults already set by init_fields
      end do

   end subroutine setup_porous_fields

end module mod_uns_fields
