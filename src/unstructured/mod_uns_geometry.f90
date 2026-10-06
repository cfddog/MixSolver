!===============================================================================
! mod_geometry.f90 -- Geometric quantities for face-based unstructured FVM
!
! Computes:
!   xf, sf, area : face centroids, outward (c0 -> c1) area vectors, areas
!   vol, xc      : cell volumes and centroids (pyramid decomposition)
!   rc0, rc1     : vectors from cell centroids to face centroids
!
! geom_stats also returns the geometric conservation law (GCL) residual,
! max over cells of |sum_f (+/-)sf| / sum_f |sf|.
!===============================================================================
module mod_uns_geometry
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_interface, only: Interface_FACE_TYPE, Interface_List, Num_Interface, &
                            PEER_UNS
   implicit none
   private
   public :: geom_t, compute_geometry, geom_stats, register_interface_zones

   type :: geom_t
      real(dp), allocatable :: xf(:,:)    ! face centroids, xf(1:3,1:nfaces)
      real(dp), allocatable :: sf(:,:)    ! face area vectors (c0 -> c1)
      real(dp), allocatable :: area(:)    ! face areas |sf|
      real(dp), allocatable :: vol(:)     ! cell volumes
      real(dp), allocatable :: xc(:,:)    ! cell centroids
      real(dp), allocatable :: rc0(:,:)   ! xf - xc0
      real(dp), allocatable :: rc1(:,:)   ! xf - xc1 (unset for boundary)
   end type geom_t

   real(dp), parameter :: THIRD  = 1.0_dp/3.0_dp
   real(dp), parameter :: QRTR   = 1.0_dp/4.0_dp

contains

   !----------------------------------------------------------------------------
   ! Compute all geometric quantities
   !----------------------------------------------------------------------------
   subroutine compute_geometry( m, g, ier )
      type(mesh_t), intent(in)  :: m
      type(geom_t), intent(out) :: g
      integer,      intent(out) :: ier

      integer  :: i, k
      real(dp) :: cen(3), s(3), a, dotv
      real(dp), allocatable :: pcen(:,:)   ! provisional interior points

      ier = 0
      allocate( g%xf(3,m%nfaces), g%sf(3,m%nfaces), g%area(m%nfaces) )
      allocate( g%vol(m%ncells), g%xc(3,m%ncells) )
      allocate( g%rc0(3,m%nfaces), g%rc1(3,m%nfaces) )

      ! ---- face centroids, area vectors, areas -------------------------------
      do i = 1, m%nfaces
         call face_props( m, i, cen, s, a )
         g%xf(:,i) = cen
         g%sf(:,i) = s
         g%area(i) = a
      end do

      ! ---- provisional interior point per cell (node average) ----------------
      allocate( pcen(3,m%ncells) )
      call prov_centers( m, pcen )

      ! ---- orient area vectors: sf must point from c0 to c1 ------------------
      do i = 1, m%nfaces
         dotv = dot_product( g%sf(:,i), g%xf(:,i) - pcen(:,m%f(i)%c0) )
         if ( dotv < 0.0_dp ) g%sf(:,i) = -g%sf(:,i)
      end do

      ! ---- cell volumes and centroids (pyramid decomposition) ----------------
      g%vol = 0.0_dp
      g%xc  = 0.0_dp
      do i = 1, m%nfaces
         call accumulate_cell( g, pcen(:,m%f(i)%c0), +1, i, m%f(i)%c0 )
         if ( m%f(i)%c1 > 0 ) &
            call accumulate_cell( g, pcen(:,m%f(i)%c1), -1, i, m%f(i)%c1 )
      end do
      do k = 1, m%ncells
         if ( g%vol(k) > 0.0_dp ) g%xc(:,k) = g%xc(:,k) / g%vol(k)
      end do

      ! ---- centroid-to-face vectors ------------------------------------------
      do i = 1, m%nfaces
         g%rc0(:,i) = g%xf(:,i) - g%xc(:,m%f(i)%c0)
         if ( m%f(i)%c1 > 0 ) &
            g%rc1(:,i) = g%xf(:,i) - g%xc(:,m%f(i)%c1)
      end do

      deallocate( pcen )
   end subroutine compute_geometry

   !----------------------------------------------------------------------------
   ! Face centroid, area vector and area from its polygon (planar assumed)
   !----------------------------------------------------------------------------
   subroutine face_props( m, i, cen, s, a )
      type(mesh_t), intent(in)  :: m
      integer,      intent(in)  :: i
      real(dp),     intent(out) :: cen(3), s(3), a

      integer  :: n, j, jn, i1
      real(dp) :: p1(3), p2(3), p3(3), st(3), at

      n = m%f(i)%nn

      ! Newell area vector: s = 1/2 * sum_j cross(p_j, p_{j+1})
      s = 0.0_dp
      do j = 1, n
         jn = mod(j, n) + 1
         p1 = m%x(:, m%f(i)%nodes(j))
         p2 = m%x(:, m%f(i)%nodes(jn))
         s = s + 0.5_dp * cross( p1, p2 )
      end do
      a = norm2( s )

      ! area-weighted centroid by fan triangulation from node 1
      cen = 0.0_dp
      at  = 0.0_dp
      i1  = m%f(i)%nodes(1)
      p1  = m%x(:, i1)
      do j = 2, n-1
         p2 = m%x(:, m%f(i)%nodes(j))
         p3 = m%x(:, m%f(i)%nodes(j+1))
         st = cross( p2 - p1, p3 - p1 ) * 0.5_dp
         cen = cen + norm2( st ) * ( p1 + p2 + p3 ) / 3.0_dp
         at  = at + norm2( st )
      end do
      if ( at > 0.0_dp ) cen = cen / at

   end subroutine face_props

   ! provisional interior point per cell: average of its nodes
   subroutine prov_centers( m, pcen )
      type(mesh_t), intent(in)  :: m
      real(dp),     intent(out) :: pcen(:,:)
      integer :: i, j, c0, c1
      integer, allocatable :: cnt(:)

      allocate( cnt(m%ncells) )
      pcen = 0.0_dp
      cnt  = 0
      do i = 1, m%nfaces
         c0 = m%f(i)%c0; c1 = m%f(i)%c1
         do j = 1, m%f(i)%nn
            pcen(:,c0) = pcen(:,c0) + m%x(:,m%f(i)%nodes(j))
            cnt(c0) = cnt(c0) + 1
            if ( c1 > 0 ) then
               pcen(:,c1) = pcen(:,c1) + m%x(:,m%f(i)%nodes(j))
               cnt(c1) = cnt(c1) + 1
            end if
         end do
      end do
      do i = 1, m%ncells
         if ( cnt(i) > 0 ) pcen(:,i) = pcen(:,i) / real(cnt(i), dp)
      end do
      deallocate( cnt )
   end subroutine prov_centers

   ! accumulate one pyramid (apex p, base = face i) into cell volume/centroid
   subroutine accumulate_cell( g, p, sgn, i, cell )
      type(geom_t), intent(inout) :: g
      real(dp),     intent(in)    :: p(3)
      integer,      intent(in)    :: sgn, i, cell

      real(dp) :: dv

      ! signed pyramid volume; outward sf (relative to the cell) gives dv > 0
      dv = THIRD * dot_product( g%sf(:,i), g%xf(:,i) - p ) * real(sgn, dp)
      g%vol(cell) = g%vol(cell) + dv
      ! pyramid centroid = (3*face_centroid + apex)/4, volume-weighted
      g%xc(:,cell) = g%xc(:,cell) + dv * ( ( 1.0_dp - QRTR )*g%xf(:,i) &
                                           + QRTR*p )
   end subroutine accumulate_cell

   !----------------------------------------------------------------------------
   ! Statistics + geometric conservation law check
   !----------------------------------------------------------------------------
   subroutine geom_stats( m, g, vol_sum, vol_min, vol_max, gcl_max, ier )
      type(mesh_t), intent(in)  :: m
      type(geom_t), intent(in)  :: g
      real(dp),     intent(out) :: vol_sum, vol_min, vol_max, gcl_max
      integer,      intent(out) :: ier

      integer  :: k, i, nbad
      real(dp) :: ssum(3), nrm, ratio

      vol_sum = sum( g%vol )
      vol_min = minval( g%vol )
      vol_max = maxval( g%vol )

      gcl_max = 0.0_dp
      do k = 1, m%ncells
         ssum = 0.0_dp
         nrm  = 0.0_dp
         do i = 1, m%nfaces
            if ( m%f(i)%c0 == k ) then
               ssum = ssum + g%sf(:,i)
               nrm  = nrm + g%area(i)
            else if ( m%f(i)%c1 == k ) then
               ssum = ssum - g%sf(:,i)
               nrm  = nrm + g%area(i)
            end if
         end do
         if ( nrm > 0.0_dp ) then
            ratio = norm2( ssum ) / nrm
            gcl_max = max( gcl_max, ratio )
         end if
      end do

      ier = 0
      nbad = count( g%vol <= 0.0_dp )
      if ( nbad > 0 ) then
         write(*,'(a,i0)') 'ERROR: non-positive cell volumes: ', nbad
         ier = 30
      end if
   end subroutine geom_stats

   ! vector cross product
   pure function cross( a, b ) result( c )
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: c(3)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross

   !----------------------------------------------------------------------------
   ! Register Fluent CAS "interface" zones into the common coupling list.
   !
   ! A zone is treated as a coupling interface when either its user_name or
   ! its cond_name (case-insensitive) contains the substring "interface".
   ! For every boundary face belonging to such a zone, one Interface_FACE_TYPE
   ! entry is appended to Interface_List with:
   !   solver    = PEER_UNS
   !   block_no  = the Fluent zone id
   !   centroid  = g%xf(:,face)  (area-weighted face centroid)
   !   normal    = unit normal of g%sf(:,face)
   !   area      = g%area(face)
   !   bbox_min/max = bounding box of the face nodes
   !
   ! This is a pure registration call: it does not modify the mesh or geometry,
   ! and it does not apply any boundary condition.  The coupling layer
   ! (phase 4) will later match these entries against the structured side's
   ! gridgen generic:8 faces.
   !----------------------------------------------------------------------------
   subroutine register_interface_zones( m, g )
      type(mesh_t), intent(in) :: m
      type(geom_t), intent(in) :: g

      integer :: zi, fi, j, niface, pos, z_nf
      type(Interface_FACE_TYPE), allocatable :: tmp(:)
      real(dp) :: nrm, bmin(3), bmax(3)
      real(dp) :: z_area, z_c(3), z_n(3), z_bmin(3), z_bmax(3)
      character(len=64) :: nm

      ! first pass: count interface faces
      niface = 0
      do zi = 1, m%nzone
         nm = trim(adjustl( lowercase(m%zone(zi)%user_name) )) // &
              ' ' // trim(adjustl( lowercase(m%zone(zi)%cond_name) ))
         if ( index( nm, 'interface' ) > 0 ) then
            do fi = 1, m%nfaces
               if ( m%f(fi)%zone == m%zone(zi)%id ) niface = niface + 1
            end do
         end if
      end do

      if ( niface == 0 ) then
         write(*,'(a)') '  No interface zones found in CAS file.'
         return
      end if

      ! append to the (possibly already populated) Interface_List
      pos = Num_Interface
      if ( .not. allocated(Interface_List) ) then
         allocate( Interface_List(niface) )
      else
         allocate( tmp(pos + niface) )
         tmp(1:pos) = Interface_List
         call move_alloc( tmp, Interface_List )
      end if

      ! second pass: fill entries
      niface = 0
      do zi = 1, m%nzone
         nm = trim(adjustl( lowercase(m%zone(zi)%user_name) )) // &
              ' ' // trim(adjustl( lowercase(m%zone(zi)%cond_name) ))
         if ( index( nm, 'interface' ) > 0 ) then
            ! per-zone aggregate diagnostics (area-weighted)
            z_area = 0.0_dp
            z_c    = 0.0_dp
            z_n    = 0.0_dp
            z_bmin =  huge(1.0_dp)
            z_bmax = -huge(1.0_dp)
            z_nf   = 0
            do fi = 1, m%nfaces
               if ( m%f(fi)%zone /= m%zone(zi)%id ) cycle
               niface = niface + 1
               z_nf   = z_nf + 1
               pos = pos + 1

               ! face centroid and area from geom
               Interface_List(pos)%centroid = g%xf(:,fi)
               nrm = norm2( g%sf(:,fi) )
               if ( nrm > 0.0_dp ) then
                  Interface_List(pos)%normal = g%sf(:,fi) / nrm
               else
                  Interface_List(pos)%normal = 0.0_dp
               end if
               Interface_List(pos)%area = g%area(fi)

               ! bounding box from face nodes
               bmin =  huge(1.0_dp)
               bmax = -huge(1.0_dp)
               do j = 1, m%f(fi)%nn
                  bmin = min( bmin, m%x(:, m%f(fi)%nodes(j)) )
                  bmax = max( bmax, m%x(:, m%f(fi)%nodes(j)) )
               end do
               Interface_List(pos)%bbox_min = bmin
               Interface_List(pos)%bbox_max = bmax

               ! identity fields
               Interface_List(pos)%solver   = PEER_UNS
               Interface_List(pos)%block_no = m%zone(zi)%id
               Interface_List(pos)%face     = 0
               Interface_List(pos)%f_no     = 0
               Interface_List(pos)%ib = 0;  Interface_List(pos)%ie = 0
               Interface_List(pos)%jb = 0;  Interface_List(pos)%je = 0
               Interface_List(pos)%kb = 0;  Interface_List(pos)%ke = 0
               Interface_List(pos)%match_state = 0
               Interface_List(pos)%peer_id     = 0

               ! face vertex coordinates (phase 4: needed for projection +
               ! bilinear weights on the peer side when this face is the
               ! interpolation host)
               Interface_List(pos)%nv = m%f(fi)%nn
               allocate( Interface_List(pos)%verts(3, m%f(fi)%nn) )
               do j = 1, m%f(fi)%nn
                  Interface_List(pos)%verts(:,j) = m%x(:, m%f(fi)%nodes(j))
               end do

               ! zone aggregates
               z_area = z_area + g%area(fi)
               z_c    = z_c + g%area(fi) * g%xf(:,fi)
               z_n    = z_n + g%sf(:,fi)
               z_bmin = min( z_bmin, bmin )
               z_bmax = max( z_bmax, bmax )
            end do
            write(*,'(a,i0,a,a)') '  Registered interface zone ', &
               m%zone(zi)%id, ' : ', trim(m%zone(zi)%user_name)
            write(*,'(a,i0,a)')    '    [uns] faces   : ', z_nf
            if ( z_area > 0.0_dp ) &
               write(*,'(a,3(1x,es12.4))') '          centroid : ', z_c/z_area
            nrm = norm2(z_n)
            if ( nrm > 0.0_dp ) &
               write(*,'(a,3(1x,es12.4))') '          net normal: ', z_n/nrm
            write(*,'(a,es12.4)')        '          area     : ', z_area
            write(*,'(a,3(1x,es12.4))') '          bbox min : ', z_bmin
            write(*,'(a,3(1x,es12.4))') '          bbox max : ', z_bmax
         end if
      end do

      Num_Interface = pos
      write(*,'(a,i0)') '  Total unstructured interface faces registered: ', niface

   end subroutine register_interface_zones

   ! case-insensitive helper (local to this module)
   pure function lowercase( s ) result( r )
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, c
      r = s
      do i = 1, len(s)
         c = iachar( s(i:i) )
         if ( c >= iachar('A') .and. c <= iachar('Z') ) &
            r(i:i) = achar( c + 32 )
      end do
   end function lowercase

end module mod_uns_geometry
