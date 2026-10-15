!===============================================================================
! mod_output.f90 -- Result output and benchmark comparison (phase 4)
!
! Contents:
!   vtk_write     -- write the solution as a VTU (VTK XML unstructured grid)
!                    through mod_vtu_writer, which emits the layout the VTK
!                    XML specification prescribes: pure node-id <connectivity>
!                    (no legacy per-cell vertex-count token), cumulative
!                    <offsets> and <types>.  Every cell is decomposed into
!                    tetrahedra by fanning each (outward-oriented) cell face
!                    from its first node and connecting the triangles to the
!                    cell centroid, which is appended as an extra point.
!                    VTK_TETRA cells are universally supported -- unlike
!                    VTK_POLYHEDRON, whose legacy face-stream encoding modern
!                    ParaView (VTK 9 based) silently drops, leaving an
!                    empty-looking dataset.
!                    With `native = .true.` (stage 17, control `out_cells =
!                    native`) the same routine instead writes ONE VTK cell per
!                    mesh cell with its real type (4 = hex -> 12, 6 = wedge ->
!                    13): the canonical VTK node order is rebuilt from the face
!                    topology and verified (every face reproduced + positive
!                    volume equal to the solver's cell volume).  Unmapped Fluent
!                    codes (1/2/3/5 = 0-D/1-D/2-D entities, 0 = unknown) are a
!                    hard error; polyhedra (7) fall back to the fan above.
!   probe_velocity-- inverse-distance interpolation of the velocity field
!                    from the K nearest cell centroids.
!   ghia_compare  -- compare the u-velocity on the vertical centerline of a
!                    lid-driven cavity against the Ghia et al. (1982)
!                    benchmark tables for Re = 100 and Re = 1000.
!
! Conventions:
!   - A tet (a,b,c,apex) has positive volume in the VTK sense when
!     (b-a)x(c-a) . (apex-a) > 0. In this .cas convention the face node
!     order turns out to be INWARD for c0 (verified by signed tet
!     volumes), hence the fan (n1,nj,nj+1) is used as stored for c0 and
!     reversed for c1, with the cell centroid as apex.
!   - VTK point ids are 0-based; mesh_t node ids are 1-based. The centroid
!     of cell c is appended as point id (nnodes + c - 1), 0-based.
!   - vtu_write takes the connectivity in the MESH numbering (1-based) and
!     performs the 0-based conversion itself; vtk_write therefore passes the
!     raw face node ids and the 1-based centroid point index (nnodes + c).
!     Passing ids that were already converted subtracts one from every node
!     id twice, which shifts the whole mesh (caught on 2026-10-15 by signed
!     tet volumes: sum|vol| stopped matching the domain volume).
!   - History: until 2026-10-15 the <connectivity> block carried the legacy
!     VTK "CELLS" vertex-count token in front of every cell (4,a,b,c,apex).
!     ParaView tolerated it (it recomputes the cell location from the cell
!     type), but Tecplot 360 read that token as a node id and rendered a
!     broken mesh.  Fixed by writing the spec-conformant layout -- see
!     worklog 2026-10-15 and mod_vtu_writer.
!===============================================================================
module mod_uns_output
   use mod_precision, only: dp, ip, pi
   use mod_vtu_writer, only: vtu_write, vtu_data, VTK_TETRA
   use mod_cell_types, only: fluent_to_vtk, fluent_code_name, &
                             cell_type_supported, cell_type_is_polyhedron
   use mod_uns_mesh
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_fields
   use mod_uns_control
   use mod_uns_bc
   implicit none
   private
   public :: vtk_write, probe_velocity, ghia_compare, &
             thermal_wall_report

contains

   !----------------------------------------------------------------------------
   ! Write the solution as a VTU file (VTK XML) with every cell tetrahedralized.
   ! The XML/byte layout is produced by mod_vtu_writer; see the module header
   ! there for the format requirements (pure node-id connectivity, cumulative
   ! offsets, per-cell types).
   !   binary = .false. (default) -> inline ASCII, Float64 (out_format = vtu)
   !   binary = .true.            -> appended raw binary, Float32
   !                                 (out_format = vtu_binary; ~4x smaller)
   !   native = .true.            -> real cell types 1:1 (out_cells = native,
   !                                 stage 17); orthogonal to binary
   !----------------------------------------------------------------------------
   subroutine vtk_write( m, conn, g, fld, filename, ier, binary, native )
      type(mesh_t),    intent(in)  :: m
      type(conn_t),    intent(in)  :: conn
      type(geom_t),    intent(in)  :: g
      type(fields_t),  intent(in)  :: fld
      character(len=*), intent(in) :: filename
      integer,         intent(out) :: ier
      logical, optional, intent(in):: binary
      logical, optional, intent(in):: native

      integer  :: j, k, c, f, nn, nsub, ntets, cen_id, npts, ncd, it
      integer  :: tri(3)
      logical  :: has_solid, has_ts, bin, nat, poly
      real(dp), allocatable :: coords(:,:)
      integer,  allocatable :: cell_conn(:), cell_off(:), cell_type(:)
      type(vtu_data), allocatable :: cdata(:)

      ier = 0
      bin = .false.
      if ( present(binary) ) bin = binary

      ! ---- stage 17: optional native (real cell type) 1:1 output ------------
      ! `native = .true.` (control `out_cells = native`) writes every mesh cell
      ! once, with its real VTK type, instead of the tetrahedral fan below.
      ! A mesh containing polyhedra (Fluent code 7) has no usable VTK
      ! encoding, so such a file falls back to the fan with one warning; every
      ! other unmapped code is a hard error inside vtk_write_native.
      ! The default path (native absent or .false.) is exactly the code below.
      nat = .false.
      if ( present(native) ) nat = native
      if ( nat ) then
         poly = .false.
         if ( allocated(m%ctype) ) then
            do c = 1, m%ncells
               if ( cell_type_is_polyhedron( m%ctype(c) ) ) then
                  poly = .true.
                  exit
               end if
            end do
         end if
         if ( poly ) then
            write(*,'(a)') 'vtk_write: WARNING - mesh contains polyhedral ' // &
               'cells (Fluent code 7), which have no usable VTK encoding; ' // &
               'falling back to the tetrahedral fan for this file.'
            nat = .false.
         end if
      end if
      if ( nat ) then
         call vtk_write_native( m, conn, g, fld, filename, ier, bin )
         return
      end if

      ! total number of tetrahedra: each face contributes (nn-2) triangles
      ! NB: cf is 0-offset CSR -- entries of cell c are cf(cf_ptr(c)+1 : cf_ptr(c+1))
      ntets = 0
      do c = 1, m%ncells
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            f = conn%cf(k)
            ntets = ntets + m%f(f)%nn - 2
         end do
      end do

      ! ---- which cell-data arrays are written (conditions unchanged) ---------
      !      solid temperature: LTNE fluid-solid exchange active (h_sf*a_sf>0)
      !      OR a solid conduction region (CZ_SOLID) present
      has_solid = .false.
      do c = 1, m%ncells
         if ( is_solid_cell( m, fld, c ) ) then
            has_solid = .true.
            exit
         end if
      end do
      has_ts = ( maxval( fld%h_sf * fld%a_sf ) > 0.0_dp .or. has_solid )

      ! ---- points: mesh nodes, then one centroid point per cell --------------
      npts = m%nnodes + m%ncells
      allocate( coords(3, npts) )
      coords(:, 1:m%nnodes)      = m%x(:, 1:m%nnodes)
      coords(:, m%nnodes+1:npts) = g%xc(:, 1:m%ncells)

      ! ---- cells: tetrahedra from face fan + cell centroid -------------------
      !      cell_conn holds the nodes of every sub-tetrahedron with the SAME
      !      1-based node numbering as mesh_t (vtu_write converts to the 0-based
      !      ids required by VTK XML itself), <offsets> the cumulative node count
      !      (4*it) and <types> the VTK cell code (10 = VTK_TETRA).  The legacy
      !      VTK "CELLS" style (a leading vertex-count token per cell) is
      !      deliberately NOT used.
      allocate( cell_conn(4*ntets), cell_off(ntets), cell_type(ntets) )
      it = 0
      do c = 1, m%ncells
         cen_id = m%nnodes + c                ! 1-based id of centroid point
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            f  = conn%cf(k)
            nn = m%f(f)%nn
            do j = 2, nn-1                    ! fan from face node 1
               if ( m%f(f)%c0 == c ) then
                  ! NB: .cas face node order is INWARD for c0 (verified by
                  ! signed tet volumes), so the fan (n1,nj,nj+1) + centroid
                  ! apex gives a positively oriented tet
                  tri = [ m%f(f)%nodes(1), m%f(f)%nodes(j), m%f(f)%nodes(j+1) ]
               else if ( m%f(f)%c1 == c ) then
                  tri = [ m%f(f)%nodes(1), m%f(f)%nodes(j+1), m%f(f)%nodes(j) ]
               else
                  ! should not happen: cf lists both owners of every face
                  write(*,'(a,i0,a,i0)') 'VTK ERROR: face ', f, &
                     ' not adjacent to cell ', c
                  ier = 2
                  return
               end if
               it = it + 1
               cell_conn(4*it-3 : 4*it) = [ tri(1), tri(2), tri(3), cen_id ]
               cell_off(it)  = 4*it
               cell_type(it) = VTK_TETRA
            end do
         end do
      end do

      ! ---- cell data: original cell values replicated on their sub-tets ------
      !      pressure / velocity / temperature are always written; the solid
      !      temperature (LTNE exchange active or a CZ_SOLID region present) and
      !      the solid volumetric heat source (phase 16, solid region only) are
      !      appended conditionally, exactly as before.
      ncd = 3 + merge(1,0,has_ts) + merge(1,0,has_solid)
      allocate( cdata(ncd) )
      cdata(1)%name  = 'pressure';    cdata(1)%ncomp = 1
      cdata(2)%name  = 'velocity';    cdata(2)%ncomp = 3
      cdata(3)%name  = 'temperature'; cdata(3)%ncomp = 1
      if ( has_ts )    cdata(4)%name   = 'temperature_solid'
      if ( has_solid ) cdata(ncd)%name = 'heat_source_solid'
      do j = 1, ncd
         allocate( cdata(j)%val(ntets, cdata(j)%ncomp) )
      end do

      it = 0
      do c = 1, m%ncells
         nsub = 0
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            nsub = nsub + m%f(conn%cf(k))%nn - 2
         end do
         do j = 1, nsub
            it = it + 1
            cdata(1)%val(it,1)   = fld%p(c)
            cdata(2)%val(it,1:3) = fld%u(:,c)
            cdata(3)%val(it,1)   = fld%T(c)
            if ( has_ts )    cdata(4)%val(it,1)   = fld%T_s(c)
            if ( has_solid ) cdata(ncd)%val(it,1) = fld%q_src(c)
         end do
      end do
      ! ---- write the file through the generic XML writer ----------------------
      call vtu_write( filename, coords, cell_conn, cell_off, cell_type, &
                      cell_data = cdata, binary = bin, ier = ier )
      if ( ier /= 0 ) return

      if ( bin ) then
         write(*,'(a)') '  wrote VTU file : ' // trim(filename) // &
                        '  (appended raw binary, Float32)'
      else
         write(*,'(a)') '  wrote VTU file : ' // trim(filename)
      end if
   end subroutine vtk_write

   !----------------------------------------------------------------------------
   ! Velocity at an arbitrary point: inverse-distance weighted interpolation
   ! from the KNN nearest cell centroids.
   !----------------------------------------------------------------------------
   subroutine probe_velocity( m, g, fld, xq, uq )
      type(mesh_t),    intent(in)  :: m
      type(geom_t),    intent(in)  :: g
      type(fields_t),  intent(in)  :: fld
      real(dp),        intent(in)  :: xq(3)
      real(dp),        intent(out) :: uq(3)

      integer,  parameter :: KNN = 8
      integer  :: c, k, j
      real(dp) :: d2, d2best(KNN), w, wsum
      integer  :: idbest(KNN)

      d2best = huge(1.0_dp)
      idbest = 0
      do c = 1, m%ncells
         d2 = sum( ( g%xc(:,c) - xq )**2 )
         do k = 1, KNN                        ! insert into sorted top-K list
            if ( d2 < d2best(k) ) then
               do j = KNN, k+1, -1
                  d2best(j) = d2best(j-1)
                  idbest(j) = idbest(j-1)
               end do
               d2best(k) = d2
               idbest(k) = c
               exit
            end if
         end do
      end do

      uq = 0.0_dp
      wsum = 0.0_dp
      do k = 1, KNN
         if ( idbest(k) == 0 ) cycle
         w = 1.0_dp / ( d2best(k) + 1.0e-12_dp )
         uq = uq + w * fld%u(:,idbest(k))
         wsum = wsum + w
      end do
      if ( wsum > 0.0_dp ) uq = uq / wsum
   end subroutine probe_velocity

   !----------------------------------------------------------------------------
   ! Compare u-velocity on the vertical centerline (x=0.5, z=0.5) with the
   ! Ghia, Ghia & Shin (1982) benchmark. Assumes unit cavity size and unit
   ! lid velocity, i.e. Re = rho / mu.
   !----------------------------------------------------------------------------
   subroutine ghia_compare( m, g, fld, re, ier )
      type(mesh_t),    intent(in)  :: m
      type(geom_t),    intent(in)  :: g
      type(fields_t),  intent(in)  :: fld
      real(dp),        intent(in)  :: re
      integer,         intent(out) :: ier

      integer,  parameter :: NPTS = 17
      integer  :: i
      real(dp) :: yg(NPTS), u100(NPTS), u1000(NPTS), uref(NPTS)
      real(dp) :: xq(3), uq(3), l2num, l2den, err

      ier = 0
      ! --- Ghia, Ghia & Shin (1982), JCP 48:387-411, Tables I & II ------
      ! u-velocity on the vertical centerline (x=0.5); both Re tables
      ! share the same 17 stations of the 129x129 benchmark grid.
      yg    = [ 1.0000_dp, 0.9766_dp, 0.9688_dp, 0.9609_dp, 0.9531_dp, &
                0.8516_dp, 0.7344_dp, 0.6172_dp, 0.5000_dp, 0.4531_dp, &
                0.2813_dp, 0.1719_dp, 0.1016_dp, 0.0703_dp, 0.0625_dp, &
                0.0547_dp, 0.0000_dp ]
      u100  = [ 1.00000_dp,  0.84123_dp,  0.78871_dp,  0.68722_dp, &
                0.60751_dp,  0.23151_dp, -0.03393_dp, -0.14320_dp, &
               -0.20533_dp, -0.21090_dp, -0.15662_dp, -0.10150_dp, &
               -0.06434_dp, -0.04775_dp, -0.04192_dp, -0.03717_dp, &
                0.00000_dp ]
      u1000 = [ 1.00000_dp,  0.65928_dp,  0.57492_dp,  0.51117_dp, &
                0.46604_dp,  0.33304_dp,  0.18719_dp,  0.05702_dp, &
               -0.06080_dp, -0.10648_dp, -0.27805_dp, -0.38289_dp, &
               -0.29730_dp, -0.22220_dp, -0.20196_dp, -0.18109_dp, &
                0.00000_dp ]

      if ( abs(re - 100.0_dp) <= 0.5_dp ) then
         uref = u100
      else if ( abs(re - 1000.0_dp) <= 5.0_dp ) then
         uref = u1000
      else
         write(*,'(a)') '  Ghia reference: no table for this Re, skipped.'
         return
      end if

      write(*,'(a)') ''
      write(*,'(a,es10.3,a)') '--- Ghia et al. (1982) comparison, Re = ', re, ' ---'
      write(*,'(a)') '      y        u_num      u_ghia      error'
      l2num = 0.0_dp
      l2den = 0.0_dp
      do i = 1, NPTS
         xq = [ 0.5_dp, yg(i), 0.5_dp ]
         call probe_velocity( m, g, fld, xq, uq )
         err = uq(1) - uref(i)
         write(*,'(f9.4,2es12.5,es12.2)') yg(i), uq(1), uref(i), err
         ! NB: the two wall rows (y=0,1) are excluded from the L2 norm --
         ! the centroid-IDW probe cannot reproduce boundary values
         if ( i > 1 .and. i < NPTS ) then
            l2num = l2num + err*err
            l2den = l2den + uref(i)*uref(i)
         end if
      end do
      write(*,'(a,es12.5)') '  relative L2 error (centerline u_x): ', &
         sqrt( l2num / max( l2den, tiny(1.0_dp) ) )
   end subroutine ghia_compare

   !----------------------------------------------------------------------------
   ! Wall heat-flux / Nusselt report for natural-convection validation.
   !
   ! For every fixed-temperature (Dirichlet) boundary face -- either a zone-
   ! level fixed-T wall or a tbc_plane face -- the discrete diffusive heat
   ! rate into the owner cell, consistent with temperature_assembly, is
   !       dQ = k_eff * area / dw * (T_face - T_P)         [W]
   ! Faces are grouped by (zone, plane direction, plane coordinate); the
   ! integrated rate Q, area A and mean flux q=Q/A are printed.  When exactly
   ! two fixed-T planes on the same axis are present, the Nusselt number
   !       Nu = (|Q|/A) * L / (k_f * |T1-T2|)
   ! is reported (k_f = ctrl%k_cond, the fluid conductivity reference).
   !----------------------------------------------------------------------------
   subroutine thermal_wall_report( m, g, ctrl, bcs, fld )
      type(mesh_t),   intent(in) :: m
      type(geom_t),   intent(in) :: g
      type(ctrl_t),   intent(in) :: ctrl
      type(bc_t),     intent(in) :: bcs
      type(fields_t), intent(in) :: fld

      integer, parameter :: MAXG = 12

      integer  :: i, c0, gi, ipl, kd, ng, idir, ig
      real(dp) :: T_face, q_face, dw, k_eff, dQ
      real(dp) :: gzone(MAXG), gdir(MAXG), gcoord(MAXG), gT(MAXG)
      real(dp) :: gQ(MAXG), gA(MAXG)
      logical  :: is_neumann, found, is_ltne
      integer  :: n1, n2
      real(dp) :: Ldist, dT, Nu

      is_ltne = ( ctrl%thermal_model == 'ltne' )
      ng = 0
      gQ = 0.0_dp
      gA = 0.0_dp

      do i = 1, m%nfaces
         gi = bcs%fgrp(i)
         if ( gi == 0 ) cycle
         c0 = m%f(i)%c0
         call bc_face_T( bcs, i, g%xf(:,i), fld%T(c0), T_face, q_face, is_neumann )
         if ( is_neumann ) cycle          ! adiabatic / fixed-flux faces skipped

         ! identify the plane this face matched (0 = zone-level fixed T)
         idir  = 0
         gcoord_loop: do ipl = 1, bcs%gb(gi)%ntbc_plane
            if ( bcs%gb(gi)%tbc_pttype(ipl) /= 1 ) cycle
            if ( abs( g%xf(bcs%gb(gi)%tbc_pdir(ipl),i) &
                      - bcs%gb(gi)%tbc_pcoord(ipl) ) <= bcs%gb(gi)%tbc_tol ) then
               idir = bcs%gb(gi)%tbc_pdir(ipl)
               exit gcoord_loop
            end if
         end do gcoord_loop

         ! find / create group
         found = .false.
         do ig = 1, ng
            if ( nint(gzone(ig)) == bcs%gb(gi)%zone .and. &
                 nint(gdir(ig))  == idir ) then
               if ( idir == 0 .or. &
                    abs( gcoord(ig) - merge( bcs%gb(gi)%tbc_pcoord(ipl), &
                                             0.0_dp, idir /= 0 ) ) &
                        < epsilon(1.0_dp) ) then
                  found = .true.; exit
               end if
            end if
         end do
         if ( .not. found ) then
            ng = ng + 1
            if ( ng > MAXG ) cycle
            gzone(ng)  = real( bcs%gb(gi)%zone, dp )
            gdir(ng)   = real( idir, dp )
            if ( idir /= 0 ) then
               gcoord(ng) = bcs%gb(gi)%tbc_pcoord(ipl)
            else
               gcoord(ng) = 0.0_dp
            end if
            gT(ng) = T_face
         end if
         kd = ig

         ! owner-cell effective conductivity (same rule as temperature_assembly)
         if ( is_ltne ) then
            k_eff = fld%porosity(c0) * ctrl%k_cond
         else
            k_eff = fld%porosity(c0)*ctrl%k_cond &
                  + (1.0_dp - fld%porosity(c0)) * fld%k_s(c0)
         end if

         dw = norm2( g%xf(:,i) - g%xc(:,c0) )
         dQ = k_eff * g%area(i) / dw * ( T_face - fld%T(c0) )
         gQ(kd) = gQ(kd) + dQ
         gA(kd) = gA(kd) + g%area(i)
      end do

      if ( ng == 0 ) return

      write(*,'(a)') ''
      write(*,'(a)') '--- Fixed-temperature wall heat balance ---'
      write(*,'(a)') '   zone  dir   coord      T_wall        Q [W]          area [m2]      q [W/m2]'
      do ig = 1, ng
         write(*,'(i6,2x,i3,2x,f7.4,2x,f9.4,2x,es12.4,2x,es12.4,2x,es12.4)') &
            nint(gzone(ig)), nint(gdir(ig)), gcoord(ig), gT(ig), &
            gQ(ig), gA(ig), gQ(ig)/max(gA(ig),tiny(1.0_dp))
      end do

      ! Nusselt number when exactly two fixed-T planes share one axis
      if ( ng == 2 .and. nint(gdir(1)) /= 0 .and. &
           nint(gdir(1)) == nint(gdir(2)) ) then
         n1 = 1; n2 = 2
         if ( gT(n2) > gT(n1) ) then; n1 = 2; n2 = 1; end if   ! n1 = hot wall
         Ldist = abs( gcoord(n1) - gcoord(n2) )
         dT    = gT(n1) - gT(n2)
         Nu = ( abs(gQ(n1)) / max(gA(n1),tiny(1.0_dp)) ) * Ldist &
              / max( ctrl%k_cond * dT, tiny(1.0_dp) )
         write(*,'(a)') ''
         write(*,'(a,es12.5)') '  hot-wall Q (into domain): ', gQ(n1)
         write(*,'(a,es12.5)') '  cold-wall Q (into domain):', gQ(n2)
         write(*,'(a,es12.5,a,es12.5)') &
            '  L = ', Ldist, '   dT = ', dT
         write(*,'(a,es12.5)') '  Nusselt number Nu (fluid k): ', Nu
         write(*,'(a)') '  (pure conduction reference: Nu = 1)'
      end if

   end subroutine thermal_wall_report

   !============================================================================
   ! Stage 17: native (real cell type) 1:1 output
   !============================================================================

   !----------------------------------------------------------------------------
   ! Write the mesh 1:1: one VTK cell per mesh cell, real VTK type code, and no
   ! centroid points.  The canonical VTK node order is rebuilt from the face
   ! topology by native_cell_nodes, which also verifies it (every face implied
   ! by the rebuilt order must match an actual face of the cell, and the
   ! divergence-theorem volume must be positive and equal to the solver's cell
   ! volume).
   ! Any unmapped Fluent code (1/2/3/5 = 0-D/1-D/2-D entities, 0 = unknown or
   ! "mixed zone marker", anything else) is a hard error: those meshes are
   ! outside this 3-D solver's scope, and a zero-volume / degenerate cell must
   ! never be written (user decision, see docs/plan.md stage 17).
   !----------------------------------------------------------------------------
   subroutine vtk_write_native( m, conn, g, fld, filename, ier, bin )
      type(mesh_t),    intent(in)  :: m
      type(conn_t),    intent(in)  :: conn
      type(geom_t),    intent(in)  :: g
      type(fields_t),  intent(in)  :: fld
      character(len=*), intent(in) :: filename
      integer,         intent(out) :: ier
      logical,         intent(in)  :: bin

      integer :: c, j, ncd, nc, nnp, nconn, code
      logical :: has_solid, has_ts
      integer, allocatable :: nodes(:)
      integer, allocatable :: cell_conn(:), cell_off(:), cell_type(:)
      real(dp), allocatable :: coords(:,:)
      real(dp) :: volc, volg, vscale
      type(vtu_data), allocatable :: cdata(:)

      ier = 0

      ! ---- the mesh reader must have provided the per-cell element types ----
      if ( .not. allocated(m%ctype) ) then
         write(*,'(a)') ''
         write(*,'(a)') 'vtk_write (out_cells = native): ERROR - this mesh has ' // &
            'no per-cell element types (mesh_t%ctype unset), so the real cell ' // &
            'types are unknown; use `out_cells = tet_fan` (the default).'
         ier = 3
         return
      end if

      ! ---- every cell must be a mapped 3-D type -----------------------------
      do c = 1, m%ncells
         code = m%ctype(c)
         if ( .not. cell_type_supported( code ) ) then
            call native_type_error( code, c, conn%cf_ptr(c+1) - conn%cf_ptr(c) )
            ier = 3
            return
         end if
      end do

      ! ---- points: mesh nodes only (no centroid points) ---------------------
      allocate( coords(3, m%nnodes) )
      coords = m%x(:, 1:m%nnodes)

      ! ---- cells: one VTK cell per mesh cell --------------------------------
      ! connectivity length: 8 nodes per hexahedron, 6 per wedge
      nconn = 0
      do c = 1, m%ncells
         nconn = nconn + merge( 8, 6, m%ctype(c) == 4 )
      end do
      allocate( cell_conn(nconn), cell_off(m%ncells), cell_type(m%ncells) )
      allocate( nodes(8) )

      nc = 0
      do c = 1, m%ncells
         call native_cell_nodes( m, conn, g, c, nodes, nnp, volc, ier )
         if ( ier /= 0 ) return
         ! the rebuilt cell must reproduce the solver's own volume: this catches
         ! a wrong outward orientation as well as a wrong node order
         volg   = g%vol(c)
         vscale = max( abs(volg), tiny(1.0_dp) )
         if ( abs( volc - volg ) > 1.0e-8_dp * vscale ) then
            write(*,'(a)') ''
            write(*,'(a)') 'vtk_write (out_cells = native): ERROR - rebuilt cell ' // &
               'volume /= computed cell volume (inverted or degenerate cell, ' // &
               'or wrong node ordering).'
            write(*,'(a,i0)')     '  cell id       : ', c
            write(*,'(a,i0)')     '  Fluent code   : ', m%ctype(c)
            write(*,'(a,es16.8)') '  rebuilt volume: ', volc
            write(*,'(a,es16.8)') '  solver volume : ', volg
            ier = 4
            return
         end if
         cell_conn(nc+1:nc+nnp) = nodes(1:nnp)
         nc = nc + nnp
         cell_off(c)  = nc
         cell_type(c) = fluent_to_vtk( m%ctype(c) )
      end do

      ! ---- cell data: one value per mesh cell (no sub-cell replication) -----
      !      same variables / same conditions as the fan path
      has_solid = .false.
      do c = 1, m%ncells
         if ( is_solid_cell( m, fld, c ) ) then
            has_solid = .true.
            exit
         end if
      end do
      has_ts = ( maxval( fld%h_sf * fld%a_sf ) > 0.0_dp .or. has_solid )

      ncd = 3 + merge(1,0,has_ts) + merge(1,0,has_solid)
      allocate( cdata(ncd) )
      cdata(1)%name  = 'pressure';    cdata(1)%ncomp = 1
      cdata(2)%name  = 'velocity';    cdata(2)%ncomp = 3
      cdata(3)%name  = 'temperature'; cdata(3)%ncomp = 1
      if ( has_ts )    cdata(4)%name   = 'temperature_solid'
      if ( has_solid ) cdata(ncd)%name = 'heat_source_solid'
      do j = 1, ncd
         allocate( cdata(j)%val(m%ncells, cdata(j)%ncomp) )
      end do
      cdata(1)%val(:,1)   = fld%p(1:m%ncells)
      cdata(2)%val(:,1:3) = transpose( fld%u(:, 1:m%ncells) )
      cdata(3)%val(:,1)   = fld%T(1:m%ncells)
      if ( has_ts )    cdata(4)%val(:,1)   = fld%T_s(1:m%ncells)
      if ( has_solid ) cdata(ncd)%val(:,1) = fld%q_src(1:m%ncells)

      ! ---- write the file through the generic XML writer --------------------
      call vtu_write( filename, coords, cell_conn, cell_off, cell_type, &
                      cell_data = cdata, binary = bin, ier = ier )
      if ( ier /= 0 ) return

      if ( bin ) then
         write(*,'(a)') '  wrote VTU file : ' // trim(filename) // &
                        '  (native cells, appended raw binary, Float32)'
      else
         write(*,'(a)') '  wrote VTU file : ' // trim(filename) // &
                        '  (native cells, 1 cell per mesh cell)'
      end if
   end subroutine vtk_write_native

   !----------------------------------------------------------------------------
   ! Canonical VTK node order of ONE mesh cell, rebuilt from its face topology,
   ! plus its divergence-theorem volume volc (> 0 required).
   !
   ! The face loops are oriented OUTWARD first: the stored .cas order is NOT
   ! guaranteed to be inward for c0 -- that is the tetrahedral fan's assumption,
   ! and it is what makes tools/validate_vtu.py --vtk report non-positive
   ! tetrahedra.  The outward direction comes from g%sf, which
   ! compute_geometry orients as c0 -> c1.
   !
   ! VTK orderings (ground truth: external/sample_code/unMesh/main.f90:138-171
   ! and .dbg/vtu_mixed_test.f90, both confirmed positive by the VTK reader in
   ! .dbg/chk_mixed.py):
   !   hexahedron (12): nodes 1-2-3-4 form one quad loop and 5-6-7-8 the
   !     opposite one, node i+4 "above" node i; both loops are traversed with
   !     their right-hand-rule normal along the stacking direction, so the
   !     loop of nodes 1-4 is the REVERSE of its outward loop;
   !   wedge (13): nodes 1-2-3 form one triangle and 4-5-6 the other, node i+3
   !     "above" node i, with the same RHR-along-the-stack rule.
   ! Verification: every face implied by the rebuilt order must be an actual
   ! face of the cell, and the volume must be positive (the caller additionally
   ! compares it with the solver's own cell volume).
   !----------------------------------------------------------------------------
   subroutine native_cell_nodes( m, conn, g, c, nodes, nnp, volc, ier )
      type(mesh_t), intent(in)  :: m
      type(conn_t), intent(in)  :: conn
      type(geom_t), intent(in)  :: g
      integer,      intent(in)  :: c
      integer,      intent(out) :: nodes(*)
      integer,      intent(out) :: nnp
      real(dp),     intent(out) :: volc
      integer,      intent(out) :: ier

      integer, parameter :: MAXF = 8
      integer  :: nf, k, f, i, nn, code, reason
      integer  :: ol(8, MAXF)      ! outward-oriented face loops
      integer  :: nno(MAXF)        ! number of nodes per face
      real(dp) :: vol

      ier = 0; nnp = 0; volc = 0.0_dp; reason = 0

      ! ---- faces of this cell, each oriented outward ------------------------
      nf = conn%cf_ptr(c+1) - conn%cf_ptr(c)
      if ( nf < 4 .or. nf > MAXF ) then
         call native_cell_error( 1, c, m, ier )
         return
      end if
      do k = 1, nf
         f  = conn%cf( conn%cf_ptr(c) + k )
         nn = m%f(f)%nn
         if ( nn < 3 .or. nn > 8 ) then
            call native_cell_error( 2, c, m, ier )
            return
         end if
         call orient_face_outward( m, g, f, c, ol(:,k), nno(k), reason )
         if ( reason /= 0 ) then
            call native_cell_error( reason, c, m, ier )
            return
         end if
      end do

      ! ---- volume from the outward-oriented faces (divergence theorem) ------
      vol = 0.0_dp
      do k = 1, nf
         do i = 2, nno(k) - 1
            vol = vol + tri_volume( m%x(:,ol(1,k)), m%x(:,ol(i,k)), &
                                    m%x(:,ol(i+1,k)) )
         end do
      end do
      if ( vol <= 0.0_dp ) then
         call native_cell_error( 4, c, m, ier )
         return
      end if
      volc = vol

      ! ---- canonical order + face verification ------------------------------
      code = m%ctype(c)
      select case ( code )
      case ( 4 )
         call hex_canonical_nodes( ol, nno, nf, nodes, reason )
         nnp = 8
      case ( 6 )
         call wedge_canonical_nodes( ol, nno, nf, nodes, reason )
         nnp = 6
      case default
         call native_cell_error( 5, c, m, ier )
         return
      end select
      if ( reason /= 0 ) then
         call native_cell_error( reason, c, m, ier )
         return
      end if
      if ( .not. implied_faces_ok( ol, nno, nf, nodes, code ) ) then
         call native_cell_error( 8, c, m, ier )
         return
      end if
   end subroutine native_cell_nodes

   !----------------------------------------------------------------------------
   ! Canonical VTK node order of a hexahedron (code 4) from its 6 outward-
   ! oriented quad loops.  reason = 0 on success, 6 = "not a hexahedron".
   !----------------------------------------------------------------------------
   subroutine hex_canonical_nodes( ol, nno, nf, nodes, reason )
      integer, intent(in)  :: ol(8,*), nno(*), nf
      integer, intent(out) :: nodes(*)
      integer, intent(out) :: reason

      integer :: k, i, r, ftop
      logical :: good

      reason = 0
      nodes(1:8) = 0
      if ( nf /= 6 ) then
         reason = 6; return
      end if
      do k = 1, nf
         if ( nno(k) /= 4 ) then
            reason = 6; return
         end if
      end do

      ! the face opposite the first quad shares no node with it
      ftop = 0
      do k = 2, nf
         good = .true.
         do i = 1, 4
            if ( any( ol(1:4,k) == ol(i,1) ) ) then
               good = .false.
               exit
            end if
         end do
         if ( good ) then
            ftop = k
            exit
         end if
      end do
      if ( ftop == 0 ) then
         reason = 6; return
      end if

      ! nodes 1-4 = REVERSE of the bottom (outward) loop: its right-hand-rule
      ! normal must point along the stack direction (the VTK convention)
      nodes(1) = ol(1,1)
      nodes(2) = ol(4,1)
      nodes(3) = ol(3,1)
      nodes(4) = ol(2,1)

      ! ---- nodes 5-8: rotate the top outward loop until ALL SIX faces implied
      !      by the numbering are actual faces of the cell (definitive test: a
      !      rotation that paired diagonal vertices would produce a face that
      !      does not exist)
      do r = 1, 4
         do i = 1, 4
            nodes(4+i) = ol( mod(r+i-2,4)+1, ftop )
         end do
         if ( implied_faces_ok( ol, nno, nf, nodes, 4 ) ) return
      end do
      reason = 6
   end subroutine hex_canonical_nodes

   !----------------------------------------------------------------------------
   ! Canonical VTK node order of a wedge / triangular prism (code 6) from its
   ! 2 triangle + 3 quad outward-oriented loops.  reason = 0 on success,
   ! 7 = "not a wedge".
   !----------------------------------------------------------------------------
   subroutine wedge_canonical_nodes( ol, nno, nf, nodes, reason )
      integer, intent(in)  :: ol(8,*), nno(*), nf
      integer, intent(out) :: nodes(*)
      integer, intent(out) :: reason

      integer :: k, i, r, ntri, nquad, f1, f2

      reason = 0
      nodes(1:6) = 0
      if ( nf /= 5 ) then
         reason = 7; return
      end if
      ntri = 0; nquad = 0; f1 = 0; f2 = 0
      do k = 1, nf
         if ( nno(k) == 3 ) then
            ntri = ntri + 1
            if ( f1 == 0 ) then
               f1 = k
            else
               f2 = k
            end if
         else if ( nno(k) == 4 ) then
            nquad = nquad + 1
         else
            reason = 7; return
         end if
      end do
      if ( ntri /= 2 .or. nquad /= 3 ) then
         reason = 7; return
      end if

      ! nodes 1-3 = REVERSE of the first triangle's outward loop (RHR normal
      ! along the stack direction towards the opposite triangle)
      nodes(1) = ol(1,f1)
      nodes(2) = ol(3,f1)
      nodes(3) = ol(2,f1)

      ! ---- nodes 4-6: rotate the other triangle's outward loop until all five
      !      implied faces are actual faces of the cell
      do r = 1, 3
         do i = 1, 3
            nodes(3+i) = ol( mod(r+i-2,3)+1, f2 )
         end do
         if ( implied_faces_ok( ol, nno, nf, nodes, 6 ) ) return
      end do
      reason = 7
   end subroutine wedge_canonical_nodes

   !----------------------------------------------------------------------------
   ! Are all faces IMPLIED by this canonical node numbering actual faces of the
   ! cell?  (node-set comparison against the outward-oriented loops).  This is
   ! the definitive validity test for the rebuilt ordering.
   !----------------------------------------------------------------------------
   logical function implied_faces_ok( ol, nno, nf, nodes, code )
      integer, intent(in) :: ol(8,*), nno(*), nf, nodes(*), code

      integer :: imp(4,6), nimp(6), nif, k, j, i, cnt, want(4)

      if ( code == 4 ) then                       ! hexahedron
         nif = 6
         nimp(1:6)        = 4
         imp(1:4,1) = [1,2,3,4]                    ! bottom
         imp(1:4,2) = [5,6,7,8]                    ! top
         imp(1:4,3) = [1,2,6,5]
         imp(1:4,4) = [2,3,7,6]
         imp(1:4,5) = [3,4,8,7]
         imp(1:4,6) = [4,1,5,8]
      else                                        ! wedge (triangular prism)
         nif = 5
         nimp(1:5)        = [3,3,4,4,4]
         imp(1:3,1) = [1,2,3]
         imp(1:3,2) = [4,5,6]
         imp(1:4,3) = [1,2,5,4]
         imp(1:4,4) = [2,3,6,5]
         imp(1:4,5) = [3,1,4,6]
      end if

      implied_faces_ok = .false.
      do k = 1, nif
         do i = 1, nimp(k)
            want(i) = nodes( imp(i,k) )
         end do
         do j = 1, nf
            if ( nno(j) /= nimp(k) ) cycle
            cnt = 0
            do i = 1, nimp(k)
               if ( any( ol(1:nno(j),j) == want(i) ) ) cnt = cnt + 1
            end do
            if ( cnt == nimp(k) ) exit             ! this face exists
         end do
         if ( j > nf ) return                      ! not found -> false
      end do
      implied_faces_ok = .true.
   end function implied_faces_ok

   !----------------------------------------------------------------------------
   ! Orient one face loop outward for cell c.  The stored .cas order is not
   ! relied on (the fan's "inward for c0" assumption does not hold everywhere):
   ! the outward direction is taken from g%sf, which compute_geometry orients
   ! as c0 -> c1.  reason = 3 when the face is degenerate (zero area) or its
   ! normal is orthogonal to the outward direction -- a degenerate cell, which
   ! this solver must not write out.
   !----------------------------------------------------------------------------
   subroutine orient_face_outward( m, g, f, c, loop, nn, reason )
      type(mesh_t), intent(in)  :: m
      type(geom_t), intent(in)  :: g
      integer,      intent(in)  :: f, c
      integer,      intent(out) :: loop(:), nn
      integer,      intent(out) :: reason

      integer  :: j, jn
      real(dp) :: s(3), out(3), ds

      reason = 0
      nn = m%f(f)%nn
      loop(1:nn) = m%f(f)%nodes(1:nn)

      ! outward direction of this face as seen from cell c
      out = g%sf(:,f)
      if ( m%f(f)%c1 == c ) out = -out

      ! Newell normal of the stored loop (same construction as face_props)
      s = 0.0_dp
      do j = 1, nn
         jn = mod(j,nn) + 1
         s = s + 0.5_dp * cross3( m%x(:,loop(j)), m%x(:,loop(jn)) )
      end do
      if ( norm2(s) <= 0.0_dp ) then
         reason = 3
         return
      end if
      ds = dot_product( s, out )
      if ( abs(ds) <= 1.0e-12_dp * norm2(s) * norm2(out) ) then
         reason = 3
         return
      end if
      if ( ds < 0.0_dp ) loop(1:nn) = loop(nn:1:-1)
   end subroutine orient_face_outward

   !----------------------------------------------------------------------------
   ! cross product of two 3-vectors
   !----------------------------------------------------------------------------
   pure function cross3( a, b ) result( c )
      real(dp), intent(in) :: a(3), b(3)
      real(dp)             :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross3

   !----------------------------------------------------------------------------
   ! Signed volume of the tetrahedron (origin, a, b, c): positive when (a,b,c)
   ! is right-handed.  Used with outward-oriented face loops to get the
   ! divergence-theorem volume of a cell.
   !----------------------------------------------------------------------------
   pure function tri_volume( a, b, c ) result( v )
      real(dp), intent(in) :: a(3), b(3), c(3)
      real(dp)             :: v

      v = dot_product( a, cross3( b, c ) ) / 6.0_dp
   end function tri_volume

   !----------------------------------------------------------------------------
   ! Error report for an unsupported Fluent cell code on the native path.
   !----------------------------------------------------------------------------
   subroutine native_type_error( code, c, nfaces )
      integer, intent(in) :: code, c, nfaces

      write(*,'(a)') ''
      write(*,'(a)') 'vtk_write (out_cells = native): ERROR - unsupported cell type'
      write(*,'(a,i0)') '  cell id       : ', c
      write(*,'(a,i0)') '  Fluent code   : ', code
      write(*,'(a,a)')  '  meaning       : ', trim(fluent_code_name(code))
      write(*,'(a,i0)') '  faces on cell : ', nfaces
      write(*,'(a)')    '  Only the 3-D cell types this solver can produce are ' // &
         'mapped (4 = hexahedron, 6 = wedge).  0-D/1-D/2-D entities (1/2/3/5) ' // &
         'and degenerate (zero-volume) cells are not allowed.'
      write(*,'(a)')    '  Use `out_cells = tet_fan` (the default) for this mesh.'
   end subroutine native_type_error

   !----------------------------------------------------------------------------
   ! Error report for a cell whose canonical VTK order cannot be rebuilt.
   !----------------------------------------------------------------------------
   subroutine native_cell_error( reason, c, m, ier )
      integer,      intent(in)  :: reason, c
      type(mesh_t), intent(in)  :: m
      integer,      intent(out) :: ier

      write(*,'(a)') ''
      write(*,'(a)') 'vtk_write (out_cells = native): ERROR - cannot rebuild the ' // &
                     'canonical VTK cell of this mesh cell.'
      write(*,'(a,i0)') '  cell id     : ', c
      write(*,'(a,i0)') '  Fluent code : ', m%ctype(c)
      write(*,'(a,a)')  '  meaning     : ', trim(fluent_code_name( m%ctype(c) ))
      select case ( reason )
      case ( 1 )
         write(*,'(a)') '  reason      : implausible number of faces on the cell'
      case ( 2 )
         write(*,'(a)') '  reason      : face with fewer than 3 or more than 8 nodes'
      case ( 3 )
         write(*,'(a)') '  reason      : zero-area face, or face normal orthogonal ' // &
                        'to the c0 -> c1 direction (degenerate cell)'
      case ( 4 )
         write(*,'(a)') '  reason      : non-positive cell volume'
      case ( 5 )
         write(*,'(a)') '  reason      : Fluent code not mapped to a VTK cell type'
      case ( 6 )
         write(*,'(a)') '  reason      : not a hexahedron (6 quads, 2 opposite needed)'
      case ( 7 )
         write(*,'(a)') '  reason      : not a wedge (2 triangles + 3 quads needed)'
      case ( 8 )
         write(*,'(a)') '  reason      : the rebuilt node order does not reproduce ' // &
                        'the cell faces'
      case default
         write(*,'(a)') '  reason      : unspecified'
      end select
      write(*,'(a)') '  (zero-volume/degenerate cells and unmapped cell types are ' // &
                     'not allowed -- see docs/plan.md stage 17)'
      ier = 20 + reason
   end subroutine native_cell_error

end module mod_uns_output
