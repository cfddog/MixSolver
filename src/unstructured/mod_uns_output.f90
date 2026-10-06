!===============================================================================
! mod_output.f90 -- Result output and benchmark comparison (phase 4)
!
! Contents:
!   vtk_write     -- write the solution as an ASCII .vtu (VTK unstructured
!                    grid, XML). Every cell is decomposed into tetrahedra by
!                    fanning each (outward-oriented) cell face from its first
!                    node and connecting the triangles to the cell centroid,
!                    which is appended as an extra point. VTK_TETRA cells are
!                    universally supported -- unlike VTK_POLYHEDRON, whose
!                    legacy face-stream encoding modern ParaView (VTK 9 based)
!                    silently drops, leaving an empty-looking dataset.
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
!===============================================================================
module mod_uns_output
   use mod_precision, only: dp, ip, pi
   use mod_uns_mesh
   use mod_uns_connectivity
   use mod_uns_geometry
   use mod_uns_fields
   use mod_uns_control
   use mod_uns_bc
   use iso_c_binding
   implicit none
   private
   public :: vtk_write, tecplot_write, probe_velocity, ghia_compare, &
             thermal_wall_report

   integer, parameter :: VTK_TETRA = 10

   !----------------------------------------------------------------------------
   ! Explicit interfaces to the TecIO (libtecio) classic 142 Fortran API, used
   ! to write a cell-centred FE-tetrahedral binary PLT file.  All scalar
   ! arguments are passed by reference (classic Fortran convention); strings
   ! are NUL-terminated.  These bind to the underscore gfortran entry points.
   !----------------------------------------------------------------------------
   interface
      integer(c_int32_t) function tecini142( Title, Variables, FName, ScratchDir, &
           FileFormat, FileType, Debug, VIsDouble ) bind(c, name="tecini142_")
         use iso_c_binding
         character(c_char), intent(in)  :: Title(*)
         character(c_char), intent(in)  :: Variables(*)
         character(c_char), intent(in)  :: FName(*)
         character(c_char), intent(in)  :: ScratchDir(*)
         integer(c_int32_t), intent(in) :: FileFormat
         integer(c_int32_t), intent(in) :: FileType
         integer(c_int32_t), intent(in) :: Debug
         integer(c_int32_t), intent(in) :: VIsDouble
      end function tecini142

      integer(c_int32_t) function teczne142( ZoneTitle, ZoneType, &
           IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces, &
           ICellMax, JCellMax, KCellMax, SolutionTime, StrandID, ParentZone, &
           IsBlock, NumFaceConnections, FaceNeighborMode, TotalNumFaceNodes, &
           NumConnectedBoundaryFaces, TotalNumBoundaryConnections, &
           PassiveVarList, ValueLocation, ShareVarFromZone, &
           ShareConnectivityFromZone ) bind(c, name="teczne142_")
         use iso_c_binding
         character(c_char), intent(in)  :: ZoneTitle(*)
         integer(c_int32_t), intent(in) :: ZoneType
         integer(c_int32_t), intent(in) :: IMxOrNumPts
         integer(c_int32_t), intent(in) :: JMxOrNumElements
         integer(c_int32_t), intent(in) :: KMxOrNumFaces
         integer(c_int32_t), intent(in) :: ICellMax
         integer(c_int32_t), intent(in) :: JCellMax
         integer(c_int32_t), intent(in) :: KCellMax
         real(c_double), intent(in)     :: SolutionTime
         integer(c_int32_t), intent(in) :: StrandID
         integer(c_int32_t), intent(in) :: ParentZone
         integer(c_int32_t), intent(in) :: IsBlock
         integer(c_int32_t), intent(in) :: NumFaceConnections
         integer(c_int32_t), intent(in) :: FaceNeighborMode
         integer(c_int32_t), intent(in) :: TotalNumFaceNodes
         integer(c_int32_t), intent(in) :: NumConnectedBoundaryFaces
         integer(c_int32_t), intent(in) :: TotalNumBoundaryConnections
         integer(c_int32_t), intent(in) :: PassiveVarList(*)
         integer(c_int32_t), intent(in) :: ValueLocation(*)
         integer(c_int32_t), intent(in) :: ShareVarFromZone(*)
         integer(c_int32_t), intent(in) :: ShareConnectivityFromZone
      end function teczne142

      integer(c_int32_t) function tecdatd142( N, FieldData ) &
           bind(c, name="tecdatd142_")
         use iso_c_binding
         integer(c_int32_t), intent(in) :: N
         real(c_double), intent(in)     :: FieldData(*)
      end function tecdatd142

      integer(c_int32_t) function tecnod142( NData ) &
           bind(c, name="tecnod142_")
         use iso_c_binding
         integer(c_int32_t), intent(in) :: NData(*)
      end function tecnod142

      integer(c_int32_t) function tecend142() bind(c, name="tecend142_")
         use iso_c_binding
      end function tecend142
   end interface

contains

   !----------------------------------------------------------------------------
   ! Write ASCII .vtu (XML, inline ascii data), cells tetrahedralized.
   !----------------------------------------------------------------------------
   subroutine vtk_write( m, conn, g, fld, filename, ier )
      type(mesh_t),    intent(in)  :: m
      type(conn_t),    intent(in)  :: conn
      type(geom_t),    intent(in)  :: g
      type(fields_t),  intent(in)  :: fld
      character(len=*), intent(in) :: filename
      integer,         intent(out) :: ier

      integer  :: iu, ios, i, j, k, c, f, nn, nsub, ntets, cen_id
      integer  :: tri(3)

      ier = 0

      ! total number of tetrahedra: each face contributes (nn-2) triangles
      ! NB: cf is 0-offset CSR -- entries of cell c are cf(cf_ptr(c)+1 : cf_ptr(c+1))
      ntets = 0
      do c = 1, m%ncells
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            f = conn%cf(k)
            ntets = ntets + m%f(f)%nn - 2
         end do
      end do

      open( newunit=iu, file=filename, status='replace', action='write', &
            iostat=ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'VTK ERROR: cannot open ' // trim(filename)
         ier = 1
         return
      end if

      write(iu,'(a)') '<?xml version="1.0"?>'
      write(iu,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" '// &
                      'byte_order="LittleEndian">'
      write(iu,'(a)') '  <UnstructuredGrid>'
      write(iu,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', &
         m%nnodes + m%ncells, '" NumberOfCells="', ntets, '">'

      ! ---- points: mesh nodes, then one centroid point per cell ----------------
      write(iu,'(a)') '      <Points>'
      write(iu,'(a)') '        <DataArray type="Float64" Name="Coordinates" '// &
                      'NumberOfComponents="3" format="ascii">'
      do i = 1, m%nnodes
         write(iu,'(3(es17.8,1x))') m%x(1,i), m%x(2,i), m%x(3,i)
      end do
      do c = 1, m%ncells
         write(iu,'(3(es17.8,1x))') g%xc(1,c), g%xc(2,c), g%xc(3,c)
      end do
      write(iu,'(a)') '        </DataArray>'
      write(iu,'(a)') '      </Points>'

      ! ---- cells: tetrahedra from face fan + cell centroid ---------------------
      write(iu,'(a)') '      <Cells>'
      write(iu,'(a)') '        <DataArray type="Int32" Name="connectivity" '// &
                      'format="ascii">'
      do c = 1, m%ncells
         cen_id = m%nnodes + c - 1            ! 0-based id of centroid point
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
                  close(iu)
                  return
               end if
               write(iu,'(5(i0,1x))') 4, tri(1)-1, tri(2)-1, tri(3)-1, cen_id
            end do
         end do
      end do
      write(iu,'(a)') '        </DataArray>'

      write(iu,'(a)') '        <DataArray type="Int32" Name="offsets" '// &
                      'format="ascii">'
      do c = 1, ntets
         write(iu,'(i0,1x)') 5*c
      end do
      write(iu,'(a)') '        </DataArray>'

      write(iu,'(a)') '        <DataArray type="UInt8" Name="types" '// &
                      'format="ascii">'
      do c = 1, ntets
         write(iu,'(i0,1x)') VTK_TETRA
      end do
      write(iu,'(a)') '        </DataArray>'
      write(iu,'(a)') '      </Cells>'

      ! ---- cell data (original cell values replicated on their sub-tets) ------
      write(iu,'(a)') '      <CellData>'
      write(iu,'(a)') '        <DataArray type="Float64" Name="pressure" '// &
                      'format="ascii">'
      do c = 1, m%ncells
         nsub = 0
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            nsub = nsub + m%f(conn%cf(k))%nn - 2
         end do
         do j = 1, nsub
            write(iu,'(es17.8,1x)') fld%p(c)
         end do
      end do
      write(iu,'(a)') '        </DataArray>'
      write(iu,'(a)') '        <DataArray type="Float64" Name="velocity" '// &
                      'NumberOfComponents="3" format="ascii">'
      do c = 1, m%ncells
         nsub = 0
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            nsub = nsub + m%f(conn%cf(k))%nn - 2
         end do
         do j = 1, nsub
            write(iu,'(3(es17.8,1x))') fld%u(:,c)
         end do
      end do
      write(iu,'(a)') '        </DataArray>'
      write(iu,'(a)') '        <DataArray type="Float64" Name="temperature" '// &
                      'format="ascii">'
      do c = 1, m%ncells
         nsub = 0
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            nsub = nsub + m%f(conn%cf(k))%nn - 2
         end do
         do j = 1, nsub
            write(iu,'(es17.8,1x)') fld%T(c)
         end do
      end do
      write(iu,'(a)') '        </DataArray>'
      ! solid temperature (LTNE): written only when fluid-solid heat exchange
      ! is active somewhere in the domain (h_sf*a_sf > 0).
      if ( maxval( fld%h_sf * fld%a_sf ) > 0.0_dp ) then
         write(iu,'(a)') '        <DataArray type="Float64" '// &
                         'Name="temperature_solid" format="ascii">'
         do c = 1, m%ncells
            nsub = 0
            do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
               nsub = nsub + m%f(conn%cf(k))%nn - 2
            end do
            do j = 1, nsub
               write(iu,'(es17.8,1x)') fld%T_s(c)
            end do
         end do
         write(iu,'(a)') '        </DataArray>'
      end if
      write(iu,'(a)') '      </CellData>'

      write(iu,'(a)') '    </Piece>'
      write(iu,'(a)') '  </UnstructuredGrid>'
      write(iu,'(a)') '</VTKFile>'
      close(iu)

      write(*,'(a)') '  wrote VTU file : ' // trim(filename)
   end subroutine vtk_write

   !----------------------------------------------------------------------------
   ! Write a binary Tecplot PLT file.  Uses the SAME tetrahedralised geometry
   ! as vtk_write (mesh nodes + one centroid point per cell; face-fan tets) so
   ! the two formats are guaranteed to represent the identical mesh.
   !
   ! Variables: X, Y, Z are node-centred; rho, u, v, w, p, T are cell-centred
   ! (one value per sub-tet, the parent-cell value replicated).  rho is the
   ! uniform constant ctrl%rho (incompressible flow); the other fields come
   ! from the collocated fld container.
   !----------------------------------------------------------------------------
   subroutine tecplot_write( m, conn, g, ctrl, fld, filename, ier )
      type(mesh_t),     intent(in)  :: m
      type(conn_t),     intent(in)  :: conn
      type(geom_t),     intent(in)  :: g
      type(ctrl_t),     intent(in)  :: ctrl
      type(fields_t),   intent(in)  :: fld
      character(len=*), intent(in)  :: filename
      integer,          intent(out) :: ier

      ! TecIO enum values (see TECIO.h / tecio.for): ZoneType FE-Tetra = 4.
      ! File format PLT, file type Full and the static-strand flags are all 0;
      ! variables are written double (VIsDouble = 1).
      !
      ! NB: the prebuilt libtecio.a bundled in lib/tecplot (UCNS3D build) has
      ! INVERTED ValueLocation semantics for the classic 142 API, verified
      ! empirically with a minimal FE-tet reproducer: array entry 1 selects a
      ! NODE-centred variable, entry 0 a CELL-centred one (opposite of the
      ! Tecplot documentation; passing NULL gives the all-nodal default).
      integer(c_int32_t), parameter :: ZT_FETETRA = 4
      integer(c_int32_t), parameter :: VL_NODE    = 1
      integer(c_int32_t), parameter :: VL_CELL    = 0

      integer, parameter :: NVAR = 9    ! X Y Z rho u v w p T

      integer  :: i, j, k, c, f, nn, ntets, npts, it, cen_id
      integer  :: tri(3)
      integer(c_int32_t) :: res
      integer(c_int32_t) :: iFmt, iFtype, iDbg, iDbl, iZero, ztype
      integer(c_int32_t) :: npts_i, ntets_i, strand, parent, block_f
      real(c_double)     :: soltime
      integer(c_int32_t) :: valueLoc(NVAR), passive(NVAR), shareVar(NVAR)
      integer(c_int32_t), allocatable :: nodmap(:)
      real(dp), allocatable :: px(:), py(:), pz(:)
      real(dp), allocatable :: vrho(:), vu(:), vv(:), vw(:), vp(:), vT(:)

      ier = 0

      ! total number of tetrahedra: each face contributes (nn-2) triangles
      ntets = 0
      do c = 1, m%ncells
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            f = conn%cf(k)
            ntets = ntets + m%f(f)%nn - 2
         end do
      end do
      npts = m%nnodes + m%ncells

      allocate( px(npts), py(npts), pz(npts) )
      allocate( nodmap(4*ntets) )
      allocate( vrho(ntets), vu(ntets), vv(ntets), vw(ntets), vp(ntets), &
                vT(ntets) )

      ! ---- node coordinates: mesh nodes first, then one centroid per cell ----
      do i = 1, m%nnodes
         px(i) = m%x(1,i)
         py(i) = m%x(2,i)
         pz(i) = m%x(3,i)
      end do
      do c = 1, m%ncells
         px(m%nnodes+c) = g%xc(1,c)
         py(m%nnodes+c) = g%xc(2,c)
         pz(m%nnodes+c) = g%xc(3,c)
      end do

      ! ---- node map (1-based) + cell-centred field per sub-tet ----------------
      it = 0
      do c = 1, m%ncells
         cen_id = m%nnodes + c            ! 1-based id of centroid point
         do k = conn%cf_ptr(c)+1, conn%cf_ptr(c+1)
            f  = conn%cf(k)
            nn = m%f(f)%nn
            do j = 2, nn-1                ! fan from face node 1
               if ( m%f(f)%c0 == c ) then
                  ! NB: .cas face node order is INWARD for c0, so the stored
                  ! fan (n1,nj,nj+1) + centroid apex gives positive volume
                  tri = [ m%f(f)%nodes(1), m%f(f)%nodes(j), m%f(f)%nodes(j+1) ]
               else if ( m%f(f)%c1 == c ) then
                  tri = [ m%f(f)%nodes(1), m%f(f)%nodes(j+1), m%f(f)%nodes(j) ]
               else
                  write(*,'(a,i0,a,i0)') 'TECPLOT ERROR: face ', f, &
                     ' not adjacent to cell ', c
                  ier = 2
                  go to 900
               end if
               it = it + 1
               nodmap(4*it-3:4*it) = [ tri(1), tri(2), tri(3), cen_id ]
               vrho(it) = ctrl%rho
               vu(it)   = fld%u(1,c)
               vv(it)   = fld%u(2,c)
               vw(it)   = fld%u(3,c)
               vp(it)   = fld%p(c)
               vT(it)   = fld%T(c)
            end do
         end do
      end do

      ! ---- initialize the PLT dataset -----------------------------------------
      iFmt   = 0      ! 0 = PLT binary
      iFtype = 0      ! 0 = Full (grid + solution)
      iDbg   = 0
      iDbl   = 1      ! all variables double precision
      res = tecini142( 'UNSSolver result'//c_null_char, &
                       'X Y Z rho u v w p T'//c_null_char, &
                       trim(filename)//c_null_char, '.'//c_null_char, &
                       iFmt, iFtype, iDbg, iDbl )
      if ( res /= 0 ) then
         write(*,'(a,a,a)') 'TECPLOT ERROR: cannot open ', &
            trim(filename), ' (tecini142)'
         ier = 1
         go to 900
      end if

      ! ---- create the FE-tetrahedral zone -------------------------------------
      npts_i  = npts
      ntets_i = ntets
      iZero   = 0
      ztype   = ZT_FETETRA
      soltime = 0.0_c_double
      strand  = 0          ! static (non-transient) zone
      parent  = 0
      block_f = 1          ! block data packing
      passive(:)  = 0
      shareVar(:) = 0
      valueLoc    = [ VL_NODE, VL_NODE, VL_NODE, VL_CELL, VL_CELL, VL_CELL, &
                      VL_CELL, VL_CELL, VL_CELL ]
      res = teczne142( 'flow'//c_null_char, ztype, npts_i, ntets_i, iZero, &
                       iZero, iZero, iZero, soltime, strand, parent, block_f, &
                       iZero, iZero, iZero, iZero, iZero, &
                       passive, valueLoc, shareVar, iZero )
      if ( res /= 0 ) then
         write(*,'(a,i0)') 'TECPLOT ERROR: teczne142 failed, code ', res
         ier = 3
         go to 900
      end if

      ! ---- field data, block mode: coordinates at nodes, rho..T at cells.
      !      Classic API requires ALL variable data BEFORE the connectivity. --
      res = tecdatd142( npts_i,  px )
      if ( res == 0 ) res = tecdatd142( npts_i,  py )
      if ( res == 0 ) res = tecdatd142( npts_i,  pz )
      if ( res == 0 ) res = tecdatd142( ntets_i, vrho )
      if ( res == 0 ) res = tecdatd142( ntets_i, vu )
      if ( res == 0 ) res = tecdatd142( ntets_i, vv )
      if ( res == 0 ) res = tecdatd142( ntets_i, vw )
      if ( res == 0 ) res = tecdatd142( ntets_i, vp )
      if ( res == 0 ) res = tecdatd142( ntets_i, vT )
      if ( res /= 0 ) then
         write(*,'(a,i0)') 'TECPLOT ERROR: tecdatd142 failed, code ', res
         ier = 5
         go to 900
      end if

      ! ---- connectivity last: whole 1-based node map, 4 nodes per tet -------
      res = tecnod142( nodmap )
      if ( res /= 0 ) then
         write(*,'(a,i0)') 'TECPLOT ERROR: tecnod142 failed, code ', res
         ier = 4
         go to 900
      end if

      res = tecend142()
      if ( res /= 0 ) then
         write(*,'(a,i0)') 'TECPLOT ERROR: tecend142 failed, code ', res
         ier = 6
         go to 900
      end if

      write(*,'(a)') '  wrote Tecplot PLT: ' // trim(filename)

 900  continue
      deallocate( px, py, pz, nodmap, vrho, vu, vv, vw, vp, vT )
   end subroutine tecplot_write

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

end module mod_uns_output
