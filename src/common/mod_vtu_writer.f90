!===============================================================================
! mod_vtu_writer.f90 -- standards-conformant VTK XML (.vtu) writer
!
! Ported from the reference sample external/sample_code/unMesh/mod_vtu.f90.
! The XML layout produced here is what the VTK XML specification prescribes and
! what that sample writes (the sample was verified to open correctly in both
! ParaView and Tecplot 360 2020):
!
!   * <connectivity> holds ONLY the 0-based node ids, back to back;
!   * <offsets> holds the cumulative node count, i.e. the nodes of cell i are
!     connectivity(offsets(i-1)+1 : offsets(i)) and the node count of cell i is
!     offset(i) - offset(i-1), which MUST equal the node count of the cell's VTK
!     type (4 for VTK_TETRA, 8 for VTK_HEXAHEDRON, ...);
!   * <types> holds the VTK cell type code of every cell (UInt8).
!
! The <PointData> / <CellData> tags are written WITHOUT attributes (no
! Scalars="..." hint).  That is legal XML/VTK, and it is required for
! compatibility with this repository: the validation scripts under
! cases/**/*.py locate the data blocks with a strict
! r'<CellData>(.*?)</CellData>' regular expression (historical MixNSSolver
! output had bare tags), so adding an attribute silently breaks every one of
! them -- exactly what happened for one build on 2026-10-15.
!
! NEVER write the legacy VTK "CELLS" convention (a leading vertex-count token
! stored inside connectivity): VTK-based readers tolerate it only by accident
! (they recompute the cell location as offsets - npts(type)), while a stricter
! reader such as Tecplot 360 reads 5 node ids for a tetrahedron and renders a
! broken mesh.  That exact defect was fixed in MixNSSolver on 2026-10-15.
!
! Conventions
!   coords(:,i)  = coordinates (x,y,z) of point i (first dimension must be 3)
!   conn(:)      = flat 1-based connectivity: all nodes of cell 1, then cell 2..
!   offset(i)    = cumulative node count through cell i = sum(nnode(1:i))
!   cell_type(i) = VTK_VERTEX / VTK_TETRA / VTK_HEXAHEDRON / ... (codes below)
!   point_data / cell_data : val(:,:) has shape (n, ncomp); the components of
!     one point/cell are stored contiguously, points/cells back to back.
!
! Both an inline ASCII file (default, full double precision) and an appended
! binary file (Float32 payload) can be written; the ASCII path is the one used
! by the solver.
!===============================================================================
module mod_vtu_writer
   use mod_precision, only: dp, sp, i8
   implicit none
   private

   public :: vtu_data, vtu_write

   ! ---- integer kind used by the appended-binary path ------------------------
   !      VTK header_type="UInt32" => 4-byte block sizes, Int32 / UInt8 data
   integer, parameter :: i1 = selected_int_kind(2)
   integer, parameter :: i4 = selected_int_kind(9)

   ! ---- VTK cell type codes (VTK file format specification) -----------------
   integer, parameter, public :: VTK_VERTEX     = 1
   integer, parameter, public :: VTK_LINE       = 3
   integer, parameter, public :: VTK_TRIANGLE   = 5
   integer, parameter, public :: VTK_QUAD       = 9
   integer, parameter, public :: VTK_TETRA      = 10
   integer, parameter, public :: VTK_HEXAHEDRON = 12
   integer, parameter, public :: VTK_WEDGE      = 13
   integer, parameter, public :: VTK_PYRAMID    = 14

   ! ---- generic data container ---------------------------------------------
   !   name  : variable name (up to 64 characters)
   !   ncomp : number of components (scalar = 1, vector = 3, tensor = 9, ...)
   !   val   : (n, ncomp) array, n = number of points (PointData) or cells
   type, public :: vtu_data
      character(len=64)     :: name  = ''
      integer               :: ncomp = 1
      real(dp), allocatable :: val(:,:)
   end type vtu_data

contains

   !============================================================================
   ! Main entry point.
   !   ier = 0 on success; ier /= 0 (plus an English message on stdout) when an
   !   argument is inconsistent or the file cannot be opened, so that callers
   !   can propagate the failure through their own ier.
   !============================================================================
   subroutine vtu_write( fname, coords, conn, offset, cell_type, &
                         point_data, cell_data, binary, ier )
      character(len=*),          intent(in)  :: fname
      real(dp),                  intent(in)  :: coords(:,:)
      integer,                   intent(in)  :: conn(:)
      integer,                   intent(in)  :: offset(:)
      integer,                   intent(in)  :: cell_type(:)
      type(vtu_data), optional,  intent(in)  :: point_data(:)
      type(vtu_data), optional,  intent(in)  :: cell_data(:)
      logical,        optional,  intent(in)  :: binary
      integer,        optional,  intent(out) :: ier

      logical :: bin

      call set_ier( ier, 0 )
      bin = .false.
      if ( present(binary) ) bin = binary

      ! ---- argument consistency checks ----
      if ( size(coords, 1) /= 3 ) then
         write(*,'(a)') 'vtu_write: ERROR - coords first dimension must be 3 (x,y,z)'
         call set_ier( ier, 1 )
         return
      end if
      if ( size(cell_type) /= size(offset) ) then
         write(*,'(a)') 'vtu_write: ERROR - size(offset) /= size(cell_type)'
         call set_ier( ier, 1 )
         return
      end if
      if ( size(offset) > 0 ) then
         if ( offset(size(offset)) /= size(conn) ) then
            write(*,'(a)') 'vtu_write: ERROR - offset(last) /= size(conn)'
            call set_ier( ier, 1 )
            return
         end if
      end if
      if ( present(point_data) ) then
         if ( .not. data_ok( 'PointData', point_data, size(coords,2) ) ) then
            call set_ier( ier, 1 )
            return
         end if
      end if
      if ( present(cell_data) ) then
         if ( .not. data_ok( 'CellData', cell_data, size(offset) ) ) then
            call set_ier( ier, 1 )
            return
         end if
      end if

      if ( bin ) then
         call write_binary( fname, coords, conn, offset, cell_type, &
                            point_data, cell_data, ier )
      else
         call write_ascii( fname, coords, conn, offset, cell_type, &
                           point_data, cell_data, ier )
      end if
   end subroutine vtu_write

   !============================================================================
   ! Inline ASCII output (default).  es24.16e3 keeps all double-precision
   ! digits; one integer per line for connectivity (as in the reference
   ! sample, which the strict Tecplot reader accepts).
   !============================================================================
   subroutine write_ascii( fname, coords, conn, offset, cell_type, &
                           point_data, cell_data, ier )
      character(len=*),         intent(in)  :: fname
      real(dp),                 intent(in)  :: coords(:,:)
      integer,                  intent(in)  :: conn(:), offset(:), cell_type(:)
      type(vtu_data), optional, intent(in)  :: point_data(:), cell_data(:)
      integer,        optional, intent(out) :: ier

      integer :: u, ios, i, j, npts, nconn, ncells, npd, ncd, nc, np
      character(len=64) :: str

      npts   = size(coords, 2)
      nconn  = size(conn)
      ncells = size(offset)
      npd    = 0
      if ( present(point_data) ) npd = size(point_data)
      ncd    = 0
      if ( present(cell_data) )  ncd = size(cell_data)

      open( newunit=u, file=trim(fname), status='replace', action='write', &
            iostat=ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'vtu_write: ERROR - cannot open ' // trim(fname)
         call set_ier( ier, 1 )
         return
      end if

      ! ---- XML header ----
      !      Name="Coordinates" on the points array is kept for backward
      !      compatibility with the existing post-processing scripts.
      write(u,'(a)') '<?xml version="1.0"?>'
      write(u,'(a)') '<VTKFile type="UnstructuredGrid" version="1.0" ' // &
                     'byte_order="LittleEndian">'
      write(u,'(a)') '  <UnstructuredGrid>'
      write(u,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', npts, &
                               '" NumberOfCells="', ncells, '">'

      ! ---- points ----
      write(u,'(a)') '      <Points>'
      write(u,'(a)') '        <DataArray type="Float64" Name="Coordinates" ' // &
                     'NumberOfComponents="3" format="ascii">'
      do i = 1, npts
         write(u,'(3(1x,es24.16e3))') coords(1,i), coords(2,i), coords(3,i)
      end do
      write(u,'(a)') '        </DataArray>'
      write(u,'(a)') '      </Points>'

      ! ---- cells: pure node-id connectivity + cumulative offsets + types ----
      write(u,'(a)') '      <Cells>'
      write(u,'(a)') '        <DataArray type="Int32" Name="connectivity" ' // &
                     'format="ascii">'
      do i = 1, nconn
         write(u,'(i0)') conn(i) - 1          ! VTK XML uses 0-based ids
      end do
      write(u,'(a)') '        </DataArray>'
      write(u,'(a)') '        <DataArray type="Int32" Name="offsets" ' // &
                     'format="ascii">'
      do i = 1, ncells
         write(u,'(i0)') offset(i)
      end do
      write(u,'(a)') '        </DataArray>'
      write(u,'(a)') '        <DataArray type="UInt8" Name="types" ' // &
                     'format="ascii">'
      do i = 1, ncells
         write(u,'(i0)') cell_type(i)
      end do
      write(u,'(a)') '        </DataArray>'
      write(u,'(a)') '      </Cells>'

      ! ---- point data ----
      ! NB: the tag is written WITHOUT attributes on purpose: the case scripts
      ! under cases/**/*.py locate the data blocks with a strict
      ! r'<CellData>(.*?)</CellData>' style regex (the historical MixNSSolver
      ! output had bare tags), so a Scalars= hint would break every one of
      ! them.  The attribute is optional in the VTK XML specification.
      if ( npd > 0 ) then
         write(u,'(a)') '      <PointData>'
         do j = 1, npd
            np = point_data(j)%ncomp
            write(str,'(i0)') np
            write(u,'(a)') '        <DataArray type="Float64" Name="' // &
                           trim(point_data(j)%name) // &
                           '" NumberOfComponents="' // trim(str) // &
                           '" format="ascii">'
            do i = 1, npts
               write(u,'(*(1x,es24.16e3))') point_data(j)%val(i,1:np)
            end do
            write(u,'(a)') '        </DataArray>'
         end do
         write(u,'(a)') '      </PointData>'
      end if

      ! ---- cell data ---- (bare tag, see the PointData note above)
      if ( ncd > 0 ) then
         write(u,'(a)') '      <CellData>'
         do j = 1, ncd
            nc = cell_data(j)%ncomp
            write(str,'(i0)') nc
            write(u,'(a)') '        <DataArray type="Float64" Name="' // &
                           trim(cell_data(j)%name) // &
                           '" NumberOfComponents="' // trim(str) // &
                           '" format="ascii">'
            do i = 1, ncells
               write(u,'(*(1x,es24.16e3))') cell_data(j)%val(i,1:nc)
            end do
            write(u,'(a)') '        </DataArray>'
         end do
         write(u,'(a)') '      </CellData>'
      end if

      write(u,'(a)') '    </Piece>'
      write(u,'(a)') '  </UnstructuredGrid>'
      write(u,'(a)') '</VTKFile>'
      close(u)

      call set_ier( ier, 0 )
   end subroutine write_ascii

   !============================================================================
   ! Appended raw binary output (Float32 payload).  Not used by the solver but
   ! kept from the reference sample so the module stays a drop-in writer.
   !============================================================================
   subroutine write_binary( fname, coords, conn, offset, cell_type, &
                            point_data, cell_data, ier )
      character(len=*),         intent(in)  :: fname
      real(dp),                 intent(in)  :: coords(:,:)
      integer,                  intent(in)  :: conn(:), offset(:), cell_type(:)
      type(vtu_data), optional, intent(in)  :: point_data(:), cell_data(:)
      integer,        optional, intent(out) :: ier

      integer :: u, ios, j, npts, nconn, ncells, npd, ncd, nc
      integer(i4) :: bsize                      ! 4-byte (UInt32) block length
      integer(i8) :: off_pts, off_conn, off_offs, off_typs, off_cur
      integer(i8), allocatable :: off_pdata(:), off_cdata(:)

      real(sp),    allocatable :: coords32(:,:), tmp32(:,:)
      integer(i4), allocatable :: conn32(:), offs32(:)
      integer(i1), allocatable :: typs8(:)

      character(len=64) :: str, strc

      npts   = size(coords, 2)
      nconn  = size(conn)
      ncells = size(offset)
      npd    = 0
      if ( present(point_data) ) npd = size(point_data)
      ncd    = 0
      if ( present(cell_data) )  ncd = size(cell_data)

      ! ---- pack the payload (Float32 / Int32 / UInt8) ----
      allocate( coords32(3, npts) )
      allocate( conn32(nconn), offs32(ncells), typs8(ncells) )
      coords32 = real(coords, sp)
      conn32   = int(conn - 1, i4)              ! VTK XML uses 0-based ids
      offs32   = int(offset, i4)
      typs8    = int(cell_type, i1)

      ! ---- appended offsets (relative to the first byte after "_") ----
      !      each block = 4-byte UInt32 length header + payload bytes; the
      !      8-byte arithmetic below cannot overflow for realistic meshes
      off_pts  = 0_i8
      off_conn = off_pts  + 4_i8 + int(npts,   i8) * 3_i8 * 4_i8
      off_offs = off_conn + 4_i8 + int(nconn,  i8) * 4_i8
      off_typs = off_offs + 4_i8 + int(ncells, i8) * 4_i8
      allocate( off_pdata(max(npd,1)), off_cdata(max(ncd,1)) )
      off_pdata = 0_i8
      off_cdata = 0_i8
      off_cur = off_typs + 4_i8 + int(ncells, i8)
      do j = 1, npd
         off_pdata(j) = off_cur
         off_cur = off_cur + 4_i8 + int(npts,i8) * int(point_data(j)%ncomp,i8) * 4_i8
      end do
      do j = 1, ncd
         off_cdata(j) = off_cur
         off_cur = off_cur + 4_i8 + int(ncells,i8) * int(cell_data(j)%ncomp,i8) * 4_i8
      end do

      ! ---- open as a byte stream ----
      open( newunit=u, file=trim(fname), status='replace', action='write', &
            form='unformatted', access='stream', iostat=ios )
      if ( ios /= 0 ) then
         write(*,'(a)') 'vtu_write: ERROR - cannot open ' // trim(fname)
         call set_ier( ier, 1 )
         return
      end if

      ! ---- XML header ----
      write(u) '<?xml version="1.0"?>' // char(10)
      write(u) '<VTKFile type="UnstructuredGrid" version="1.0" ' // &
               'byte_order="LittleEndian" header_type="UInt32">' // char(10)
      write(u) '  <UnstructuredGrid>' // char(10)
      write(u,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', npts, &
                               '" NumberOfCells="', ncells, '">' // char(10)

      write(u) '      <Points>' // char(10)
      write(str,'(i0)') off_pts
      write(u) '        <DataArray type="Float32" Name="Coordinates" ' // &
               'NumberOfComponents="3" format="appended" offset="' // &
               trim(str) // '"/>' // char(10)
      write(u) '      </Points>' // char(10)

      write(u) '      <Cells>' // char(10)
      write(str,'(i0)') off_conn
      write(u) '        <DataArray type="Int32" Name="connectivity" ' // &
               'format="appended" offset="' // trim(str) // '"/>' // char(10)
      write(str,'(i0)') off_offs
      write(u) '        <DataArray type="Int32" Name="offsets" ' // &
               'format="appended" offset="' // trim(str) // '"/>' // char(10)
      write(str,'(i0)') off_typs
      write(u) '        <DataArray type="UInt8" Name="types" ' // &
               'format="appended" offset="' // trim(str) // '"/>' // char(10)
      write(u) '      </Cells>' // char(10)

      if ( npd > 0 ) then
         write(u) '      <PointData>' // char(10)
         do j = 1, npd
            write(strc,'(i0)') point_data(j)%ncomp
            write(str ,'(i0)') off_pdata(j)
            write(u) '        <DataArray type="Float32" Name="' // &
                     trim(point_data(j)%name) // '" NumberOfComponents="' // &
                     trim(strc) // '" format="appended" offset="' // &
                     trim(str) // '"/>' // char(10)
         end do
         write(u) '      </PointData>' // char(10)
      end if

      if ( ncd > 0 ) then
         write(u) '      <CellData>' // char(10)
         do j = 1, ncd
            write(strc,'(i0)') cell_data(j)%ncomp
            write(str ,'(i0)') off_cdata(j)
            write(u) '        <DataArray type="Float32" Name="' // &
                     trim(cell_data(j)%name) // '" NumberOfComponents="' // &
                     trim(strc) // '" format="appended" offset="' // &
                     trim(str) // '"/>' // char(10)
         end do
         write(u) '      </CellData>' // char(10)
      end if

      write(u) '    </Piece>' // char(10)
      write(u) '  </UnstructuredGrid>' // char(10)

      ! ---- appended payload ----
      write(u) '  <AppendedData encoding="raw">' // char(10) // '_'

      bsize = u32_head( int(npts,   i8) * 12_i8 );  write(u) bsize; write(u) coords32
      bsize = u32_head( int(nconn,  i8) *  4_i8 );  write(u) bsize; write(u) conn32
      bsize = u32_head( int(ncells, i8) *  4_i8 );  write(u) bsize; write(u) offs32
      bsize = u32_head( int(ncells, i8) *  1_i8 );  write(u) bsize; write(u) typs8

      do j = 1, npd
         nc = point_data(j)%ncomp
         allocate( tmp32(nc, npts) )
         tmp32 = real(transpose(point_data(j)%val), sp)  ! interleave components
         bsize = u32_head( int(npts, i8) * int(nc, i8) * 4_i8 )
         write(u) bsize; write(u) tmp32
         deallocate(tmp32)
      end do

      do j = 1, ncd
         nc = cell_data(j)%ncomp
         allocate( tmp32(nc, ncells) )
         tmp32 = real(transpose(cell_data(j)%val), sp)
         bsize = u32_head( int(ncells, i8) * int(nc, i8) * 4_i8 )
         write(u) bsize; write(u) tmp32
         deallocate(tmp32)
      end do

      write(u) char(10) // '  </AppendedData>' // char(10)
      write(u) '</VTKFile>' // char(10)

      close(u)

      call set_ier( ier, 0 )
   end subroutine write_binary

   !============================================================================
   ! Internal helpers
   !============================================================================
   ! Convert a byte count into the 4-byte (UInt32) block length header.  A block
   ! larger than 4 GB cannot be represented by a 32-bit header, so stop instead
   ! of writing a corrupt file.
   integer(i4) function u32_head( nbytes ) result( h )
      integer(i8), intent(in) :: nbytes

      if ( nbytes < 0_i8 .or. nbytes > 4294967295_i8 ) then
         write(*,'(a,i0)') 'vtu_write: ERROR - block size exceeds the ' // &
                           'UInt32 header range: ', nbytes
         error stop
      end if
      h = int(nbytes, i4)
   end function u32_head

   ! Check one data array: it must be allocated with shape (n, ncomp).
   logical function data_ok( label, data, n ) result( ok )
      character(len=*), intent(in) :: label
      type(vtu_data),   intent(in) :: data(:)
      integer,          intent(in) :: n

      integer :: j
      character(len=16) :: s

      ok = .true.
      do j = 1, size(data)
         write(s,'(i0)') j
         if ( .not. allocated(data(j)%val) ) then
            write(*,'(a)') 'vtu_write: ERROR - ' // trim(label) // '(' // &
                           trim(s) // ')%val is not allocated'
            ok = .false.
            return
         end if
         if ( data(j)%ncomp < 1 ) then
            write(*,'(a)') 'vtu_write: ERROR - ' // trim(label) // '(' // &
                           trim(s) // ')%ncomp < 1'
            ok = .false.
            return
         end if
         if ( size(data(j)%val, 1) /= n ) then
            write(*,'(a,i0,a,i0)') 'vtu_write: ERROR - ' // trim(label) // '(' // &
                 trim(s) // ')%val rows = ', size(data(j)%val, 1), ' /= ', n
            ok = .false.
            return
         end if
         if ( size(data(j)%val, 2) /= data(j)%ncomp ) then
            write(*,'(a,i0,a,i0)') 'vtu_write: ERROR - ' // trim(label) // '(' // &
                 trim(s) // ')%val cols = ', size(data(j)%val, 2), &
                 ' /= ncomp ', data(j)%ncomp
            ok = .false.
            return
         end if
      end do
   end function data_ok

   ! Store the error code when the caller asked for one (no-op otherwise).
   subroutine set_ier( ier, code )
      integer, optional, intent(out) :: ier
      integer,          intent(in)  :: code

      if ( present(ier) ) ier = code
   end subroutine set_ier

end module mod_vtu_writer
