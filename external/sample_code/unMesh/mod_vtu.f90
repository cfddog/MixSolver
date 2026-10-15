module vtu_writer
    implicit none
    private

    ! 实数精度 —— 必须 public: 调用方需要用 real(dp) 声明数据
    integer, parameter, public :: dp = kind(1.d0)

    ! -------- VTK 单元类型码 --------
    integer, parameter, public :: VTK_VERTEX     = 1
    integer, parameter, public :: VTK_LINE       = 3
    integer, parameter, public :: VTK_TRIANGLE   = 5
    integer, parameter, public :: VTK_QUAD       = 9
    integer, parameter, public :: VTK_TETRA      = 10
    integer, parameter, public :: VTK_HEXAHEDRON = 12
    integer, parameter, public :: VTK_WEDGE      = 13
    integer, parameter, public :: VTK_PYRAMID    = 14

    ! -------- 通用数据容器 --------
    !   name  : 变量名 (<=64 字符)
    !   ncomp : 分量数 (标量=1, 矢量=3, 张量=9 ...)
    !   val   : (n, ncomp), n = 节点数或单元数
    type, public :: vtu_data
        character(len=64)    :: name  = ''
        integer              :: ncomp = 1
        real(dp), allocatable :: val(:,:)
    end type

    public :: vtu_write

contains

! ==================================================================
!  主入口
!     coords    : (3, npts) 坐标
!     conn      : 扁平连接表, 1-based 节点编号
!     offset    : (ncells), 第 i 个单元结束时的累计节点数
!                 (即 sum(nnode(1:i)))
!     cell_type : (ncells), VTK 类型码
!     point_data, cell_data : 可选数据数组
!     binary    : 可选, 默认 .false. (ASCII)
! ==================================================================
    subroutine vtu_write(fname, coords, conn, offset, cell_type, &
                         point_data, cell_data, binary)
        character(len=*),                   intent(in) :: fname
        real(dp),                           intent(in) :: coords(:,:)
        integer,                            intent(in) :: conn(:)
        integer,                            intent(in) :: offset(:)
        integer,                            intent(in) :: cell_type(:)
        type(vtu_data), optional,           intent(in) :: point_data(:)
        type(vtu_data), optional,           intent(in) :: cell_data(:)
        logical, optional,                  intent(in) :: binary

        logical :: bin

        bin = .false.
        if (present(binary)) bin = binary

        ! ---- 输入检查 ----
        if (size(coords, 1) /= 3) then
            print *, 'vtu_write: ERROR - coords 第一维必须是 3 (x,y,z)'
            return
        end if

        if (size(cell_type) /= size(offset)) then
            print *, 'vtu_write: ERROR - size(offset) /= size(cell_type)'
            return
        end if

        ! 最后一个 offset = 连接表总长度 (累积节点数)
        if (size(offset) > 0) then
            if (offset(size(offset)) /= size(conn)) then
                print *, 'vtu_write: ERROR - offset(末) /= size(conn):', &
                         offset(size(offset)), size(conn)
                return
            end if
        end if

        if (present(point_data)) then
            if (.not. data_ok('PointData', point_data, size(coords, 2))) return
        end if
        if (present(cell_data)) then
            if (.not. data_ok('CellData', cell_data, size(offset))) return
        end if

        if (bin) then
            call write_binary(fname, coords, conn, offset, cell_type, &
                              point_data, cell_data)
        else
            call write_ascii(fname, coords, conn, offset, cell_type, &
                             point_data, cell_data)
        end if
    end subroutine vtu_write

! ==================================================================
!  ASCII 输出
! ==================================================================
    subroutine write_ascii(fname, coords, conn, offset, cell_type, &
                           point_data, cell_data)
        character(len=*),          intent(in) :: fname
        real(dp),                  intent(in) :: coords(:,:)
        integer,                   intent(in) :: conn(:), offset(:), cell_type(:)
        type(vtu_data), optional,  intent(in) :: point_data(:), cell_data(:)

        integer :: u, i, j, k, npts, nconn, ncells, npd, ncd, nc
        character(len=32)  :: str

        npts   = size(coords, 2)
        nconn  = size(conn)
        ncells = size(offset)
        npd    = 0; if (present(point_data)) npd = size(point_data)
        ncd    = 0; if (present(cell_data))  ncd = size(cell_data)

        open(newunit=u, file=trim(fname), status='replace', action='write')

        write(u,'(a)') '<?xml version="1.0"?>'
        write(u,'(a)') '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">'
        write(u,'(a)') '  <UnstructuredGrid>'
        write(u,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', npts, &
                                  '" NumberOfCells="', ncells, '">'

        ! ---- Points ----
        write(u,'(a)') '      <Points>'
        write(u,'(a)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
        do i = 1, npts
            write(u,'(3(1x,es24.16e3))') coords(1,i), coords(2,i), coords(3,i)
        end do
        write(u,'(a)') '        </DataArray>'
        write(u,'(a)') '      </Points>'

        ! ---- Cells ----
        write(u,'(a)') '      <Cells>'

        write(u,'(a)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
        do i = 1, nconn
            write(u,'(i0)') conn(i) - 1
        end do
        write(u,'(a)') '        </DataArray>'

        write(u,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
        do i = 1, ncells
            write(u,'(i0)') offset(i)
        end do
        write(u,'(a)') '        </DataArray>'

        write(u,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
        do i = 1, ncells
            write(u,'(i0)') cell_type(i)
        end do
        write(u,'(a)') '        </DataArray>'

        write(u,'(a)') '      </Cells>'

        ! ---- PointData ----
        if (npd > 0) then
            write(u,'(a)') '      <PointData Scalars="' // trim(point_data(1)%name) // '">'
            do j = 1, npd
                nc = point_data(j)%ncomp
                write(str,'(i0)') nc
                write(u,'(a)') '        <DataArray type="Float64" Name="' // &
                    trim(point_data(j)%name) // '" NumberOfComponents="' // &
                    trim(str) // '" format="ascii">'
                do i = 1, npts
                    write(u,'(*(1x,es24.16e3))') (point_data(j)%val(i,k), k=1,nc)
                end do
                write(u,'(a)') '        </DataArray>'
            end do
            write(u,'(a)') '      </PointData>'
        end if

        ! ---- CellData ----
        if (ncd > 0) then
            write(u,'(a)') '      <CellData Scalars="' // trim(cell_data(1)%name) // '">'
            do j = 1, ncd
                nc = cell_data(j)%ncomp
                write(str,'(i0)') nc
                write(u,'(a)') '        <DataArray type="Float64" Name="' // &
                    trim(cell_data(j)%name) // '" NumberOfComponents="' // &
                    trim(str) // '" format="ascii">'
                do i = 1, ncells
                    write(u,'(*(1x,es24.16e3))') (cell_data(j)%val(i,k), k=1,nc)
                end do
                write(u,'(a)') '        </DataArray>'
            end do
            write(u,'(a)') '      </CellData>'
        end if

        write(u,'(a)') '    </Piece>'
        write(u,'(a)') '  </UnstructuredGrid>'
        write(u,'(a)') '</VTKFile>'

        close(u)
    end subroutine write_ascii

! ==================================================================
!  Binary (appended raw) 输出
! ==================================================================
    subroutine write_binary(fname, coords, conn, offset, cell_type, &
                            point_data, cell_data)
        character(len=*),          intent(in) :: fname
        real(dp),                  intent(in) :: coords(:,:)
        integer,                   intent(in) :: conn(:), offset(:), cell_type(:)
        type(vtu_data), optional,  intent(in) :: point_data(:), cell_data(:)

        integer :: u, j, npts, nconn, ncells, npd, ncd, nc
        integer(4) :: bsize                  ! 4 字节 (UInt32) 块长度头
        integer(8) :: off_pts, off_conn, off_offs, off_typs, off_cur
        integer(8), allocatable :: off_pdata(:), off_cdata(:)

        real(4),    allocatable :: coords32(:,:)
        integer(4), allocatable :: conn32(:), offs32(:)
        integer(1), allocatable :: typs8(:)
        real(4),    allocatable :: tmp32(:,:)

        character(len=512) :: line
        character(len=32)  :: str
        character(len=32)  :: strc

        npts   = size(coords, 2)
        nconn  = size(conn)
        ncells = size(offset)
        npd    = 0; if (present(point_data)) npd = size(point_data)
        ncd    = 0; if (present(cell_data))  ncd = size(cell_data)

        ! ---- 打包数据 (转 Float32 / Int32 / UInt8) ----
        allocate(coords32(3, npts));  coords32 = real(coords, 4)
        allocate(conn32(nconn));      conn32   = int(conn - 1, 4)
        allocate(offs32(ncells));     offs32   = int(offset, 4)
        allocate(typs8(ncells));      typs8    = int(cell_type, 1)

        ! ---- 计算 appended 偏移量 (相对 '_' 之后第一个字节) ----
        !      每个块 = 4 字节 UInt32 长度头 + 数据字节; 用 8 字节算术防溢出
        off_pts  = 0_8
        off_conn = off_pts  + 4_8 + int(npts,8)*3_8*4_8
        off_offs = off_conn + 4_8 + int(nconn,8)*4_8
        off_typs = off_offs + 4_8 + int(ncells,8)*4_8

        allocate(off_pdata(max(npd,1)), off_cdata(max(ncd,1)))
        off_pdata = 0_8; off_cdata = 0_8

        ! 第一个点/单元数据块: 紧跟 types 块; 若无点数据, 单元数据直接从这里开始
        off_cur = off_typs + 4_8 + int(ncells,8)

        do j = 1, npd
            off_pdata(j) = off_cur
            off_cur = off_cur + 4_8 + int(npts,8)*int(point_data(j)%ncomp,8)*4_8
        end do

        do j = 1, ncd
            off_cdata(j) = off_cur
            off_cur = off_cur + 4_8 + int(ncells,8)*int(cell_data(j)%ncomp,8)*4_8
        end do

        ! ---- 以 stream 打开 ----
        open(newunit=u, file=trim(fname), status='replace', action='write', &
             form='unformatted', access='stream')

        ! ---- XML 头 ----
        write(u) '<?xml version="1.0"?>' // char(10)
        write(u) '<VTKFile type="UnstructuredGrid" version="1.0" ' // &
                 'byte_order="LittleEndian" header_type="UInt32">' // char(10)
        write(u) '  <UnstructuredGrid>' // char(10)
        write(line,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', npts, &
                                     '" NumberOfCells="',  ncells, '">'
        write(u) trim(line) // char(10)

        write(u) '      <Points>' // char(10)
        write(str,'(i0)') off_pts
        write(u) '        <DataArray type="Float32" NumberOfComponents="3" ' // &
                 'format="appended" offset="' // trim(str) // '"/>' // char(10)
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

        if (npd > 0) then
            write(u) '      <PointData Scalars="' // trim(point_data(1)%name) // '">' // char(10)
            do j = 1, npd
                nc = point_data(j)%ncomp
                write(strc,'(i0)') nc
                write(str ,'(i0)') off_pdata(j)
                write(u) '        <DataArray type="Float32" Name="' // &
                    trim(point_data(j)%name) // '" NumberOfComponents="' // &
                    trim(strc) // '" format="appended" offset="' // &
                    trim(str) // '"/>' // char(10)
            end do
            write(u) '      </PointData>' // char(10)
        end if

        if (ncd > 0) then
            write(u) '      <CellData Scalars="' // trim(cell_data(1)%name) // '">' // char(10)
            do j = 1, ncd
                nc = cell_data(j)%ncomp
                write(strc,'(i0)') nc
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

        ! ---- Appended 数据 ----
        write(u) '  <AppendedData encoding="raw">' // char(10) // '_'

        bsize = u32_head(int(npts,8)*12_8);      write(u) bsize;  write(u) coords32
        bsize = u32_head(int(nconn,8)*4_8);      write(u) bsize;  write(u) conn32
        bsize = u32_head(int(ncells,8)*4_8);     write(u) bsize;  write(u) offs32
        bsize = u32_head(int(ncells,8));         write(u) bsize;  write(u) typs8

        do j = 1, npd
            nc = point_data(j)%ncomp
            allocate(tmp32(nc, npts))
            tmp32 = real(transpose(point_data(j)%val), 4)   ! 交错分量
            bsize = u32_head(int(npts,8)*int(nc,8)*4_8)
            write(u) bsize;  write(u) tmp32
            deallocate(tmp32)
        end do

        do j = 1, ncd
            nc = cell_data(j)%ncomp
            allocate(tmp32(nc, ncells))
            tmp32 = real(transpose(cell_data(j)%val), 4)
            bsize = u32_head(int(ncells,8)*int(nc,8)*4_8)
            write(u) bsize;  write(u) tmp32
            deallocate(tmp32)
        end do

        write(u) char(10) // '  </AppendedData>' // char(10)
        write(u) '</VTKFile>' // char(10)

        close(u)
    end subroutine write_binary

! ==================================================================
!  内部工具
! ==================================================================
!  由 8 字节字节数得到 4 字节 (UInt32) 长度头
!  (单块 > 4 GB 时 32 位头表示不了, 报错终止而不是写坏文件)
    integer(4) function u32_head(nbytes) result(h)
        integer(8), intent(in) :: nbytes

        if (nbytes < 0_8 .or. nbytes > 4294967295_8) then
            print *, 'vtu_write: ERROR - 数据块字节数超出 UInt32 头范围:', nbytes
            error stop
        end if
        h = int(nbytes, 4)
    end function u32_head

! ------------------------------------------------------------------
!  检查数据数组: 必须已分配, 且形状为 (n, ncomp)
! ------------------------------------------------------------------
    logical function data_ok(label, data, n) result(ok)
        character(len=*), intent(in) :: label
        type(vtu_data),   intent(in) :: data(:)
        integer,          intent(in) :: n

        integer :: j
        character(len=16) :: s

        ok = .true.
        do j = 1, size(data)
            write(s,'(i0)') j
            if (.not. allocated(data(j)%val)) then
                print *, 'vtu_write: ERROR - '//trim(label)//'('//trim(s)//')%val 未分配'
                ok = .false.
                return
            end if
            if (data(j)%ncomp < 1) then
                print *, 'vtu_write: ERROR - '//trim(label)//'('//trim(s)//')%ncomp < 1'
                ok = .false.
                return
            end if
            if (size(data(j)%val, 1) /= n) then
                print *, 'vtu_write: ERROR - '//trim(label)//'('//trim(s)//')%val 行数', &
                         size(data(j)%val, 1), ' /= ', n
                ok = .false.
                return
            end if
            if (size(data(j)%val, 2) /= data(j)%ncomp) then
                print *, 'vtu_write: ERROR - '//trim(label)//'('//trim(s)//')%val 列数', &
                         size(data(j)%val, 2), ' /= ncomp ', data(j)%ncomp
                ok = .false.
                return
            end if
        end do
    end function data_ok

end module vtu_writer