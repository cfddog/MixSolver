! ==================================================================
!  test_vtu : vtu_writer 模块的使用示例 / 自检程序
!
!  编译:  gfortran -O2 mod_vtu.f90 main.f90 -o test   (或者 make)
!  运行:  ./test
!  输出:  out_ascii.vtu   / out_binary.vtu    结构化六面体网格
!         mixed_ascii.vtu / mixed_binary.vtu  混合单元类型网格
!
!  约定 (详见 mod_vtu.f90)
!     coords(:,i)  = 第 i 个节点坐标 (x, y, z)
!     conn(:)      = 扁平 1-based 连接表: 第 1 个单元的全部节点, 然后第 2 个 ...
!     offset(i)    = 到第 i 个单元为止的累计节点数 (= sum(nnode(1:i)))
!     cell_type(i) = VTK_TETRA / VTK_HEXAHEDRON / ...
! ==================================================================
program test
    use vtu_writer
    implicit none

    call run_structured_hex()   ! 算例 1: 结构化六面体网格
    call run_mixed_cells()      ! 算例 2: 混合单元类型 (hex/wedge/tet/pyramid)

contains

! ------------------------------------------------------------------
!  算例 1: nx*ny*nz 个 VTK_HEXAHEDRON 组成的单位立方体网格
! ------------------------------------------------------------------
    subroutine run_structured_hex()
        integer, parameter :: nx = 5, ny = 4, nz = 3
        integer, parameter :: npts   = (nx+1)*(ny+1)*(nz+1)
        integer, parameter :: ncells = nx*ny*nz

        real(dp) :: coords(3, npts)
        integer  :: conn(8*ncells), offset(ncells), cell_type(ncells)
        type(vtu_data) :: pdata(2), cdata(2)

        integer  :: i, j, k, n, c
        real(dp) :: xc, yc, zc

        ! ---- 1. 节点坐标 (单位立方体均匀网格) ----
        n = 0
        do k = 0, nz
            do j = 0, ny
                do i = 0, nx
                    n = n + 1
                    coords(1,n) = dble(i)/dble(nx)
                    coords(2,n) = dble(j)/dble(ny)
                    coords(3,n) = dble(k)/dble(nz)
                end do
            end do
        end do

        ! ---- 2. 连接表 (VTK 六面体节点顺序: 底面 4 点, 然后顶面 4 点) ----
        c = 0
        do k = 0, nz-1
            do j = 0, ny-1
                do i = 0, nx-1
                    c = c + 1
                    conn(8*c-7) = hexid(i  , j  , k  , nx, ny)
                    conn(8*c-6) = hexid(i+1, j  , k  , nx, ny)
                    conn(8*c-5) = hexid(i+1, j+1, k  , nx, ny)
                    conn(8*c-4) = hexid(i  , j+1, k  , nx, ny)
                    conn(8*c-3) = hexid(i  , j  , k+1, nx, ny)
                    conn(8*c-2) = hexid(i+1, j  , k+1, nx, ny)
                    conn(8*c-1) = hexid(i+1, j+1, k+1, nx, ny)
                    conn(8*c  ) = hexid(i  , j+1, k+1, nx, ny)
                    offset(c)    = 8*c
                    cell_type(c) = VTK_HEXAHEDRON
                end do
            end do
        end do

        ! ---- 3. 节点数据: 标量 rho + 矢量 velocity ----
        pdata(1)%name  = 'rho'
        pdata(1)%ncomp = 1
        allocate(pdata(1)%val(npts, 1))
        pdata(1)%val(:,1) = 1.0_dp + coords(1,:) + 2.0_dp*coords(2,:) + 3.0_dp*coords(3,:)

        pdata(2)%name  = 'velocity'
        pdata(2)%ncomp = 3
        allocate(pdata(2)%val(npts, 3))
        pdata(2)%val(:,1) = coords(2,:) - coords(3,:)   ! u = y - z
        pdata(2)%val(:,2) = coords(3,:) - coords(1,:)   ! v = z - x
        pdata(2)%val(:,3) = coords(1,:) - coords(2,:)   ! w = x - y

        ! ---- 4. 单元数据 (用单元中心): 标量 pressure + 矢量 cell_velocity ----
        cdata(1)%name  = 'pressure'
        cdata(1)%ncomp = 1
        cdata(2)%name  = 'cell_velocity'
        cdata(2)%ncomp = 3
        allocate(cdata(1)%val(ncells, 1), cdata(2)%val(ncells, 3))

        do c = 1, ncells
            xc = 0.0_dp;  yc = 0.0_dp;  zc = 0.0_dp
            do n = 8*(c-1)+1, 8*c
                xc = xc + coords(1, conn(n))
                yc = yc + coords(2, conn(n))
                zc = zc + coords(3, conn(n))
            end do
            xc = xc/8.0_dp;  yc = yc/8.0_dp;  zc = zc/8.0_dp
            cdata(1)%val(c,1) = 1.0_dp + xc + 2.0_dp*yc + 3.0_dp*zc
            cdata(2)%val(c,1) = yc - zc
            cdata(2)%val(c,2) = zc - xc
            cdata(2)%val(c,3) = xc - yc
        end do

        ! ---- 5. 输出 (ASCII 与二进制各一份) ----
        call vtu_write('out_ascii.vtu',  coords, conn, offset, cell_type, &
                       point_data=pdata, cell_data=cdata, binary=.false.)
        call vtu_write('out_binary.vtu', coords, conn, offset, cell_type, &
                       point_data=pdata, cell_data=cdata, binary=.true.)

        print '(a,i0,a,i0,a)', '算例1 结构化六面体: ', npts, ' 节点 / ', ncells, &
                               ' 单元 -> out_ascii.vtu, out_binary.vtu'
    end subroutine run_structured_hex

    ! 结构化网格节点编号 (i,j,k 从 0 开始) -> 1-based 节点号
    integer function hexid(i, j, k, nx, ny)
        integer, intent(in) :: i, j, k, nx, ny
        hexid = 1 + i + (nx+1)*(j + (ny+1)*k)
    end function hexid

! ------------------------------------------------------------------
!  算例 2: 混合单元类型 —— 同一个文件里同时有
!          六面体 / 三棱柱(WEDGE) / 四面体 / 四棱锥(PYRAMID)
!          用来验证变长度连接表的 offset 与类型码
! ------------------------------------------------------------------
    subroutine run_mixed_cells()
        integer, parameter :: npts = 23, ncells = 4

        real(dp) :: coords(3, npts)
        integer  :: conn(npts), offset(ncells), cell_type(ncells)
        type(vtu_data) :: pdata(1), cdata(2)

        integer  :: c, n, n1, n2, ncell
        real(dp) :: xc, yc, zc

        ! ---- 1. 节点坐标 (各单元在空间上分开, 便于目视检查) ----
        ! 单元 1: 六面体 (节点 1-8)
        coords(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        coords(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        coords(:, 3) = [1.0_dp, 1.0_dp, 0.0_dp]
        coords(:, 4) = [0.0_dp, 1.0_dp, 0.0_dp]
        coords(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        coords(:, 6) = [1.0_dp, 0.0_dp, 1.0_dp]
        coords(:, 7) = [1.0_dp, 1.0_dp, 1.0_dp]
        coords(:, 8) = [0.0_dp, 1.0_dp, 1.0_dp]
        ! 单元 2: 三棱柱 (节点 9-14), VTK 顺序 = 下三角 3 点 + 上三角 3 点
        coords(:, 9) = [2.0_dp, 0.0_dp, 0.0_dp]
        coords(:,10) = [3.0_dp, 0.0_dp, 0.0_dp]
        coords(:,11) = [2.5_dp, 1.0_dp, 0.0_dp]
        coords(:,12) = [2.0_dp, 0.0_dp, 1.0_dp]
        coords(:,13) = [3.0_dp, 0.0_dp, 1.0_dp]
        coords(:,14) = [2.5_dp, 1.0_dp, 1.0_dp]
        ! 单元 3: 四面体 (节点 15-18)
        coords(:,15) = [0.0_dp, 0.0_dp, 2.0_dp]
        coords(:,16) = [1.0_dp, 0.0_dp, 2.0_dp]
        coords(:,17) = [0.0_dp, 1.0_dp, 2.0_dp]
        coords(:,18) = [0.0_dp, 0.0_dp, 3.0_dp]
        ! 单元 4: 四棱锥 (节点 19-23), VTK 顺序 = 底面 4 点 + 顶点
        coords(:,19) = [2.0_dp, 0.0_dp, 2.0_dp]
        coords(:,20) = [3.0_dp, 0.0_dp, 2.0_dp]
        coords(:,21) = [3.0_dp, 1.0_dp, 2.0_dp]
        coords(:,22) = [2.0_dp, 1.0_dp, 2.0_dp]
        coords(:,23) = [2.5_dp, 0.5_dp, 3.0_dp]

        ! ---- 2. 连接表: 各单元节点连续编号 -> conn(n) = n, offset = 累计节点数 ----
        do n = 1, npts
            conn(n) = n
        end do
        offset    = [8, 14, 18, 23]
        cell_type = [VTK_HEXAHEDRON, VTK_WEDGE, VTK_TETRA, VTK_PYRAMID]

        ! ---- 3. 节点数据 ----
        pdata(1)%name  = 'z'
        pdata(1)%ncomp = 1
        allocate(pdata(1)%val(npts, 1))
        do n = 1, npts
            pdata(1)%val(n,1) = coords(3,n)
        end do

        ! ---- 4. 单元数据: 单元号 (标量) + 单元中心 (矢量) ----
        cdata(1)%name  = 'cell_id'
        cdata(1)%ncomp = 1
        cdata(2)%name  = 'centroid'
        cdata(2)%ncomp = 3
        allocate(cdata(1)%val(ncells, 1), cdata(2)%val(ncells, 3))

        n1 = 1                                          ! 第 1 个单元的首节点号
        do c = 1, ncells
            n2 = offset(c)
            ncell = n2 - n1 + 1
            xc = 0.0_dp;  yc = 0.0_dp;  zc = 0.0_dp
            do n = n1, n2
                xc = xc + coords(1, conn(n))
                yc = yc + coords(2, conn(n))
                zc = zc + coords(3, conn(n))
            end do
            cdata(1)%val(c,1) = dble(c)
            cdata(2)%val(c,:) = [xc, yc, zc]/dble(ncell)
            n1 = n2 + 1
        end do

        ! ---- 5. 输出 ----
        call vtu_write('mixed_ascii.vtu',  coords, conn, offset, cell_type, &
                       point_data=pdata, cell_data=cdata, binary=.false.)
        call vtu_write('mixed_binary.vtu', coords, conn, offset, cell_type, &
                       point_data=pdata, cell_data=cdata, binary=.true.)

        print '(a,i0,a,i0,a)', '算例2 混合单元:     ', npts, ' 节点 / ', ncells, &
                               ' 单元 -> mixed_ascii.vtu, mixed_binary.vtu'
    end subroutine run_mixed_cells


end program test
