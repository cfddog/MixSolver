!===============================================================================
! mod_struct_init.f90 -- structured solver: startup, mesh/flow initialization,
! control-parameter reading.
! Encapsulates sub_init + sub_read_parameter (phase 2b, batch B4).
!===============================================================================

module mod_struct_init
   use precision_EC
   ! Configuration/state owned solely by the startup path (read_parameter,
   ! init, init_flow_zero).  Relocated from Global_Var (phase 2b s3) after
   ! verifying no other module references them.
   integer, save :: Iflag_init, Kstep_init_smooth, NUM_THREADS, IF_Walldist, &
                    Pre_Step_Mesh(3)
   ! AoS: angle of slide (sideslip), read via control_ec namelist
   real(PRE_EC), save :: AoS
contains

! 初始化：包括创建数据结构及赋初值
! 对于多重网格，根据上级网格的信息，创建各级网格
!  检查网格是否适用于多重网格 
!  单方向网格数= 2*K+1 可用2重网格，=4*K+1 可用3重网格，=8*K+1 可用4重网格 ...  
!---------------------------------------------------------------------------------
!------------------------------------------------------------------------------     
  subroutine init
   use Global_var
   use mod_struct_fdm, only: init_FDM
   use mod_struct_mpi, only: partation, update_Mesh_Center, Update_coordinate_buffer_onemesh
   use mod_struct_grid, only: Comput_Goemetric_var, Output_mesh_debug, register_bc_interfaces
   use mod_struct_io, only: read_main_Mesh, read_inc, comput_dist_wall
   implicit none
   integer :: i,j,k,m,nx1,ny1,nz1,Num_Block1,ksub,Kmax_grid
   real(PRE_EC),allocatable,dimension(:,:,:):: xc,yc,zc
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
 !--------------------------------------------------------------------
 ! initial of const variables
   Ralfa(1)=1.d0 ;  Ralfa(2)=3.d0/4.d0 ; Ralfa(3)=1.d0/3.d0
   Rbeta(1)=1.d0 ;  Rbeta(2)=1.d0/4.d0 ; Rbeta(3)=2.d0/3.d0
   Rgamma(1)=0.d0;  Rgamma(2)=1.d0/4.d0; Rgamma(3)=2.d0/3.d0
   Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma)
!--------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
   call partation                         ! 区域分割 （确定每块所属的进程）
   allocate( Mesh(Num_Mesh) )             ! 主数据结构： “网格” （其成员是“网格块”）。 Mesh(1)为多重网格中最细的网格，Mesh(2),Mesh(3)为粗、更粗的网格。
   call Creat_main_Mesh                   ! 创建主网格(多重网格中最细的网格) (从网格文件Mesh3d.x)
   call read_main_Mesh                    ! 读入主网格
   call read_inc    !读网格连接信息 (bc3d.inc)
   call Update_coordinate_buffer_onemesh(1)
   call Comput_Goemetric_var(1)
   call update_Mesh_Center(1)    ! 更新中心点的坐标的Ghost 值 （周期条件使用）

   call register_bc_interfaces   ! register gridgen generic:8 coupling faces (2b s4b)

   if(IF_Debug==1) call Output_mesh_debug                  ! 输出含虚网格的网格
  
    if(IF_Walldist ==  1)  then
	   call comput_dist_wall   ! 计算(或读取)到壁面的距离
    else
	   if(my_id .eq. 0) print*, " Need not read wall_dist.dat "
	endif

   if(Num_Mesh .ge. 2) then
     call Creat_Mesh(1,2)                 ! 根据1号网格（最细网格）信息，创建2号网格（粗网格）
   endif
   if(Num_Mesh .ge. 3) then
     call Creat_Mesh(2,3)                 ! 根据2号网格信息（粗网格）， 创建3号网格（最粗网格）
   endif
   
   do m=1,Num_Mesh
     call set_BcK(m)   ! 设定边界指示符 (2013-11)
   enddo 

! !!! FVM_FDM FVM_FDM !!!!         hybrid Finite-Difference/ Finite-Valume Method   
   call init_FDM
! !!!----------------------------------------------------------------------------
  end   


!--------------------------------------------------------------------------------------
!   创建数据结构： 最细网格 （储存几何量及守恒变量）
  subroutine Creat_main_Mesh
   use Global_var
   implicit none
   integer:: m,Num_Cell,ierr
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   
!   print*, "-----------------------------"
! ---------node Coordinates----------------------------------------  
!  网格文件：PLOT3D格式；   
   MP=>Mesh(1)
   MP%NVAR=NVAR   ! 最密网格上的变量数目
   MP%Num_Block=Num_Block                      ! 本mpi进程包含的块数
   allocate(MP%Block(Num_block))               ! 创建“网格块”
   call set_size_blocks                        ! 设定每块的大小
   call allocate_mem_Blocks(1)                 ! 给每块的成员数组开辟内存

!  设定参数初值
   allocate(MP%Res_max(NVAR),MP%Res_rms(NVAR))   ! 最大残差与均方根残差 
   MP%Kstep=0
   MP%tt=0.d0
   Num_Cell=0   ! 网格点数
    do m=1,Num_Block
    B => MP%Block(m)
    Num_Cell=Num_Cell+(B%nx-1)*(B%ny-1)*(B%nz-1)
	
	B%IF_OverLimit=0                       ! 限制流场标志 （降低精度等）

    enddo
    call MPI_ALLREDUCE(Num_Cell,MP%Num_Cell,1,MPI_INTEGER,MPI_SUM,Struct_Comm,ierr)
	if (my_id .eq. 0) then
	  print*, "creat main mesh OK, Num_Cell=",MP%Num_Cell
    endif
  end   subroutine Creat_main_Mesh


!---------------------------------------------
! 设定每块的大小(B%nx,B%ny,B%nz), 根据网格文件Mesh3d.x 
  subroutine set_size_blocks
   use Global_var
   implicit none
   integer:: m,mb
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   integer:: NB,NB1,k,ierr
!   integer,allocatable,dimension(:):: NI,NJ,NK
  
   MP=>Mesh(1)
   NB=Total_Block
   allocate(bNI(NB),bNJ(NB),bNK(NB))         ! 块的维数（全局）
  
  if(my_id .eq. 0) then      ! 主进程进行读写操作

   if( Mesh_File_Format .eq. 1) then   ! 格式文件
     open(99,file="Mesh3d.x")
     read(99,*) NB1   ! Block number
     read(99,*) (bNI(k), bNJ(k), bNK(k), k=1,NB)
    else                                ! 无格式文件
     open(99,file="Mesh3d.x",form="unformatted")
     read(99) NB1                       ! 总块数
      if(NB1 .ne. NB) then 
	     print*, "Warning !!! Block number Error !!!"
         print*, "please check 'partation.dat' ..."
		 stop
	  endif
	 read(99) (bNI(k), bNJ(k), bNK(k), k=1,NB)
    endif
    close(99)
   endif
   
   call MPI_bcast(bNI,NB,MPI_Integer,0,  Struct_Comm,ierr)
   call MPI_bcast(bNJ,NB,MPI_Integer,0,  Struct_Comm,ierr)
   call MPI_bcast(bNK,NB,MPI_Integer,0,  Struct_Comm,ierr)
   
   do m=1,MP%Num_Block   ! 本进程包含的块数
    B=> MP%Block(m)       ! 本块
    mb= my_blocks(m)      ! 块号
    B%block_no=mb         ! 块号  
    B%mpi_id=my_id        ! 本块的进程号
    B%nx=bNI(mb)
    B%ny=bNJ(mb)
    B%nz=bNK(mb)
    B%IFLAG_FVM_FDM=Method_FVM   ! 默认有限体积法 
	B%IF_OverLimit=0

   enddo
!   deallocate(NI,NJ,NK)
!   print*, " define size ok ...", my_id

  end subroutine set_size_blocks

!------读入网格信息 (如采用多重网格，则为最密的网格)----------------------------------
  



! 根据上级网格信息，创建新网格m2 
  subroutine Creat_Mesh(m1,m2)
   use Global_Var
   use mod_struct_mpi, only: update_Mesh_Center, Update_coordinate_buffer_onemesh
   use mod_struct_grid, only: Comput_Goemetric_var
   implicit none
   integer:: NB,NVAR1,m,m1,m2,ksub,nx,ny,nz,i,j,k,i1,j1,k1,Bsub,Num_Cell,ierr
   Type (Block_TYPE),pointer:: B1,B2
   TYPE (BC_MSG_TYPE),pointer:: Bc1,Bc2
   Type (Mesh_TYPE),pointer:: MP1,MP2
   MP1=>Mesh(m1)             ! 上一级网格 （细网格）
   Mp2=>Mesh(m2)             ! 本级网格   （粗网格）
   
   MP2%NVAR=5                ! 粗网格上的变量数目

   NB=MP1%Num_Block
   MP2%Num_Block=NB          !  网格m2与m1 块数相同
   MP2%Mesh_no=m2            ! 网格号
   MP2%Num_Cell=0      
   NVAR1=MP2%NVAR    ! NVAR1=5  粗网格不使用湍流模型
   allocate(MP2%Res_max(NVAR1),MP2%Res_rms(NVAR1)) 
   allocate(MP2%Block(NB))   ! 在MP2中创建数据结构：“块”
     Num_Cell=0
   do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)
	 B2%Block_no=B1%block_no
     nx=(B1%nx-1)/2+1        ! 粗网格的点数
	 ny=(B1%ny-1)/2+1
	 nz=(B1%nz-1)/2+1
	 B2%nx=nx
	 B2%ny=ny
	 B2%nz=nz    
	 Num_Cell=Num_Cell+(nx-1)*(ny-1)*(nz-1)          ! 统计MP2的总网格单元数
     B2%IFLAG_FVM_FDM=Method_FVM   ! 默认有限体积法 
 	 B2%IF_OverLimit=0           ! 物理量超限，降低精度
  
    enddo
    call MPI_ALLREDUCE(Num_Cell,MP2%Num_Cell,1,MPI_INTEGER,MPI_SUM,Struct_Comm,ierr)


!    创建几何量及物理量
!--------------------------------------------------
     call allocate_mem_Blocks(m2)           ! 开辟内存
!------------------------------------------------
    do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)
!   设定坐标信息（根据粗、细网格的对应关系）
     do k=1,B2%nz     
	   do j=1,B2%ny
	     do i=1,B2%nx
	       i1=2*i-1 ; j1=2*j-1 ;k1=2*k-1
	       B2%x(i,j,k)=B1%x(i1,j1,k1)         !粗网格与细网格的对应关系 （隔一个点设置一个粗网格点）
           B2%y(i,j,k)=B1%y(i1,j1,k1)
		   B2%z(i,j,k)=B1%z(i1,j1,k1)
	     enddo
	   enddo
	 enddo
	 enddo
!-----------------------------------------------  
    do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)

!    创建连接信息
     Bsub=B1%subface        ! 子面数
     B2%subface=Bsub
     allocate(B2%bc_msg(Bsub))
     do ksub=1, Bsub
	   Bc1=> B1%bc_msg(ksub)    ! 上一级网格的连接信息
	   Bc2=> B2%bc_msg(ksub)    ! 本级网格的连接信息
      
	   Bc2%ib=(Bc1%ib-1)/2+1   ! 粗、细网格下标的对应关系
	   Bc2%ie=(Bc1%ie-1)/2+1   ! 粗、细网格下标的对应关系
	   Bc2%jb=(Bc1%jb-1)/2+1
	   Bc2%je=(Bc1%je-1)/2+1
	   Bc2%kb=(Bc1%kb-1)/2+1
	   Bc2%ke=(Bc1%ke-1)/2+1
 	   
	   Bc2%ib1=(Bc1%ib1-1)/2+1   ! 粗、细网格下标的对应关系
	   Bc2%ie1=(Bc1%ie1-1)/2+1   ! 粗、细网格下标的对应关系
	   Bc2%jb1=(Bc1%jb1-1)/2+1
	   Bc2%je1=(Bc1%je1-1)/2+1
	   Bc2%kb1=(Bc1%kb1-1)/2+1
	   Bc2%ke1=(Bc1%ke1-1)/2+1
      
	   Bc2%bc=Bc1%bc            ! 边界条件 （-1 为内边界）
	   Bc2%face=Bc1%face        ! 面类型(1-6分别代表 i-,j-,k-,i+,j+,k+)
	   Bc2%f_no=Bc1%f_no        ! 子面号
	   Bc2%nb1=Bc1%nb1          ! 连接块
	   Bc2%face1=Bc1%face1      ! 连接面的类型(1-6)
	   Bc2%f_no1=Bc1%f_no1      ! 连接面的子面号
	   Bc2%L1=Bc1%L1            ! 连接方式描述
	   Bc2%L2=Bc1%L2
	   Bc2%L3=Bc1%L3
     enddo
   enddo

   call Update_coordinate_buffer_onemesh(m2)
   call Comput_Goemetric_var(m2)

   call update_Mesh_Center(m2)

   Mesh(m2)%Kstep=0
   Mesh(m2)%tt=0.d0
  
  end  subroutine Creat_Mesh

! ------------------------------------------------------------------------------

  subroutine Init_flow
    use Global_var
    use mod_struct_io, only: read_flow_data
    implicit none
    integer:: i,j,k,m,m1,NVAR1
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B	 
  
   if(Iflag_init .le. 0) then
	 call init_flow_zero                   ! 从均匀场（自由流或静止流场）算起 （先从粗网格计算，再插值到细网格）
   else
	 call read_flow_data
   endif
 
 !    n及n-1时刻流场 （初始时刻设为相同）, 双时间步LU-SGS使用 （仅支持单重网格）
    if(Time_Method .eq. Time_dual_LU_SGS) then             
      MP=> Mesh(1)
      NVAR1=MP%NVAR
      do m=1,MP%Num_Block
      B => MP%Block(m)                
	   do k=-1,B%nz+1
        do j=-1,B%ny+1
         do i=-1,B%nx+1
          do m1=1,NVAR1
           B%Un(m1,i,j,k)=B%U(m1,i,j,k)
		   B%Un1(m1,i,j,k)=B%U(m1,i,j,k)
		  enddo
		 enddo
		enddo
	   enddo
	  enddo		   
    endif   
  end subroutine Init_flow




!--------------------------------------------------------------------------------
! 用来流初始化； 多重网格情况下，从最粗网格开始计算（然后插值到细网格）  
  subroutine init_flow_zero
   use Global_var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   use mod_struct_solver, only: NS_Time_advance, output_Res, prolong_U
   use mod_struct_io, only: smoothing_oneMesh
   implicit none
   real(PRE_EC):: d0,u0,v0,w0,p0,T0,vx,tmp,pin0
   integer:: i,j,k,m,step,nMesh,n
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B	 
!-------------------------------------------------------------------
   if(my_id .eq. 0) then      
    if( Iflag_init .eq. Init_By_FreeStream)then    ! 用自由来流初始化
      print*, "Initial by Free-stream flow ......"
    else  if( Iflag_init .eq. Init_By_Zeroflow) then
      print*, "Initial by Zero flow ......"                                          ! 用静止流场初始化
    endif
   endif

  
   MP=> Mesh(Num_Mesh)   ! Mesh(Num_Mesh) 是最粗的网格
   do m=1,MP%Num_Block
     B => MP%Block(m)                
     ! External-flow initialization (free-stream or zero flow)

	  d0=1.d0
      p0=1.d0/(gamma*Ma*Ma)
     if( Iflag_init .eq. Init_By_FreeStream)then    ! 用自由来流初始化
       u0=cos(A_alfa)*cos(A_beta)
       v0=sin(A_alfa)*cos(A_beta) 
       w0=sin(A_beta)
     else  if( Iflag_init .eq. Init_By_Zeroflow) then                                          ! 用静止流场初始化
 	   u0=0.d0
       v0=0.d0
       w0=0.d0
     endif


	   do k=1-LAP,B%nz+LAP-1
       do j=1-LAP,B%ny+LAP-1
       do i=1-LAP,B%nx+LAP-1

           B%U(1,i,j,k)=d0
           B%U(2,i,j,k)=d0*u0
           B%U(3,i,j,k)=d0*v0
           B%U(4,i,j,k)=d0*w0
           B%U(5,i,j,k)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0+w0*w0)

           if(MP%NVAR .eq. 6) then
! see:       http://turbmodels.larc.nasa.gov/spalart.html
		     B%U(6,i,j,k)=5.d0                 ! 设定为层流粘性系数的5倍 （0.98c以后版本）

		   else if (MP%NVAR .eq. 7) then
		     B%U(6,i,j,k)=10.d0*Kt_Inf   ! 湍动能 （初值设置为来流的10倍）
			 B%U(7,i,j,k)=Wt_Inf         ! 比耗散率
           endif
       enddo
       enddo
	   enddo
     

   enddo

   call Boundary_condition_onemesh(Num_Mesh)     ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(Num_Mesh)          ! 同步各块的交界区

!  Initial smoothing    ! 初始光顺 (迭代Kstep_Init_Smooth步）
   do n=1,Kstep_init_smooth
     if(my_id .eq. 0 .and. mod(n,10) .eq. 0) print*, "Initial smoothing",n
	 call smoothing_oneMesh(Num_Mesh,Smooth_2nd)     
   enddo
   
   
   
!-----------------------------------------------------------------
!------------------------------------------------------
!   准备初值的过程
!   从最粗网格计算，逐级插值到细网格
   do nMesh=Num_Mesh,1,-1 
     do step=1, Pre_Step_Mesh(nMesh)   
       call NS_Time_advance(nMesh)
       if(mod(step,Kstep_show) .eq. 0) call output_Res(nMesh)
     enddo
!     call output (nMesh)
     if(nMesh .gt. 1) then
       call prolong_U(nMesh,nMesh-1,1)                  ! 把nMesh重网格上的物理量插值到上一重网格; flag=1 插值U本身
	   call Boundary_condition_onemesh(nMesh-1)         ! 边界条件 （设定Ghost Cell的值）
	   call update_buffer_onemesh(nMesh-1)              ! 同步各块的交界区 
       print*, " Prolong  to mesh ", nMesh-1, "   OK"           
     endif
   enddo

  end subroutine init_flow_zero 


!----------------------------------------------------


  subroutine allocate_mem_Blocks(nMesh)
   use Global_var
   implicit none
   integer:: nMesh,m,nx,ny,nz,NVAR1
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B

    MP=>Mesh(nMesh)
    NVAR1=MP%NVAR

   do m=1,MP%Num_Block
	 B=>MP%Block(m)

     nx=B%nx ; ny= B%ny  ; nz= B%nz

	  	 
!   申请内存   (x,y,z) 节点坐标； (xc,yc,zc)网格中心坐标; Vol 控制体体积； U, Un 守恒变量
     allocate( B%x(0:nx+1,0:ny+1,0:nz+1), B%y(0:nx+1,0:ny+1,0:nz+1),B%z(0:nx+1,0:ny+1,0:nz+1))   ! 格点坐标
  
!  格心坐标， LAP层虚网格  （便于使用周向周期条件）          
   	 allocate( B%xc(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1) , &
	           B%yc(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1) , &
			   B%zc(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1)  )   ! LAP 层虚网格   
	 
	 
	 allocate( B%U(NVAR1,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1) )   ! LAP 层虚网格   ! bug is removed
	 allocate( B%Un(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))  ! 双层 Ghost Cell
     allocate( B%Vol(nx-1,ny-1,nz-1)) 
     allocate( B%Res(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))        !  残差
	 allocate( B%dt(-1:nx+1,-1:ny+1,-1:nz+1))               !  时间步长
	 allocate( B%deltU(5,-1:nx+1,-1:ny+1,-1:nz+1))          !  两时间步U的差值 （多重网格使用，从粗网格插值而来，5个变量）；
	 allocate( B%dU(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))         !  两时间步U的差值 （LU-SGS中使用）
     allocate( B%Si(nx,ny,nz), B%Sj(nx,ny,nz), B%Sk(nx,ny,nz) )  ! 表面积
     allocate( B%ni1(nx,ny,nz),B%ni2(nx,ny,nz),B%ni3(nx,ny,nz), & 
               B%nj1(nx,ny,nz),B%nj2(nx,ny,nz),B%nj3(nx,ny,nz), &
               B%nk1(nx,ny,nz),B%nk2(nx,ny,nz),B%nk3(nx,ny,nz))
	 allocate( B%dw(nx-1,ny-1,nz-1))         ! 到壁面的距离
!  Jocabian变换系数，用来计算粘性项中的导数 
     allocate(B%ix1(nx,ny,nz),B%iy1(nx,ny,nz),B%iz1(nx,ny,nz), &
	          B%jx1(nx,ny,nz),B%jy1(nx,ny,nz),B%jz1(nx,ny,nz), &
              B%kx1(nx,ny,nz),B%ky1(nx,ny,nz),B%kz1(nx,ny,nz))
     allocate(B%ix2(nx,ny,nz),B%iy2(nx,ny,nz),B%iz2(nx,ny,nz), &
	          B%jx2(nx,ny,nz),B%jy2(nx,ny,nz),B%jz2(nx,ny,nz), &
              B%kx2(nx,ny,nz),B%ky2(nx,ny,nz),B%kz2(nx,ny,nz))
     allocate(B%ix3(nx,ny,nz),B%iy3(nx,ny,nz),B%iz3(nx,ny,nz), &
	          B%jx3(nx,ny,nz),B%jy3(nx,ny,nz),B%jz3(nx,ny,nz), &
              B%kx3(nx,ny,nz),B%ky3(nx,ny,nz),B%kz3(nx,ny,nz))
     allocate(B%ix0(nx,ny,nz),B%iy0(nx,ny,nz),B%iz0(nx,ny,nz), &
	          B%jx0(nx,ny,nz),B%jy0(nx,ny,nz),B%jz0(nx,ny,nz), &
              B%kx0(nx,ny,nz),B%ky0(nx,ny,nz),B%kz0(nx,ny,nz))
!-------------------------------------------------------
      allocate(B%dtime_mesh(nx-1,ny-1,nz-1))   ! 时间步长因子 （由网格质量决定）
               B%dtime_mesh(:,:,:)=1.d0                ! 初值
	 
	 if(Time_Method .eq. Time_Dual_LU_SGS) then             ! 双时间LU_SGS使用 n-1时刻的物理量
          allocate(B%Un1(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))
	 endif

!-------------------------------------------------------

     if(If_viscous .eq. 1) then 
      allocate(B%mu(-1:nx+1,-1:ny+1,-1:nz+1))   ! 层流粘性系数
	  allocate(B%mu_t(-1:nx+1,-1:ny+1,-1:nz+1))    ! 湍流粘性系数
	  B%mu(:,:,:)=1.d0/Re
	  B%mu_t(:,:,:)=0.d0
     endif
!------表面力  (计算总体气动力时使用)------------------------------------
     allocate( B%Surf1(ny,nz,3),B%Surf2(nx,nz,3), B%Surf3(nx,ny,3),  &
               B%Surf4(ny,nz,3),B%Surf5(nx,nz,3), B%Surf6(nx,ny,3)   )          
	 B%Surf1(:,:,:)=0.d0; B%Surf2(:,:,:)=0.d0; B%Surf3(:,:,:)=0.d0
     B%Surf4(:,:,:)=0.d0; B%Surf5(:,:,:)=0.d0; B%Surf6(:,:,:)=0.d0
!---------------------------------------------------------------------------
! --------变量清零（总能、密度设置为1，其余清零）----------------------------
     B%x(:,:,:)=0.d0; B%y(:,:,:)=0.d0; B%z(:,:,:)=0.d0; B%xc(:,:,:)=0.d0; B%yc(:,:,:)=0.d0; B%zc(:,:,:)=0.d0 
     B%U(1,:,:,:)=1.d0; B%U(2,:,:,:)=0.d0; B%U(3,:,:,:)=0.d0; B%U(4,:,:,:)=0.d0; B%U(5,:,:,:)=1.d0
	 
	 if(NVAR1 .eq. 6) then
	    B%U(6,:,:,:)=1.d0
	 else if (NVAR1 .eq. 7) then
	    B%U(6,:,:,:)=0.d0
		B%U(7,:,:,:)=1.d0
	 endif
	   B%Res(:,:,:,:)=0.d0
   
    if( nMesh .ne. 1) then
       allocate( B%QF(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))         ! 强迫函数
	    B%QF(:,:,:,:)=0.d0                                    ! 强迫函数初始化为0
    endif
 
 !----边界指示符 (1 物理边界，0内边界)
     allocate(B%BcI(ny-1,nz-1,2),B%BcJ(nx-1,nz-1,2),B%BcK(nx-1,ny-1,2))
     B%BcI(:,:,:)=0
	 B%BcJ(:,:,:)=0
	 B%BcK(:,:,:)=0
    
   enddo

  end subroutine allocate_mem_Blocks


 ! 设定边界指示符 (0 物理边界， 1 内边界), 用于 高阶格式 （是否启用边界格式）
   subroutine set_BcK(nm)  
   use Global_var
   implicit none
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer :: i,j,k,m,nm,ksub
   integer:: ib,ie,jb,je,kb,ke
     do m=1,Mesh(nm)%Num_Block
       B=>Mesh(nm)%Block(m)
       B%BcI(:,:,:)=0
	   B%BcJ(:,:,:)=0
	   B%BcK(:,:,:)=0
     do  ksub=1,B%subface
     Bc=> B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      
     if(Bc%bc >=0 ) then   ! 非内边界
       if(Bc%face .eq. 1 ) then   
         B%BcI(jb:je-1,kb:ke-1,1)=1
	   else if(Bc%face .eq. 2) then
         B%BcJ(ib:ie-1,kb:ke-1,1)=1
	   else if(Bc%face .eq. 3) then
         B%BcK(ib:ie-1,jb:je-1,1)=1
       else if(Bc%face .eq. 4 ) then   
         B%BcI(jb:je-1,kb:ke-1,2)=1
	   else if(Bc%face .eq. 5) then
         B%BcJ(ib:ie-1,kb:ke-1,2)=1
	   else if(Bc%face .eq. 6) then
         B%BcK(ib:ie-1,jb:je-1,2)=1
       endif
	 endif
	enddo
    enddo
   end
!------------------------------------------------

!  -----------------------读取流动参数及控制变量----------------------
  subroutine read_parameter(ctlfile)
   use Global_var
   implicit none
   ! Optional path to the control file (namelist format).  When absent the
   ! historical default "control.ec" in the current directory is used.
   character(len=*), intent(in), optional :: ctlfile
   character(len=:), allocatable :: cfile
   logical ext1
!-----------------------------------------------------------------------------
   if ( present(ctlfile) ) then
      cfile = trim(ctlfile)
   else
      cfile = "control.ec"
   end if

   call set_default_parameter ! 设置参数默认值

   if(my_id .eq. 0) then
     inquire(file=cfile,exist=ext1)
      if(ext1) then           ! 优先读取control.ec (Namelist 格式控制文件)
       call read_parameter_ec(cfile)
      else
       print*, "Can not find '" // cfile // "', stop !"
	   stop
	  endif


   endif

   call bcast_para      ! 广播至全部进程
   
   call set_const_para
!--------------------------------------------------------------------

  end subroutine read_parameter
   
!--------------------------------------------------------------------    
! 设置参数的默认值  
  subroutine set_default_parameter 
   use Global_var
   implicit none
    Ma=1.d0        ! Mach number     
    Re= 1000.d0    ! Reynolds number
    gamma=1.4d0    ! 

	AOA=0.d0       ! Angle of attack    
	AOS=0.d0       ! Angle of Slide
	P_outlet=-1.d0   ! Outlet pressure (<0 extrapolation)  
	t_end=100.d0     ! End time (non-dimensional)
	Kstep_save=1000  ! Save data per xxx steps
    Iflag_turbulence_model=0   ! turbulence model (0 none, 1 BL, 2 SA, 3 SST)
	Iflag_init=0  ! 0 从初始值（均匀来流）开始计算； 1  续算； -1 从0 流场开始计算
    If_viscous=1  ! 0 无粘； 1 有粘
    Iflag_local_dt=1   ! 0 全局步长；  1 局部时间步长
    dt_global=0.01     ! Global time step
    CFL=1.d0           ! CFL number 
    dtmax=10.d0       !Limit of maximum time step 
	dtmin=1.d-9       !Limit of minimum time step
	Time_Method=0     ! Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
    If_Residual_smoothing=0   ! 0 Do not need smoothing
    w_LU=1.d0         ! factor in LU-SGS
    If_dtime_mesh=1   ! decrease time step when grid quality is not good
    Iflag_Scheme= 5    ! 0 UD1;  1 NND2 ; 2 UD3 ; 3 MUSCL2U ;  4 MUSCL2C; 5 MUSCL3; 6 OMUSCL2; 7 WENO5 ; 8 UD5; 9 WENO7 
    Iflag_Flux=5      ! 1 Steger_Warming; 2 HLL, 3 HLLC,  4 Roe, 5  Van_Leer, 6 Ausm
    IFlag_Reconstruction=0   !  0 Original, 1 Conservative, 2 Characteristic
    Mesh_File_Format=0    ! 0 unformatted, 1 formatted
    Kstep_show=1          ! show per xxx steps
    Kstep_average=0       ! do not average
    Kstep_smooth=-1       ! Smoothing flow per xxx steps (<0 donot smooth)
    Kstep_init_smooth=0   ! Smoothing flow at initial time
    Num_Mesh=1            ! Number of mesh (1 single-grid), 2,3 multi-grid
    T_inf=288.15d0        ! Reference temperature  (in viscous coefficient)
    Twall=-1.d0           ! Wall temperature (in K degree) (<0 adabitic)
    Kt_inf=1.d-5          ! initial of Kt for SST model
    Wt_inf=0.01d0         ! Initial of Wt for SST model
    IF_Debug=0            ! 1 Debug mode
    NUM_THREADS=1         ! Threads for OpenMP
    Step_Inner_Limit=20   ! Inner time advance limit (Dual-time)
    Res_Inner_Limit=1.d-10 ! Inner Resdial limit (Dual-time)
    MUT_MAX=-1.d0          ! Limit for vt/vs (<0 no limit)
	Bound_Scheme= Scheme_MUSCL2C         ! Boundary scheme (Default: MUSCL2C)
    Pre_Step_Mesh(1:3)=0   ! Pre step for multi-grid
    Ref_S=1.d0             ! Ref. area
	Ref_L=1.d0             ! Ref. length
	Centroid(1:3)=0.d0  ! Centroid coordinate  
    Cood_Y_UP=1   !       默认Y轴垂直向上
	  IFLAG_LIMIT_FLOW=1         ! 限制流场（密度、速度、压力），设定为1 ， 2016-10-21
	Pdebug(1:4)=1
    PrL=0.7d0          ! Linear Prandtl number
	PrT=0.9d0          ! Turbulent Prandtl number
    Ldmin=1.d-6
	Ldmax=1000.d0
	Lpmin=1.d-6
	Lpmax=1000.d0
	Lumax=1000.d0
	LSAmax=1000.d0
	CP1_NSA=0.2d0   ! for New SA
	CP2_NSA=100.d0
	Periodic_dX=0.d0   ! 周期边界的几何增量
	Periodic_dY=0.d0 
	Periodic_dZ=0.d0

    Iflag_savefile=0       ! 默认写入flow3d.dat
    IF_Scheme_Positivity=1     ! 检查插值过程中压力、密度是否非负，否则使用1阶迎风；
end

!------read parameter (Namelist type)---------------- 
  subroutine read_parameter_ec(ctlfile)
   use Global_var
   implicit none
   character(len=*), intent(in), optional :: ctlfile

 	namelist /control_ec/ Ma, Re, AoA, AoS, p_outlet, t_end, &
	    gamma, PrL, PrT, &
	    Kstep_save, &
	    Iflag_turbulence_model,Iflag_init,If_viscous,  &
        Iflag_local_dt,dt_global,CFL,dtmax,dtmin,Time_Method, &
		If_Residual_smoothing,w_LU,If_dtime_mesh,  &
        Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction, &
		Mesh_File_Format,Kstep_show,Kstep_average,Kstep_smooth,Kstep_init_smooth,  &
        Num_Mesh,T_inf,Twall,Kt_inf,Wt_inf,IF_Debug,NUM_THREADS,  &
		Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, Bound_Scheme, &
        Pre_Step_Mesh,Ref_S,Ref_L,Centroid,Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug, &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
        IF_Scheme_Positivity, &
		Periodic_dX, Periodic_dY, Periodic_dZ, &
		Iflag_savefile


   if ( present(ctlfile) ) then
	open(99,file=trim(ctlfile))
   else
	open(99,file="control.ec")
   end if
	read(99,nml=control_ec)
    close(99)
 
 !---- convert parameters ----------------------


 !---output paramters----------------------------
    open(99,file="output_para.out")
	write(99,*) "-------OpenCFD-EC (Ver 1.1t), (c) Li Xinliang, lixl@imech.ac.cn--"
	write(99,*) "-------------------------------------------"
	write(99,*) "Ma=", Ma,  "  Re= ", Re , "gamma=", gamma
	write(99,*) "PrL=", PrL, "PrT=",PrT

	write(99,*) "A_alfa=", A_alfa, " A_beta=", A_beta 
	write(99,*) "P_outlet=", P_outlet
	write(99,*) "t_end=", t_end, " Iflag_local_dt=", Iflag_local_dt
	write(99,*) "dt_global=",dt_global 
	write(99,*) "Time_Method=", Time_Method,  " CFL= ", CFL 
	write(99,*) "dtmax=", dtmax, " dtmin=", dtmin 
    write(99,*) "Iflag_init=", Iflag_init," If_viscous=", If_viscous
	write(99,*) "Iflag_turbulence_model=",Iflag_turbulence_model
    write(99,*) "Iflag_Scheme= ",Iflag_Scheme, " Iflag_Flux=", Iflag_Flux
	write(99,*) "IFlag_Reconstruction=", IFlag_Reconstruction, "Bound_Scheme= ", Bound_Scheme
    write(99,*) "Kstep_save= ",Kstep_save,  " Kstep_show=", Kstep_show, "Kstep_average=",Kstep_average
	write(99,*) "Kstep_smooth= ", Kstep_smooth, " Kstep_init_smooth=",Kstep_init_smooth
    write(99,*) "If_Residual_smoothing=", If_Residual_smoothing, " If_dtime_mesh=",If_dtime_mesh
	write(99,*) "w_LU=",w_LU
    write(99,*) "Mesh_File_Format=",Mesh_File_Format, " Num_Mesh=",Num_Mesh
    write(99,*) "T_inf=", T_inf, " Twall=", Twall
	write(99,*) "Kt_inf=",Kt_inf,  " Wt_inf=",Wt_inf
    write(99,*) "Step_Inner_Limit=",Step_Inner_Limit, " MUT_MAX=", MUT_MAX
	write(99,*) "Pre_Step_Mesh(:)=", Pre_Step_Mesh(1:Num_Mesh)
    write(99,*) "IF_Debug=",IF_Debug
	write(99,*) "Ref_S=", Ref_S, "Ref_L=", Ref_L
	write(99,*) "Centroid=", Centroid(1:3)
	write(99,*) "Cood_Y_UP=",Cood_Y_UP
    write(99,*) "Periodic_dX, dY, dZ=", Periodic_dX,Periodic_dY,Periodic_dZ
	write(99,*) "IFLAG_LIMIT_FLOW=",IFLAG_LIMIT_FLOW
	write(99,*) "Ldmin,Ldmax,Lpmin,Lpmax,Lumax=", Ldmin,Ldmax,Lpmin,Lpmax,Lumax
	write(99,*) "CP1_NSA,CP2_NSA=",CP1_NSA,CP2_NSA
	write(99,*) "IF_Scheme_Positivity =", IF_Scheme_Positivity
	write(99,*) "Iflag_savefile=", Iflag_savefile
	write(99,*) "NUM_THREADS=",NUM_THREADS
    write(99,*) "--------------------------------------------"

    close(99)

 !----------------------------------------------
 end
!--------------------------------------------------
 

!------------------------------------------------------------------
  
  
   subroutine bcast_para
   use Global_var
   implicit none
    integer:: Ipara(100),ierr
    real(PRE_EC):: rpara(100)
    Ipara=0
	rpara=0.d0
!----
	rpara(1)=Ma 
	rpara(2)=Re
	rpara(3)=AoA
	rpara(4)=AoS
	rpara(5)=p_outlet
	rpara(6)=t_end
	rpara(7)=dt_global
	rpara(8)=CFL
	rpara(9)=dtmax
	rpara(10)=dtmin
	rpara(11)=w_LU
	rpara(12)=T_inf
	rpara(13)=Twall
	rpara(14)=Kt_inf
	rpara(15)=Wt_inf
    rpara(16)=Res_Inner_Limit
    rpara(17)=MUT_MAX
    rpara(18)=Ref_S
	rpara(19)=Ref_L
	rpara(20:22)=Centroid(1:3)
	rpara(23)=Ldmin
	rpara(24)=Ldmax
	rpara(25)=Lpmin
	rpara(26)=Lpmax
	rpara(27)=Lumax
	rpara(28)=LSAmax
	rpara(29)=CP1_NSA
	rpara(30)=CP2_NSA
	rpara(31)=gamma
	rpara(32)=PrL
	rpara(33)=PrT
    rpara(36)=Periodic_dX
	rpara(37)=Periodic_dY
	rpara(38)=Periodic_dZ




    Ipara(1)=Kstep_save
	Ipara(2)=Iflag_turbulence_model
	Ipara(3)=Iflag_init
	Ipara(4)=If_viscous
    Ipara(5)=Iflag_local_dt
	Ipara(6)=Time_Method
	Ipara(7)=If_Residual_smoothing
	Ipara(8)=If_dtime_mesh
    Ipara(9)=Iflag_Scheme
	Ipara(10)=Iflag_Flux
	Ipara(11)=IFlag_Reconstruction
	Ipara(12)=Mesh_File_Format
	Ipara(13)=Kstep_show
	Ipara(14)=Kstep_smooth
	Ipara(15)=Kstep_init_smooth
    Ipara(16)=Num_Mesh
	Ipara(17)=IF_Debug
	Ipara(18)=NUM_THREADS
    Ipara(19)=Step_Inner_Limit
    Ipara(20)=Bound_Scheme
    Ipara(21)=Cood_Y_UP
	Ipara(22)=IFLAG_LIMIT_Flow
    Ipara(23:26)=Pdebug(1:4)
	Ipara(28)=IF_Scheme_Positivity
   	Ipara(30)=Kstep_average
    Ipara(31)=Iflag_savefile

	 call MPI_bcast(rpara,100,OCFD_DATA_TYPE,0,  Struct_Comm,ierr)
	 call MPI_bcast(Ipara,100,MPI_Integer,0,  Struct_Comm,ierr)

	Ma=rpara(1) 
	Re=rpara(2)
	AoA=rpara(3)
	AoS=rpara(4)
	p_outlet=rpara(5)
	t_end=rpara(6)
	dt_global=rpara(7)
	CFL=rpara(8)
	dtmax=rpara(9)
	dtmin=rpara(10)
	w_LU=rpara(11)
	T_inf=rpara(12)
	Twall=rpara(13)
	Kt_inf=rpara(14)
	Wt_inf=rpara(15)
    Res_Inner_Limit=rpara(16)
    MUT_MAX=rpara(17)
    Ref_S=rpara(18)
	Ref_L=rpara(19)
	Centroid(1:3)=rpara(20:22)
	Ldmin=rpara(23)
	Ldmax=rpara(24)
	Lpmin=rpara(25)
	Lpmax=rpara(26)
	Lumax=rpara(27)
	LSAmax=rpara(28)
	CP1_NSA=rpara(29)
	CP2_NSA=rpara(30)
	gamma=rpara(31)
	PrL=rpara(32)
	PrT=rpara(33)
    Periodic_dX=rpara(36)
	Periodic_dY=rpara(37)
	Periodic_dZ=rpara(38)



    Kstep_save=Ipara(1)
	Iflag_turbulence_model=Ipara(2)
	Iflag_init=Ipara(3)
	If_viscous=Ipara(4)
    Iflag_local_dt=Ipara(5)
	Time_Method=Ipara(6)
	If_Residual_smoothing=Ipara(7)
	If_dtime_mesh=Ipara(8)
    Iflag_Scheme=Ipara(9)
	Iflag_Flux=Ipara(10)
	IFlag_Reconstruction=Ipara(11)
	Mesh_File_Format=Ipara(12)
	Kstep_show=Ipara(13)
	Kstep_smooth=Ipara(14)
	Kstep_init_smooth=Ipara(15)
    Num_Mesh=Ipara(16)
	IF_Debug=Ipara(17)
	NUM_THREADS=Ipara(18)
    Step_Inner_Limit=Ipara(19)
    Bound_Scheme=Ipara(20)
    Cood_Y_UP=Ipara(21)
	IFLAG_LIMIT_FLOW=Ipara(22)
    Pdebug(1:4)=Ipara(23:26)
	IF_Scheme_Positivity=Ipara(28)
   	Kstep_average=Ipara(30)
    Iflag_savefile=Ipara(31)


    call MPI_bcast(Pre_Step_Mesh,Num_Mesh,MPI_Integer,0,  Struct_Comm,ierr)
    end  subroutine bcast_para
  
  
      
! 设定常数
   subroutine set_const_para 
   use Global_var
   implicit none

   Twall=Twall/T_inf                         ! wall temperature
   Lpmin=Lpmin/(gamma*Ma*Ma)                 ! rato of free-stream pressure
   Lpmax=Lpmax/(gamma*Ma*Ma)



   if(Bound_Scheme== Scheme_none)  Bound_Scheme=Iflag_Scheme    ! 如不使用边界格式，则与内点格式一致   
   if(If_viscous .eq. 0) Iflag_turbulence_model=Turbulence_NONE !    求解无粘方程，不采用湍流模型
    
   if(Iflag_turbulence_model .eq. Turbulence_SA .or. & 
      Iflag_turbulence_model .eq. Turbulence_NewSA) then
     NVAR=6                       ! SA模型，总共6个变量
   else if (Iflag_turbulence_model .eq. Turbulence_SST) then
     NVAR=7                       ! SST 模型，总共7个变量
   else
     NVAR=5                       ! 5个变量
   endif

  if(Iflag_turbulence_model .eq. 0) then
   IF_Walldist=0                   ! 无需该数据
  else
   IF_Walldist=1
  endif
 
!--------------------------------------------------------------------------------
  AoA=AoA*PI/180.d0   ! Angle of attack
  AoS=AoS*PI/180.d0   ! Angle of Slide
 if(Cood_Y_UP ==1) then   ! Y 轴垂直向上 or Z轴垂直向上 
   A_alfa=AoA             ! Y轴向上， A_alfa为攻角
   A_beta=AoS
 else
   A_alfa=AoS  
   A_beta=AoA
 endif
   
  
   Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma) 
   Cp=Cv*gamma 
!--------------------------------------------------------------------

!  计算到壁面的距离
!   if(Iflag_turbulence_model .eq. Turbulence_SA .or. Iflag_turbulence_model .eq. Turbulence_SST) then
!     inquire(file="wall_dist.dat",exist=file_exist)
!	 if( .not. file_exist) then
!	   call  comput_wall_dist  
!     endif
!   endif
    
	if(Num_Mesh .ne. 1) then
	   if(Time_Method .eq. Time_LU_SGS .or. Time_Method .eq. Time_Dual_LU_SGS) then
	    print*, "In this version, LU_SGS (or Dual_time_LU_SGS) DO NOT Support Multi-Grid !!!"
		stop
	   endif
	endif
  end

end module mod_struct_init
