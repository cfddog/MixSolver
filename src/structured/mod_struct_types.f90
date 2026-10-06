!===============================================================================
! mod_struct_types.f90 -- OpenCFD-EC core derived types
!   BC_MSG_TYPE : boundary/interface connection message
!   Block_TYPE  : per-block geometry and flow-field data
!   Mesh_TYPE   : collection of blocks forming one mesh level
! Split verbatim (phase 2a) from pristine OpenCFD-EC sub_modules.f90.
!===============================================================================
 module Mod_Type_Def
   use precision_EC
   implicit none

    TYPE BC_MSG_TYPE              ! 边界链接信息
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
	 integer:: L1,L2,L3                     ! 子面号，连接顺序描述符
   END TYPE BC_MSG_TYPE

!------------------------------------网格块--------------------------------------
   TYPE Block_TYPE                                 ! 数据结构：网格块 ；包含几何变量及物理变量的信息 
     integer::  Block_no,mpi_id           ! 块号；所属的进程号
	 integer::  nx,ny,nz                   ! 网格数nx,ny,nz
	 integer::  subface                   ! 子面数
 !   几何量  
     real(PRE_EC),pointer,dimension(:,:,:):: x,y,z     ! coordinates of vortex, 网格节点坐标
     real(PRE_EC),pointer,dimension(:,:,:):: xc,yc,zc  ! coordinates of cell center, 网格中心坐标 
     real(PRE_EC),pointer,dimension(:,:,:):: Vol,Si,Sj,Sk ! Volume and surface area, 控制体的体积，i,j,k方向控制体边界面的面积 
	 real(PRE_EC),pointer,dimension(:,:,:):: ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3  !  i,j,k方向三个控制面的法方向
!      Jocabian 变换系数	 
	 real(PRE_EC),pointer,dimension(:,:,:):: ix1,iy1,iz1,jx1,jy1,jz1,kx1,ky1,kz1
	 real(PRE_EC),pointer,dimension(:,:,:):: ix2,iy2,iz2,jx2,jy2,jz2,kx2,ky2,kz2
	 real(PRE_EC),pointer,dimension(:,:,:):: ix3,iy3,iz3,jx3,jy3,jz3,kx3,ky3,kz3
	 real(PRE_EC),pointer,dimension(:,:,:):: ix0,iy0,iz0,jx0,jy0,jz0,kx0,ky0,kz0
!     物理量
	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U,Un,Un1    ! 守恒变量 (本时间步及前一、二个时间步的值), conversation variables 
     real(PRE_EC),pointer,dimension(:,:,:,:) :: Res         ! 残差 （净通量）
     real(PRE_EC),pointer,dimension(:,:,:):: dt          ! (局部)时间步长
     real(PRE_EC),pointer,dimension(:,:,:,:) :: QF      ! 强迫函数 (多重网格法中粗网格使用)
     real(PRE_EC),pointer,dimension(:,:,:,:) :: deltU   ! 守恒变量的差值, dU=U(n+1)-U(n)  多重网格使用
     real(PRE_EC),pointer,dimension(:,:,:,:) :: DU      ! U(n+1)-U(n)    LU-SGS中使用
	 real(PRE_EC),pointer,dimension(:,:,:):: dw         ! turbulent viscous ; distance to the wall  (used in SA model)
	 real(PRE_EC),pointer,dimension(:,:,:):: surf1,surf2,surf3,surf4,surf5,surf6  ! 边界处的通量，供计算气动力、热使用；
	 real(PRE_EC),pointer,dimension(:,:,:):: mu,mu_t    ! 层流粘性系数和湍流粘性系数
     real(PRE_EC),pointer,dimension(:,:,:):: dtime_mesh   ! 时间步长因子 （根据网格质量情况）

	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U_average ! 时间平均场 （d,u,v,w,T）
	 
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! 边界链接信息 
     integer,pointer,dimension(:,:,:):: BcI,BcJ,BcK     ! 边界指示符 （物理边界 or 内边界）
   	 
	 integer:: IFLAG_FVM_FDM               ! 差分 or 有限体积
	 integer:: IF_OverLimit                          ! 物理量超限（如负温度）， 需要降低精度（1阶）
	End TYPE Block_TYPE  

!---------------------------网格 -------------------------------------------------------- 
!  (如单重网格，只有1套；如多重网格，可以有多套) 
  
   TYPE Mesh_TYPE                     ! 数据结构“网格”； 包含几何变量及物理变量信息
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep          ! 网格编号 (1号为最细网格，2号为粗网格， 3号为更粗网格...)，网格块数，网格数目(本进程中), 时间步 
     integer:: NVAR       ! 变量的数目；基本变量5个+ 0，1或2个附件变量 （BL模型0个，SA模型1个，SST模型2个）； 粗网格不使用湍流模型
	 real(PRE_EC)::  tt                   !  推进的时间
	 real(PRE_EC),pointer,dimension(:)::  Res_max,Res_rms                 ! 最大残差，均方根残差, 推进的时间
	 TYPE (Block_TYPE),pointer,dimension(:):: Block       ! “网格块”  （从属于“网格”）

!                                                       控制参数，用于控制数值方法、通量技术、湍流模型等    
!             这些控制参数从属于“网格”，不同“网格”可以采用不同的计算方法、湍流模型等。	 （例如，粗网格用低精度方法，粗网格不使用湍流模型,...）
!  If_dtime_mesh     是否根据网格情况降低局部时间步长（局部时间步长法有效）	
    integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction, Bound_Scheme
   End TYPE Mesh_TYPE
  
  end module Mod_Type_Def
