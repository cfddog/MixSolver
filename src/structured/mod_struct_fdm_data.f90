!===============================================================================
! mod_struct_fdm_data.f90 -- finite-difference-method (FDM) grid data types
! Split verbatim (phase 2a) from pristine OpenCFD-EC sub_modules.f90.
!===============================================================================
   module FDM_data
    use precision_EC
    implicit none
    TYPE FDM_Block_TYPE                              ! 数据结构：网格块 ；包含几何变量及物理变量的信息 
    real(PRE_EC), pointer,dimension(:,:,:):: ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac   ! Jocabian变换系数
    End TYPE FDM_Block_TYPE  
  
   TYPE FDM_Mesh_TYPE                     ! 数据结构“网格”； 包含几何变量及物理变量信息
	 TYPE (FDM_Block_TYPE),pointer,dimension(:):: Block       ! “网格块”  （从属于“网格”）
   End TYPE FDM_Mesh_TYPE 
  
   TYPE (FDM_Mesh_TYPE),pointer,dimension(:):: FDM_Mesh       ! 主数据 “网格”

   ! Flux splitting and reconstruction selectors for embedded FDM blocks;
   ! relocated from Global_Var (phase 2b s3): only mod_struct_fdm uses them.
   integer, save :: FD_Flux, FD_scheme

   end  module FDM_data
