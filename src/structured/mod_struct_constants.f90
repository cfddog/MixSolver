!===============================================================================
! mod_struct_constants.f90 -- OpenCFD-EC named constants (scheme/flux/BC codes)
! Split verbatim (phase 2a) from pristine OpenCFD-EC sub_modules.f90.
! NOTE: PI/Lim_Zero overlap conceptually with common/mod_constants; they are
! retained unchanged in 2a and will be reconciled in the phase-5 units work.
!===============================================================================
  module const_var
   use precision_EC
   implicit none
   real(PRE_EC),parameter::  PI=3.1415926535897932d0, Lim_Zero=1.d-20   ! 小于该值认为0
   integer,parameter:: LAP=4                          ! 虚网格的数目 （使用3阶格式该值不小于2；如使用5阶WENO, 该值不小于3; 如果使用WENO7, 则该值不小于4）
   integer,parameter:: Scheme_UD1=0, Scheme_NND2=1, Scheme_UD3=2,Scheme_MUSCL2U=3,Scheme_MUSCL2C=4,   &
                       Scheme_MUSCL3=5,Scheme_OMUSCL2=6,Scheme_WENO5=7,Scheme_UD5=8, Scheme_WENO7=9
   integer,parameter:: Scheme_CD2=20, Scheme_none=-1         ! 不使用（边界）格式
   integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_Van_Leer=5,Flux_Ausm=6
   integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2


!---------------边界条件-------------------------
!  integer,parameter:: BC_Wall=-10, BC_Farfield=-20, BC_Periodic=-30,BC_Symmetry=-40,BC_Outlet=-22
   integer,parameter::  BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6 
   integer,parameter::  BC_Periodic=501, BC_Extrapolate=401      ! (扩展) 与Griggen .inp文件的定义可能有所区别，请注意
!                     周期性边界条件并非设置成BC_Peridodic=501, 而是设置BC_PeriodicL=-2, BC_PeriodicR=-3

!   周期边界按照内边界处理， 但几何变量须专门处理。 叶轮机模式下，周期边界也须特殊处理   

!  -1 内边界 （非物理边界）；  -2  左周期边界； -3 右周期边界 （区分左、右便于处理几何坐标）
   integer,parameter::  BC_Inner=-1, BC_PeriodicL=-2, BC_PeriodicR=-3   
   
   integer,parameter::  BC_Zero=0       ! 无边界条件

!                 用户自定义的边界条件，要求代码 >=900  
   integer,parameter::  BC_USER_FixedInlet=901, BC_USER_Inlet_time=902       !给定入口流动； 给定入口时间序列
   integer,parameter:: BC_USER_Blow_Suction_Wall=903    ! 吹吸扰动壁面

   integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
   integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2, Turbulence_SST=3,Turbulence_NewSA=21
   integer,parameter:: Init_continue=1, Init_By_FreeStream=0, Init_By_Zeroflow=-1,  Smooth_2nd=0,Smooth_4th=1 
!   real(PRE_EC), parameter::  Density_LIMIT=1.d-4,Temperature_LIMIT=1.d-4,Pressure_LIMIT=1.d-4

   integer,parameter:: Method_FVM=0, Method_FDM=1        ! 差分、有限体积
   integer,parameter:: FD_WENO5=1,FD_WENO7=2,FD_OMP6=3  ! 差分法采用的数值格式
   integer,parameter:: FD_Steger_Warming=1,FD_Van_Leer=2


  end module const_var
