!===============================================================================
! mod_struct_global.f90 -- OpenCFD-EC global control variables and mesh handles
! 19 of the 31 ported files 'use Global_Var'.  Kept verbatim in phase 2a;
! passing data through arguments instead is deferred to the phase-2b
! modularisation.
! Split verbatim (phase 2a) from pristine OpenCFD-EC sub_modules.f90.
!===============================================================================
  module Global_Var
   use mpi                  ! MPI_COMM_WORLD default for Struct_Comm
   use const_var        ! 常量
   use mod_type_def     ! 边界连接
   implicit none


!---------------------------------------------------------------------------------------------
! global variables                                       各子程序均可见的全局变量
!----------------------------------------------------------------------------

   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh                          ! 主数据 “网格”
   integer,save:: Num_Mesh,NVAR , Total_block, Num_block                      ! 网格的套数 ，变量数， 总网格块数, 本mpi进程的网格块数  
   integer,pointer,dimension(:):: bNi,bNj,bNk                            ! 各块的维数（全局）， bNi(k)为（全局）第k块的nx
   integer,save::  Kstep_save, Iflag_turbulence_model,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method, &
	  Kstep_show,If_viscous,If_Residual_smoothing,Mesh_File_Format,IF_Debug, &
	  Kstep_smooth,If_dtime_mesh, Step_Inner_Limit, &
      Bound_Scheme, &                         ! 边界格式
      IFLAG_LIMIT_FLOW, &                     ! 是否需要限制压力增长率
      IF_Scheme_Positivity,          &    ! 检查插值过程中压力、密度是否非负，否则使用1阶迎风；
      Kstep_average,                 &          ! 时均统计的步数间隔， 0 不统计
      Iflag_savefile                       ! 0 保存到flow3d.dat, 1 保存到flow3d-xxxxxxx.dat

   integer,save:: KRK=0            ! Runge-Kutta方法的 子步
 ! global parameter (for all Meshes )                     流动参数, 对全体“网格”都适用

   real(PRE_EC),save:: Ma,Re,gamma,Cp,Cv,t_end,P_OUTLET,&
                       A_alfa,A_beta,PrL,PrT,T_inf,Twall,w_LU,Kt_inf,Wt_inf , Res_Inner_Limit, MUT_MAX,AoA
   real(PRE_EC),save:: Periodic_dX,Periodic_dY,Periodic_dZ

 
 ! 全局控制参数，控制数值方法、通量技术及湍流模型等 （有些只对最细网格有效）
   real(PRE_EC),save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin    ! RK方法中的常量，与时间步长有关的量
   integer,save:: Cood_Y_UP             ! 1 Y轴垂直向上， 0 Z轴垂直向上
   integer,save:: Pdebug(4)             ! debug 使用， 输出某块的某一点的值
   real(PRE_EC),save:: Ref_S,Ref_L, Centroid(3)                           ! （计算六分量使用）参考面积, 参考长度, 矩心坐标
   real(PRE_EC),save:: Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax                ! 对密度、压力、速度、及SA模型中变量的限制条件 (最小与最大值)
   real(PRE_EC),save:: CP1_NSA,CP2_NSA       ! parameters in New SA model
 !-----------mpi data -----------------------------------------------------------
   integer:: my_id,Total_proc                   ! my_id (本进程号), Total_proc 总进程数
   ! Group communicator for all struct-internal MPI calls.  Defaults to
   ! MPI_COMM_WORLD so the standalone structured solver is unchanged; the
   ! coupling driver (mod_struct_driver) overwrites it with the STRUCT_GROUP
   ! sub-communicator so struct collectives never wait on uns ranks.
   integer:: Struct_Comm = MPI_COMM_WORLD
   integer,pointer,dimension(:):: B_Proc, B_n    ! B_proc(m) m块所在的进程号; B_n(m) m块所在进程中的内部编号
   integer,pointer,dimension(:):: my_Blocks       ! 本进程包含的块号
  end module Global_Var  
