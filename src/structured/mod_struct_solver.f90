!===============================================================================
! mod_struct_solver.f90 -- structured solver: residual / time advance /
! turbulence models (SA, NewSA, SST, BL) / limiter / filtering.
! SCC(6) cycle (Residual <-> time_advance <-> turbulence) forces all of these
! into one module: Fortran forbids cyclic use between modules.
! The legacy filting_Var helper module is kept verbatim above mod_struct_solver.
! Encapsulates sub_Residual + sub_time_advance + sub_turbulence_{SA,NewSA,SST,BL}
! + sub_limitflow + sub_filtering (phase 2b, batch B3).
!===============================================================================

  module filting_Var
   use precision_EC
   real(PRE_EC), save,pointer,dimension(:,:,:,:)::  f,f0 ! 变量
  end module Filting_Var

module mod_struct_solver

contains

 !-----The core subroutines: Comput inviscous and viscous flux ------------------------------
 !     OpenCFD-EC 3D 
 !     Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS. lixl@imech.ac.cn
 !     Code by Li Xinliang and Leng Yan
 !     Ver 0.43 2010-11-28
 !     Ver 0.50 2010-12-9
 !     Ver 0.74 2011-12-29
 !     Ver 0.8  2012-5-5
 !     Ver 0.97a 2013-5-4:  viscous flux code modified, corner points is not used
 !     Ver 1.01  2013-11-13:  boundary scheme can be used 
 !     Ver 1.16  2017-7-11:   如果物理量超限，则本块降为1阶迎风，且关闭粘性项；
  
! 计算残差（网格的全部块）
  Subroutine Comput_Residual_one_mesh(nMesh)
   use Global_Var
   use Flow_Var
   use mod_struct_time, only: comput_Lijk, Residual_smoothing, comput_dt, du_LU_SGS
   use mod_struct_fdm,  only: Residual_FDM
   implicit none
   integer:: nMesh,NVAR1,mBlock,nx,ny,nz,i,j,k,m,KL,IR,JR,KR
   integer,save:: Iflag1=0
   real(PRE_EC):: Sfac,Sfac1
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B

!---------------------------------------------  
  if(Time_Method .eq. Time_Dual_LU_SGS) then
	 if(Iflag1 .eq. 0) then
	   Iflag1=1
	   Sfac=1.d0/(3.d0*dt_global)  ! 第一个时间步， Un=U(n-1), 时间精度1阶
	   Sfac1=1.d0/dt_global
	 else
	   Sfac=1.d0/(2.d0*dt_global) 
	   Sfac1=3.d0/(2.d0*dt_global)       
	 endif
  else
     Sfac=0.d0
	 Sfac1=0.d0
  endif

!----------------------------------------------- 
   MP=>Mesh(nMesh)
   MP%Res_max(:) =0.d0  ! 最大残差
   MP%Res_rms(:) =0.d0  ! 均方根残差            
   NVAR1=MP%NVAR
!----------------------------------------------
   do mBlock=1,MP%Num_Block
     B => MP%Block(mBlock)                  ! 第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny; nz=B%nz
     KL=1-LAP
	 IR=nx+LAP-1
	 JR=ny+LAP-1
	 KR=nz+LAP-1
	 allocate(d(KL:IR,KL:JR,KL:KR),uu(KL:IR,KL:JR,KL:KR),v(KL:IR,KL:JR,KL:KR),  &
              w(KL:IR,KL:JR,KL:KR), T(KL:IR,KL:JR,KL:KR), &
              cc(KL:IR,KL:JR,KL:KR),p(KL:IR,KL:JR,KL:KR))
    
	
	 allocate(Flux(NVAR1,nx,ny,nz))                            ! 通量(流体方程)
     allocate(Lci(nx,ny,nz),Lcj(nx,ny,nz),Lck(nx,ny,nz),Lvi(nx,ny,nz),Lvj(nx,ny,nz),Lvk(nx,ny,nz))  ! 谱半径（无粘、粘性）

!------------------------------------------------------------------------------------     
	 call limit_vt(nMesh,mBlock)    ! 对SA，SST方程的物理量(vt,Kt,Wt)进行限制

	 call comput_duvtpc(nMesh,mBlock)                      ! 计算基本量 d,u,v,T,p,cc , vt,kt,Wt

!---------------------------------------------------------------------------------------
!  求解N-S方程的核心模块： 计算残差（右端项）  
     call Residual (nMesh,mBlock)                          ! 计算一个网格块的残差（右端项） ; 第nMesh 重网格的第mBlock块
!----------------------------------------------
!    OpenCFD-SEC： 采用差分法计算该块的残差
! !!! FVM-FDM --------------
   	 if(B%IFLAG_FVM_FDM  .eq. Method_FDM) call Residual_FDM(nMesh,mBlock)  ! 用有限体积法计算残差（补丁程序）   !!!SEC!!!
!---FVM-FDM-------------





!  双时间步长法，添加附加残差     
	 if(Time_Method .eq. Time_Dual_LU_SGS) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,m)
	  do k=1,nz-1
	  do j=1,ny-1
      do i=1,nx-1
      do m=1,NVAR1
       B%Res(m,i,j,k)=B%Res(m,i,j,k)-(3.d0*B%U(m,i,j,k)-4.d0*B%Un(m,i,j,k)+B%Un1(m,i,j,k))*B%vol(i,j,k)*Sfac
      enddo
	  enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO
     endif



     call comput_Lijk(nMesh,mBlock)                        ! 计算谱半径  (Blazek's book, p189-190), 残差光顺，局部时间步长及LU-SGS中均使用该值
	 if(If_Residual_smoothing .eq. 1 ) then
	     call Residual_smoothing(nMesh,mBlock)             ! 残差光顺
	 endif
	 call comput_dt(nMesh,mBlock)                          ! 计算(当地) 时间步长

     if(Time_Method .eq. Time_LU_SGS .or. Time_Method .eq. Time_Dual_LU_SGS) then
        call du_LU_SGS(nMesh,mBlock,Sfac1)                      ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
     endif

	 deallocate(d,uu,v,w,T,cc,p,Flux)
     deallocate(Lci,Lcj,Lck,Lvi,Lvj,Lvk)


   enddo    
   

  end Subroutine Comput_Residual_one_mesh


!------------------------------------------------------------------------------
! 计算残差 (残差=右端项=净流量); 是求解N-S方程的核心模块; 
! Flux (= inviscous flux + viscous flux )

    subroutine Residual(nMesh,mBlock)
     Use Global_Var
     Use Flow_Var
     implicit none
     integer:: nx,ny,nz,i,j,k,m,nMesh,mBlock !,Scheme,IFlux,Reconstruction
     Type (Block_TYPE),pointer:: B
     TYPE (Mesh_TYPE),pointer:: MP
     MP=> Mesh(nMesh)
     B => MP%Block(mBlock)
     
     if(If_viscous .eq. 1 ) then
      call get_viscous(nMesh,mBlock)                   ! 计算层流粘性系数  
!  湍流模型
     if(MP%Iflag_turbulence_model .eq. Turbulence_BL) then
       call  turbulence_model_BL(nMesh,mBlock)               !  BL模型
     else if(MP%Iflag_turbulence_model .eq. Turbulence_SA) then
       call Turbulence_model_SA(nMesh,mBlock)                  ! SA模型
     else if(MP%Iflag_turbulence_model .eq. Turbulence_NewSA) then
       call Turbulence_model_NewSA(nMesh,mBlock)                  ! New SA模型 (By Li XL, He ZW and Li L)
   
     else if(MP%Iflag_turbulence_model .eq. Turbulence_SST) then
	   call  Turbulence_model_SST(nMesh,mBlock)

     else
       B%mu_t(:,:,:)=0.d0                                ! 层流
     endif 

!    对湍流粘性系数进行限制	 
 	 call limit_mut(nMesh,mBlock)

     endif
     
!-----------------无粘项-------------------------------------
     call flux_inviscous_i(nMesh,mBlock)
     call flux_inviscous_j(nMesh,mBlock)
     call flux_inviscous_k(nMesh,mBlock)
    
! 三个方向的粘性通量
    if(If_viscous .eq. 1) then
     call flux_viscous_i(nMesh,mBlock)
     call flux_viscous_j(nMesh,mBlock)
     call flux_viscous_k(nMesh,mBlock)
    endif
 
    end  
!---------------------------------------------------------------------------------------------------------
!  i方向的无粘通量, 即穿过以(i,j+1/2,k+1/2) (即(I-1/2,J,K)) 为中心侧面的通量 
!  支持OpenMP, since ver0.71 
   subroutine flux_inviscous_i(nMesh,mBlock)
    Use Global_Var
    Use Flow_Var
    use mod_struct_flux, only: Flux_steger_warming_1Da, Flux_Roe_1D, &
                               Flux_Van_Leer_1Da, Flux_Ausmpw_1Da, Flux_HLL_HLLC_1D
    implicit none
    real(PRE_EC),dimension(5):: UL,UR,QL,QR,Flux0
    real(PRE_EC):: U0(1-LAP:LAP,5)
    integer:: mBlock,NVAR1,i,j,k,m,nx1,ny1,nz1,ksub,nMesh
    integer:: Scheme,IFlux,Reconstruction
    real(PRE_EC):: s1x,s1y,s1z,s2x,s2y,s2z,s3x,s3y,s3z,s0,t10,t1x,t1y,t1z  ! 法方向(s1)及切方向(s2,s3)
    real(PRE_EC):: d0,uu0,v0,w0,p0,E0,un
    Type (Block_TYPE),pointer:: B
    TYPE (Mesh_TYPE),pointer:: MP
    Type (BC_MSG_TYPE),pointer:: Bc
    
    MP=> Mesh(nMesh)
    B => MP%Block(mBlock)
    NVAR1=MP%NVAR
	
!	Scheme=MP%Iflag_Scheme         ! 数值格式 （不同网格上采用不同格式）
    
	IFlux=MP%Iflag_Flux            ! 通量技术
    Reconstruction=MP%IFlag_Reconstruction  ! 重构方式
    nx1=B%nx ; ny1=B%ny ; nz1=B%nz


! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量
!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED (MP,B,NVAR1,IFlux,Reconstruction,nx1,ny1,nz1,gamma,Flux,d,uu,v,w,p)
!$OMP DO   
	do k=1,nz1
      do j=1,ny1
        do i=1,nx1
          do m=1,NVAR1
            Flux(m,i,j,k)=0.d0            ! 初始化
          enddo
        enddo
      enddo
    enddo
!$OMP END DO
 
!$OMP DO   
	do k=1,nz1-1 
      do j=1,ny1-1
        do i=1,nx1

!      设定边界格式 (区分内点格式及边界点格式)
		    Scheme=MP%Iflag_Scheme
		   
		   if(B%BcI(j,k,1)==1) then            ! 是否为物理边界
		     if(i .eq. 1) then
			   Scheme=Scheme_CD2    ! 第1个界面采用2阶中心格式 （利用虚网格，保持通量守恒）  
             else if (i < LAP) then
		       Scheme=MP%Bound_Scheme            ! 边界格式
             endif
           endif
		 
		   if(B%BcI(j,k,2)==1) then      ! 是否为物理边界
			if(i .eq. nx1 ) then
			 Scheme=Scheme_CD2        ! 最右边界，采用2阶中心
			else if (i > nx1-LAP+1) then
		     Scheme=MP%Bound_Scheme   ! 边界格式
            endif
           endif

          if(B%IF_OverLimit == 1)  Scheme=Scheme_UD1            ! 物理量超限（如出现负温度）， 该块采用1阶迎风
		 		 

!-----侧面的法方向和两个切方向-----------------------------------------------------------------------------
!  为了节省内存，本程序不存储界面的法方向及切方向，在使用时计算 (会增加些计算量)
          s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)            ! 归一化的法方向
          t1x=B%x(i,j+1,k)-B%x(i,j,k+1); t1y=B%y(i,j+1,k)-B%y(i,j,k+1); t1z=B%z(i,j+1,k)-B%z(i,j,k+1)   ! 对角线1
          t10=1.d0/(sqrt(t1x*t1x+t1y*t1y+t1z*t1z))  
          s2x=t1x*t10; s2y=t1y*t10; s2z=t1z*t10                             ! 归一化的切方向1 （对角线1）
          s3x=s1y*s2z-s1z*s2y; s3y=s1z*s2x-s1x*s2z ; s3z= s1x*s2y-s1y*s2x   ! 切方向2 （法方向叉乘切方向1）
!--------------------------------------------------------------------------------------------
! 内点，采用2,3阶格式重构  (I-1/2,J,K) 点的值；
! 理论手册及教科书中，通常写重构出(I+1/2,J,K)点的值。 本程序中重构(I-1/2,J,K)点的值为了方便 
! 需要使用4个点的值： I-2, I-1, I, I+1; 其中左值(UL)使用I-2,I-1,I点的值重构； 右值(UR)使用I-1,I,I+1点的值重构  
          if(Reconstruction .eq. Reconst_Original) then
! 使用原始变量重构 U0(:,:) 存储的是4个点上的密度、速度、压力
                U0(:,1)=d(i-LAP:i+LAP-1,j,k) 
			    U0(:,2)=uu(i-LAP:i+LAP-1,j,k) 
			    U0(:,3)=v(i-LAP:i+LAP-1,j,k)
			    U0(:,4)=w(i-LAP:i+LAP-1,j,k)
			    U0(:,5)=p(i-LAP:i+LAP-1,j,k)
            call Reconstuction_original(U0,UL,UR,gamma,Scheme)          ! 数值格式
          else if (Reconstruction .eq. Reconst_Conservative) then
! 使用守恒变量重构 U0(:,:) 存储的是4个点上的守恒变量（质量密度、动量密度和能量密度）
            do m=1,5
               U0(:,m)=B%U(m,i-LAP:i+LAP-1,j,k)
            enddo
            call Reconstuction_conservative(U0,UL,UR,gamma,Scheme)
          else
! 采用特征变量重构
            do m=1,5
              U0(:,m)=B%U(m,i-LAP:i+LAP-1,j,k)
            enddo
            call Reconstuction_Characteristic(U0,UL,UR,gamma,Scheme)
          endif   
!-------重构结束，得到(I-1/2,J,K)点上的左、右值UL,UR  (程序中为原始变量)------------
!----检查 UL, UR 中的密度、压力是否为负 -------------------------
! 如果为负，则使用1阶迎风        
	  if(IF_Scheme_Positivity == 1) then
		if(UL(1) .le. Lim_Zero  .or. UL(5) .le. Lim_Zero) then      ! Lim_Zero=1.d-20 
          UL(1)=d(i-1,j,k)            ! 1阶迎风
		  UL(2)=uu(i-1,j,k)
		  UL(3)=v(i-1,j,k)
		  UL(4)=w(i-1,j,k)
		  UL(5)=p(i-1,j,k)
		endif

       if(UR(1) .le. Lim_Zero .or. UR(5) .le. Lim_Zero) then
          UR(1)=d(i,j,k)         ! 1阶迎风
		  UR(2)=uu(i,j,k)
		  UR(3)=v(i,j,k)
		  UR(4)=w(i,j,k)
		  UR(5)=p(i,j,k)
       endif
     endif




!----------------------------------------------------------------
!-----将UL,UR进行坐标旋转变换，得到坐标系(s1,s2,s3) (法方向，切方向1，切方向2) 中的守恒变量QL,QR
!  标量（密度、压力）保持不变； 向量（速度）投影到新的坐标轴方向
!  左值
          QL(1)=UL(1)                           ! 密度（标量，坐标旋转保持不变）
          QL(2)=UL(2)*s1x+UL(3)*s1y+UL(4)*s1z   ! 法向速度
          QL(3)=UL(2)*s2x+UL(3)*s2y+UL(4)*s2z   ! 切方向1的速度分量
          QL(4)=UL(2)*s3x+UL(3)*s3y+UL(4)*s3z   ! 切方向2的速度分量
          QL(5)=UL(5)                           ! 压力（标量） 
! 右值 （变量含义与左值相同）
          QR(1)=UR(1)                            
          QR(2)=UR(2)*s1x+UR(3)*s1y+UR(4)*s1z    
          QR(3)=UR(2)*s2x+UR(3)*s2y+UR(4)*s2z    
          QR(4)=UR(2)*s3x+UR(3)*s3y+UR(4)*s3z    
          QR(5)=UR(5)                             
!---------通量分裂 (FVS或FDS), 根据左、右值计算穿过界面的通量(扩展1维问题）-------------
          if(IFlux .eq. Flux_Steger_Warming ) then
            call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)     ! Steger-Warming FVS方法 
          else  if(IFlux .eq. Flux_Roe ) then
            call Flux_Roe_1D(QL,QR,Flux0,gamma)                 ! Roe FDS方法
          else  if(IFlux .eq. Flux_Van_Leer ) then
            call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)           ! Van Leer FVS方法
          else  if(IFlux .eq. Flux_Ausm ) then          
            call Flux_Ausmpw_1Da(QL,QR,Flux0,gamma)             ! AUSM + 方法
          else
            call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,IFlux)      !HLL/HLLC FDS方法（近似Riemann解）
          endif
! ---将(s1,s2,s3)坐标系下的通量Flux0 进行变换，得到(x,y,z)坐标系下的通量； 标量保持不变，向量进行投影
! ---变换后，乘以面积，就得到穿过(I-1/2,J,K)点 （即 (i,j+1/2,k+1/2)点）所在侧面的通量	  
          Flux(1,i,j,k)=-Flux0(1)*B%Si(i,j,k)                        ! 质量通量 （标量，不随坐标旋转而变化）
          Flux(2,i,j,k)=-(Flux0(2)*s1x+Flux0(3)*s2x+Flux0(4)*s3x)*B%Si(i,j,k)  ! x方向的动量通量  （投影到x方向）
          Flux(3,i,j,k)=-(Flux0(2)*s1y+Flux0(3)*s2y+Flux0(4)*s3y)*B%Si(i,j,k)  ! y方向的动量通量
          Flux(4,i,j,k)=-(Flux0(2)*s1z+Flux0(3)*s2z+Flux0(4)*s3z)*B%Si(i,j,k)  ! z方向的动量通量
          Flux(5,i,j,k)=-Flux0(5)*B%Si(i,j,k)                        ! 能量通量 （标量）
        enddo
      enddo
    enddo
!$OMP END DO

    
!----------------Residual -------------------------------
!$OMP DO   
	do k=1,nz1-1
      do j=1,ny1-1
        do i=1,nx1-1
          do m=1,5
            B%Res(m,i,j,k)=Flux(m,i+1,j,k)-Flux(m,i,j,k)        
          enddo
        enddo
      enddo
    enddo
!$OMP END DO

!$OMP END PARALLEL

   end subroutine flux_inviscous_i
!--------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
! 穿过j方向界面((I,J-1/2,K)点所在界面)的无粘通量 
   subroutine flux_inviscous_j(nMesh,mBlock)
    Use Global_Var
    Use Flow_Var
    use mod_struct_flux, only: Flux_steger_warming_1Da, Flux_Roe_1D, &
                               Flux_Van_Leer_1Da, Flux_Ausmpw_1Da, Flux_HLL_HLLC_1D
    implicit none
    real(PRE_EC),dimension(5):: UL,UR,QL,QR,Flux0
    real(PRE_EC):: U0(1-LAP:LAP,5)
    integer:: mBlock,i,j,k,m,nx1,ny1,nz1,ksub,nMesh
    integer:: Scheme,IFlux,Reconstruction
    real(PRE_EC):: s1x,s1y,s1z,s2x,s2y,s2z,s3x,s3y,s3z,s0,t10,t1x,t1y,t1z  ! 法方向(s1)及切方向(s2,s3)
    real(PRE_EC):: d0,uu0,v0,w0,p0,E0,un
    Type (Block_TYPE),pointer:: B
    TYPE (Mesh_TYPE),pointer:: MP
    Type (BC_MSG_TYPE),pointer:: Bc
    MP=> Mesh(nMesh)
    B => MP%Block(mBlock)
    Scheme=MP%Iflag_Scheme                  ! 数值格式 （不同网格上采用不同格式）
    IFlux=MP%Iflag_Flux                     ! 通量技术
    Reconstruction=MP%IFlag_Reconstruction  ! 重构方式
    nx1=B%nx ; ny1=B%ny ; nz1=B%nz
    Flux=0.d0                               ! 初始化
   
!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(MP,B,IFlux,Reconstruction,nx1,ny1,nz1,gamma,Flux,d,uu,v,w,p)

!$OMP DO   
    do k=1,nz1-1 
      do j=1,ny1
		do i=1,nx1-1   
  
  
  !      设定边界格式 (区分内点格式及边界点格式)

		   Scheme=MP%Iflag_Scheme
		    if(B%BcJ(i,k,1)==1) then            ! 是否为物理边界
             if(j .eq. 1) then
			   Scheme=Scheme_CD2    ! 第1个界面采用2阶中心格式 （利用虚网格，保持通量守恒）  
			 else if( j < LAP ) then                  ! 靠近左边界
               Scheme=MP%Bound_Scheme            ! 边界格式
			 endif
            endif
		  
		    if(B%BcJ(i,k,2)==1) then
             if(j .eq. ny1) then
			   Scheme=Scheme_CD2    ! 第1个界面采用2阶中心格式 （利用虚网格，保持通量守恒）  
		     else if( j > ny1-LAP+1 ) then      ! 靠近右边界
               Scheme=MP%Bound_Scheme            ! 边界格式
			 endif
            endif

           if(B%IF_OverLimit == 1)  Scheme=Scheme_UD1            ! 物理量超限（如出现负温度）， 该块采用1阶迎风



          s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向
          t1x=B%x(i+1,j,k+1)-B%x(i,j,k); t1y=B%y(i+1,j,k+1)-B%y(i,j,k) ; t1z=B%z(i+1,j,k+1)-B%z(i,j,k)  ! 对角线1
          t10=1.d0/(sqrt(t1x*t1x+t1y*t1y+t1z*t1z))  
          s2x=t1x*t10; s2y=t1y*t10; s2z=t1z*t10     ! 归一化的切方向1 （对角线1）
          s3x=s1y*s2z-s1z*s2y; s3y=s1z*s2x-s1x*s2z ; s3z= s1x*s2y-s1y*s2x   ! 切方向2 （法方向叉乘切方向1: s1*s2）
!--------------------------------------------------------------------------------------------
! 内点，采用2,3阶格式重构  (I,J-1/2,K) 点的值；
          if(Reconstruction .eq. Reconst_Original) then
! 使用原始变量重构 U0(:,:) 存储的是4个点上的密度、速度、压力
             U0(:,1)=d(i,j-LAP:j+LAP-1,k)
			 U0(:,2)=uu(i,j-LAP:j+LAP-1,k) 
			 U0(:,3)=v(i,j-LAP:j+LAP-1,k)
			 U0(:,4)=w(i,j-LAP:j+LAP-1,k)
			 U0(:,5)=p(i,j-LAP:j+LAP-1,k)
            call Reconstuction_original(U0,UL,UR,gamma,Scheme)          ! 数值格式
          else if (Reconstruction .eq. Reconst_Conservative) then
! 使用守恒变量重构 U0(:,:) 存储的是4个点上的守恒变量（质量密度、动量密度和能量密度）
            do m=1,5
              U0(:,m)=B%U(m,i,j-LAP:j+LAP-1,k)
            enddo
            call Reconstuction_conservative(U0,UL,UR,gamma,Scheme)
          else
! 采用特征变量重构
            do m=1,5
              U0(:,m)=B%U(m,i,j-LAP:j+LAP-1,k)
            enddo
            call Reconstuction_Characteristic(U0,UL,UR,gamma,Scheme)
          endif
!-------重构结束，得到(I,J-1/2,K)点上的左、右值UL,UR  （原始变量：密度、速度、压力）------------

!----检查 UL, UR 中的密度、压力是否为负 -------------------------
! 如果为负，则使用1阶迎风        
	  if(IF_Scheme_Positivity == 1) then
		if(UL(1) .le. Lim_Zero  .or. UL(5) .le. Lim_Zero) then      ! Lim_Zero=1.d-20 
          UL(1)=d(i,j-1,k)            ! 1阶迎风
		  UL(2)=uu(i,j-1,k)
		  UL(3)=v(i,j-1,k)
		  UL(4)=w(i,j-1,k)
		  UL(5)=p(i,j-1,k)
		endif

       if(UR(1) .le. Lim_Zero .or. UR(5) .le. Lim_Zero) then
          UR(1)=d(i,j,k)         ! 1阶迎风
		  UR(2)=uu(i,j,k)
		  UR(3)=v(i,j,k)
		  UR(4)=w(i,j,k)
		  UR(5)=p(i,j,k)
       endif
     endif






!-----将UL,UR进行坐标旋转变换，得到坐标系(s1,s2,s3) (法方向，切方向1，切方向2) 中的守恒变量QL,QR
!  左值
          QL(1)=UL(1)                           ! 密度（标量，坐标旋转保持不变）
          QL(2)=UL(2)*s1x+UL(3)*s1y+UL(4)*s1z   ! 法方向的速度分量  （投影到法方向）
          QL(3)=UL(2)*s2x+UL(3)*s2y+UL(4)*s2z   ! 切方向1的速度分量
          QL(4)=UL(2)*s3x+UL(3)*s3y+UL(4)*s3z   ! 切方向2的速度分量
          QL(5)=UL(5)                           ! 压力（标量） 
! 右值 （各变量含义同上）
          QR(1)=UR(1)                           
          QR(2)=UR(2)*s1x+UR(3)*s1y+UR(4)*s1z   
          QR(3)=UR(2)*s2x+UR(3)*s2y+UR(4)*s2z   
          QR(4)=UR(2)*s3x+UR(3)*s3y+UR(4)*s3z   
          QR(5)=UR(5)                          
!---------通量分裂 (FVS或FDS), 根据左、右值计算穿过界面的通量(扩展1维问题）-------------
          if(IFlux .eq. Flux_Steger_Warming ) then
            call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)     ! Steger-Warming FVS方法 
          else  if(IFlux .eq. Flux_Roe ) then
            call Flux_Roe_1D(QL,QR,Flux0,gamma)                 ! Roe FDS方法
          else  if(IFlux .eq. Flux_Van_Leer ) then
            call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)           ! Van Leer FVS方法
          else  if(IFlux .eq. Flux_Ausm ) then
            call Flux_Ausmpw_1Da(QL,QR,Flux0,gamma)   
          else 
            call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,IFlux)      !HLL/HLLC FDS方法（近似Riemann解）
          endif
! ---将(s1,s2,s3)坐标系下的通量Flux0 进行变换，得到(x,y,z)坐标系下的通量； 标量保持不变，向量进行投影	  
          Flux(1,i,j,k)=-Flux0(1)*B%Sj(i,j,k)                                  ! 质量通量 （标量，不随坐标旋转而变化）
          Flux(2,i,j,k)=-(Flux0(2)*s1x+Flux0(3)*s2x+Flux0(4)*s3x)*B%Sj(i,j,k)  ! x方向的动量通量  （投影到x方向）
          Flux(3,i,j,k)=-(Flux0(2)*s1y+Flux0(3)*s2y+Flux0(4)*s3y)*B%Sj(i,j,k)  ! y方向的动量通量
          Flux(4,i,j,k)=-(Flux0(2)*s1z+Flux0(3)*s2z+Flux0(4)*s3z)*B%Sj(i,j,k)  ! z方向的动量通量
          Flux(5,i,j,k)=-Flux0(5)*B%Sj(i,j,k)                                  ! 能量通量 （标量） 
        enddo
      enddo
    enddo
!$OMP END DO   

!----------------Residual -------------------------------
!$OMP DO   
    do k=1,nz1-1
      do j=1,ny1-1
        do i=1,nx1-1
          do m=1,5
            B%Res(m,i,j,k)=B%Res(m,i,j,k)+Flux(m,i,j+1,k)-Flux(m,i,j,k)           
          enddo
        enddo
      enddo
    enddo
!$OMP END DO  
!$OMP END PARALLEL
 
   end subroutine flux_inviscous_j
!--------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
! 穿过k方向界面((I,J,K-1/2)点所在界面)的无粘通量 
   subroutine flux_inviscous_k(nMesh,mBlock)
    Use Global_Var
    Use Flow_Var
    use mod_struct_flux, only: Flux_steger_warming_1Da, Flux_Roe_1D, &
                               Flux_Van_Leer_1Da, Flux_Ausmpw_1Da, Flux_HLL_HLLC_1D
    implicit none
    real(PRE_EC),dimension(5):: UL,UR,QL,QR,Flux0
    real(PRE_EC):: U0(1-LAP:LAP,5)
    integer:: mBlock,i,j,k,m,nx1,ny1,nz1,ksub,nMesh
    integer:: Scheme,IFlux,Reconstruction
    real(PRE_EC):: s1x,s1y,s1z,s2x,s2y,s2z,s3x,s3y,s3z,s0,t10,t1x,t1y,t1z  ! 法方向(s1)及切方向(s2,s3)
    real(PRE_EC):: d0,uu0,v0,w0,p0,E0,un
    Type (Block_TYPE),pointer:: B
    TYPE (Mesh_TYPE),pointer:: MP
    Type (BC_MSG_TYPE),pointer:: Bc
    MP=> Mesh(nMesh)
    B => MP%Block(mBlock)
    Scheme=MP%Iflag_Scheme                  ! 数值格式 （不同网格上采用不同格式）
    IFlux=MP%Iflag_Flux                     ! 通量技术
    Reconstruction=MP%IFlag_Reconstruction  ! 重构方式
    nx1=B%nx ; ny1=B%ny ; nz1=B%nz
    
	Flux=0.d0                               ! 初始化

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(MP,B,IFlux,Reconstruction,nx1,ny1,nz1,gamma,Flux,d,uu,v,w,p)

!$OMP  DO  
    do k=1,nz1 
      do j=1,ny1-1
        do i=1,nx1-1

 !      设定边界格式 (区分内点格式及边界点格式)

		   Scheme=MP%Iflag_Scheme
		   if(B%BcK(i,j,1)==1) then            ! 是否为物理边界
              if(k .eq. 1) then
 			   Scheme=Scheme_CD2    ! 第1个界面采用2阶中心格式 （利用虚网格，保持通量守恒）  
			  else if( k < LAP ) then                  ! 靠近左边界
               Scheme=MP%Bound_Scheme            ! 边界格式
			  endif
           endif
	
		   if(B%BcK(i,j,2)==1) then
             if(k .eq. nz1) then
			   Scheme=Scheme_CD2    ! 第1个界面采用2阶中心格式 （利用虚网格，保持通量守恒）  
		     else if( k > nz1-LAP+1 ) then      ! 靠近右边界
               Scheme=MP%Bound_Scheme            ! 边界格式
			 endif
           endif
           
		   if(B%IF_OverLimit == 1)  Scheme=Scheme_UD1            ! 物理量超限（如出现负温度）， 该块采用1阶迎风


          s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向
          t1x=B%x(i+1,j+1,k)-B%x(i,j,k); t1y=B%y(i+1,j+1,k)-B%y(i,j,k) ; t1z=B%z(i+1,j+1,k)-B%z(i,j,k)  ! 对角线1
          t10=1.d0/(sqrt(t1x*t1x+t1y*t1y+t1z*t1z))  
          s2x=t1x*t10; s2y=t1y*t10; s2z=t1z*t10     ! 归一化的切方向1 （对角线1）
          s3x=s1y*s2z-s1z*s2y; s3y=s1z*s2x-s1x*s2z ; s3z= s1x*s2y-s1y*s2x   ! 切方向2 （法方向叉乘切方向1: s1*s2）
!--------------------------------------------------------------------------------------------
! 内点，采用2,3阶格式重构  (I,J-1/2,K) 点的值；
          if(Reconstruction .eq. Reconst_Original) then
! 使用原始变量重构 U0(:,:) 存储的是4个点上的密度、速度、压力
              U0(:,1)=d(i,j,k-LAP:k+LAP-1) 
			  U0(:,2)=uu(i,j,k-LAP:k+LAP-1) 
			  U0(:,3)=v(i,j,k-LAP:k+LAP-1)
			  U0(:,4)=w(i,j,k-LAP:k+LAP-1)
			  U0(:,5)=p(i,j,k-LAP:k+LAP-1)
            call Reconstuction_original(U0,UL,UR,gamma,Scheme)          ! 数值格式
          else if (Reconstruction .eq. Reconst_Conservative) then
! 使用守恒变量重构 U0(:,:) 存储的是4个点上的守恒变量（质量密度、动量密度和能量密度）
            do m=1,5
              U0(:,m)=B%U(m,i,j,k-LAP:k+LAP-1)
            enddo
            call Reconstuction_conservative(U0,UL,UR,gamma,Scheme)
          else
! 采用特征变量重构
            do m=1,5
              U0(:,m)=B%U(m,i,j,k-LAP:k+LAP-1)
            enddo
            call Reconstuction_Characteristic(U0,UL,UR,gamma,Scheme)
          endif  
!-------重构结束，得到(I,J-1/2,K)点上的左、右值UL,UR  （原始变量：密度、速度、压力）------------

!----检查 UL, UR 中的密度、压力是否为负 -------------------------
! 如果为负，则使用1阶迎风        
	  if(IF_Scheme_Positivity == 1) then
		if(UL(1) .le. Lim_Zero  .or. UL(5) .le. Lim_Zero) then      ! Lim_Zero=1.d-20 
          UL(1)=d(i,j,k-1)            ! 1阶迎风
		  UL(2)=uu(i,j,k-1)
		  UL(3)=v(i,j,k-1)
		  UL(4)=w(i,j,k-1)
		  UL(5)=p(i,j,k-1)
		endif

       if(UR(1) .le. Lim_Zero .or. UR(5) .le. Lim_Zero) then
          UR(1)=d(i,j,k)         ! 1阶迎风
		  UR(2)=uu(i,j,k)
		  UR(3)=v(i,j,k)
		  UR(4)=w(i,j,k)
		  UR(5)=p(i,j,k)
       endif
     endif



!-----将UL,UR进行坐标旋转变换，得到坐标系(s1,s2,s3) (法方向，切方向1，切方向2) 中的守恒变量QL,QR
!  左值
          QL(1)=UL(1)                           ! 密度（标量，坐标旋转保持不变）
          QL(2)=UL(2)*s1x+UL(3)*s1y+UL(4)*s1z   ! 法方向的速度分量  （投影到法方向）
          QL(3)=UL(2)*s2x+UL(3)*s2y+UL(4)*s2z   ! 切方向1的速度分量
          QL(4)=UL(2)*s3x+UL(3)*s3y+UL(4)*s3z   ! 切方向2的速度分量
          QL(5)=UL(5)                           ! 压力（标量） 
! 右值 （各变量含义同上）
          QR(1)=UR(1)                           
          QR(2)=UR(2)*s1x+UR(3)*s1y+UR(4)*s1z   
          QR(3)=UR(2)*s2x+UR(3)*s2y+UR(4)*s2z   
          QR(4)=UR(2)*s3x+UR(3)*s3y+UR(4)*s3z   
          QR(5)=UR(5)                          
!---------通量分裂 (FVS或FDS), 根据左、右值计算穿过界面的通量(扩展1维问题）-------------
          if(IFlux .eq. Flux_Steger_Warming ) then
            call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)     ! Steger-Warming FVS方法 
          else  if(IFlux .eq. Flux_Roe ) then
            call Flux_Roe_1D(QL,QR,Flux0,gamma)                ! Roe FDS方法
          else  if(IFlux .eq. Flux_Van_Leer ) then
            call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)          ! Van Leer FVS方法
          else  if(IFlux .eq. Flux_Ausm ) then
            call Flux_Ausmpw_1Da(QL,QR,Flux0,gamma)   
          else
            call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,IFlux)    !HLL/HLLC FDS方法（近似Riemann解）
          endif
! ---将(s1,s2,s3)坐标系下的通量Flux0 进行变换，得到(x,y,z)坐标系下的通量； 标量保持不变，向量进行投影
	  
          Flux(1,i,j,k)=-Flux0(1)*B%Sk(i,j,k)                                  ! 质量通量 （标量，不随坐标旋转而变化）
          Flux(2,i,j,k)=-(Flux0(2)*s1x+Flux0(3)*s2x+Flux0(4)*s3x)*B%Sk(i,j,k)  ! x方向的动量通量  （投影到x方向）
          Flux(3,i,j,k)=-(Flux0(2)*s1y+Flux0(3)*s2y+Flux0(4)*s3y)*B%Sk(i,j,k)  ! y方向的动量通量
          Flux(4,i,j,k)=-(Flux0(2)*s1z+Flux0(3)*s2z+Flux0(4)*s3z)*B%Sk(i,j,k)  ! z方向的动量通量
          Flux(5,i,j,k)=-Flux0(5)*B%Sk(i,j,k)                                  ! 能量通量 （标量）
        enddo
      enddo
    enddo
!$OMP END DO  

!---------------------------------------------------------------
!  边界处理
!$OMP  DO  
    do  ksub=1,B%subface
      Bc=> B%bc_msg(ksub)
      if(Bc%bc .gt. 0  .and. (Bc%face .eq. 3 .or. Bc%face .eq. 6)) then   ! 非内边界，且为k-或k+边界
        k=Bc%kb
        do j= Bc%jb, Bc%je-1
          do i= Bc%ib, Bc%ie-1
            d0=(d(i,j,k)+d(i,j,k-1))*0.5d0 
            uu0=(uu(i,j,k)+uu(i,j,k-1))*0.5d0 
            v0=(v(i,j,k)+v(i,j,k-1))*0.5d0 
            w0=(w(i,j,k)+w(i,j,k-1))*0.5d0
            p0=(p(i,j,k)+p(i,j,k-1))*0.5d0
            E0=p0/(gamma-1.d0)+d0*(uu0*uu0+v0*v0+w0*w0)*0.5d0
            un=uu0*B%nk1(i,j,k)+v0*B%nk2(i,j,k)+w0*B%nk3(i,j,k)          ! 法向速度       
            Flux(1,i,j,k)=-d0*un*B%Sk(i,j,k)
            Flux(2,i,j,k)=-(d0*un*uu0+p0*B%nk1(i,j,k))*B%Sk(i,j,k)
            Flux(3,i,j,k)=-(d0*un*v0+p0*B%nk2(i,j,k))*B%Sk(i,j,k)
            Flux(4,i,j,k)=-(d0*un*w0+p0*B%nk3(i,j,k))*B%Sk(i,j,k)
            Flux(5,i,j,k)=-(E0+p0)*un*B%Sk(i,j,k)
          enddo
        enddo
      endif
    enddo
!$OMP END DO  

!----------------Residual -------------------------------
!$OMP  DO  
    do k=1,nz1-1
      do j=1,ny1-1
        do i=1,nx1-1
          do m=1,5
            B%Res(m,i,j,k)=B%Res(m,i,j,k)+Flux(m,i,j,k+1)-Flux(m,i,j,k)           
          enddo
        enddo
      enddo
    enddo
!$OMP END DO  
!$OMP END PARALLEL

   end subroutine flux_inviscous_k

!----------------粘性通量---------------------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------------------------------------------
!-----------i-方向的粘性通量-----------------------------------------------------------------  
! 为了保证稳定性，避免使用计算域角部的网格点（如 (0,0,k), (0,j,0) 等网格点处的值）

   subroutine flux_viscous_i(nMesh,mBlock)
    Use Global_Var
    Use Flow_Var
    implicit none
    integer:: mBlock,i,j,k,m,nx,ny,nz,nMesh
    real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz,s1x,s1y,s1z,ui,vi,wi,Ti,uj,vj,wj,Tj,uk,vk,wk,Tk
    real(PRE_EC):: ux,uy,uz,vx,vy,vz,wx,wy,wz,Tx,Ty,Tz
    real(PRE_EC):: s11,s12,s13,s22,s23,s33,u1,v1,w1,E1,E2,E3
    real(PRE_EC):: mu0,k0
    real(PRE_EC):: tmp1,tmp2
    real(PRE_EC):: ui1,ui2,uj1,uj2,uk1,uk2,vi1,vi2,vj1,vj2,vk1,vk2, &
	               wi1,wi2,wj1,wj2,wk1,wk2,Ti1,Ti2,Tj1,Tj2,Tk1,Tk2

    Type (Block_TYPE),pointer:: B
    TYPE (Mesh_TYPE),pointer:: MP
    MP=> Mesh(nMesh)
    B => MP%Block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
    tmp1=4.d0/3.d0; tmp2=2.d0/3.d0


! revised 2017-7-11
if(B%IF_OverLimit == 1 ) then     ! 物理量超限， 采用1阶迎风计算无粘性， 不计算粘性项

!$OMP  PARALLEL DO  
   do k=1,nz-1 
   do j=1,ny-1
     B%Surf1(j,k,1)=0.d0
     B%Surf1(j,k,2)=0.d0
     B%Surf1(j,k,3)=0.d0
     B%Surf4(j,k,1)=0.d0
     B%Surf4(j,k,2)=0.d0
     B%Surf4(j,k,3)=0.d0
   enddo
   enddo
!$OMP END  PARALLEL DO  

	return             ! 跳过粘性项计算
endif

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(MP,B,nx,ny,nz,tmp1,tmp2,Cp,PrL,PrT,uu,v,w,T,flux)

!$OMP  DO  
    do k=1,nz-1 
      do j=1,ny-1
        do i=1,nx
          s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)  ! 归一化的法方向
! 物理量对于计算坐标（下标）的导数;  (I-1/2,J,K)点
! 注意，物理量存储在网格中心，与坐标（存储在网格节点）有区别
! 导数均采用二阶中心差分计算，并利用周围网格点导数的平均
! e.g. uj(I-1/2,J,K)=0.5*(uj(I,J,K)+uj(I-1,J,K));  uj(I,J,K)=(uu(I,J+1,k)-uu(I,J-1,K))*0.5
          ui=uu(i,j,k)-uu(i-1,j,k)              ! du/di (I-1/2,J,K)= u(I,J,K)-u(I-1,J,K)
          vi=v(i,j,k)-v(i-1,j,k)
          wi=w(i,j,k)-w(i-1,j,k)
          Ti=T(i,j,k)-T(i-1,j,k)
  
  ! 为了保证稳定性，避免使用计算域角部的网格点（如 (0,0,k), (0,j,0) 等网格点处的值）
      
		  if( (i==1 .or. i==B%nx) .and. (j==1 .or. j==B%ny-1) ) then
           
		   uj1=0.d0; vj1=0.d0; wj1=0.d0; Tj1=0.d0
		   uj2=0.d0; vj2=0.d0; wj2=0.d0; Tj2=0.d0
		   		 
		  else 
		   uj1=0.5d0*(uu(i,j-1,k)+uu(i-1,j-1,k))
		   vj1=0.5d0*(v(i,j-1,k)+v(i-1,j-1,k))
		   wj1=0.5d0*(w(i,j-1,k)+w(i-1,j-1,k))
		   Tj1=0.5d0*(T(i,j-1,k)+T(i-1,j-1,k))
		   uj2=0.5d0*(uu(i,j+1,k)+uu(i-1,j+1,k))
		   vj2=0.5d0*(v(i,j+1,k)+v(i-1,j+1,k))
		   wj2=0.5d0*(w(i,j+1,k)+w(i-1,j+1,k))
		   Tj2=0.5d0*(T(i,j+1,k)+T(i-1,j+1,k))
          endif
        
          if( (i==1 .or. i==B%nx) .and. (k==1 .or. k==B%nz-1) ) then

 		   uk1=0.d0; vk1=0.d0; wk1=0.d0; Tk1=0.d0
		   uk2=0.d0; vk2=0.d0; wk2=0.d0; Tk2=0.d0
		 
		  else
            uk1=0.5d0*(uu(i,j,k-1)+uu(i-1,j,k-1))
            vk1=0.5d0*(v(i,j,k-1)+v(i-1,j,k-1))
            wk1=0.5d0*(w(i,j,k-1)+w(i-1,j,k-1))
            Tk1=0.5d0*(T(i,j,k-1)+T(i-1,j,k-1))
            uk2=0.5d0*(uu(i,j,k+1)+uu(i-1,j,k+1))
            vk2=0.5d0*(v(i,j,k+1)+v(i-1,j,k+1))
            wk2=0.5d0*(w(i,j,k+1)+w(i-1,j,k+1))
            Tk2=0.5d0*(T(i,j,k+1)+T(i-1,j,k+1))
          endif
 		   uj=0.5d0*(uj2-uj1)
 		   vj=0.5d0*(vj2-vj1)
 		   wj=0.5d0*(wj2-wj1)
 		   Tj=0.5d0*(Tj2-Tj1)
 		   uk=0.5d0*(uk2-uk1)
		   vk=0.5d0*(vk2-vk1)
		   wk=0.5d0*(wk2-wk1)
		   Tk=0.5d0*(Tk2-Tk1)

          ix=B%ix1(i,j,k); iy=B%iy1(i,j,k); iz=B%iz1(i,j,k)
          jx=B%jx1(i,j,k); jy=B%jy1(i,j,k); jz=B%jz1(i,j,k)
          kx=B%kx1(i,j,k); ky=B%ky1(i,j,k); kz=B%kz1(i,j,k)

!----对物理坐标的偏导数----------------------------------------------
          ux=ui*ix+uj*jx+uk*kx
          vx=vi*ix+vj*jx+vk*kx
          wx=wi*ix+wj*jx+wk*kx
          Tx=Ti*ix+Tj*jx+Tk*kx

          uy=ui*iy+uj*jy+uk*ky
          vy=vi*iy+vj*jy+vk*ky
          wy=wi*iy+wj*jy+wk*ky
          Ty=Ti*iy+Tj*jy+Tk*ky

          uz=ui*iz+uj*jz+uk*kz
          vz=vi*iz+vj*jz+vk*kz
          wz=wi*iz+wj*jz+wk*kz
          Tz=Ti*iz+Tj*jz+Tk*kz

!---粘性应力----------------------------------------------------------
    ! (I-1/2,J,k)点, 即 (i,j+1/2,k+1/2)点 处的粘性系数
         mu0=0.5d0*(B%mu(i,j,k)+B%mu_t(i,j,k) + B%mu(i-1,j,k)+B%mu_t(i-1,j,k))  
         k0=0.5d0*Cp*(B%mu(i,j,k)/PrL + B%mu_t(i,j,k)/PrT + B%mu(i-1,j,k)/PrL +B%mu_t(i-1,j,k)/PrT)   ! (I-1/2,J,K) 点的热传导系数
!-----------------------------------------------------------
         s11=mu0*(tmp1*ux-tmp2*(vy+wz))    ! tmp1=4.d0/3.d0; tmp2=2.d0/3.d0
         s12=mu0*(uy+vx)
         s13=mu0*(uz+wx)
         s22=mu0*(tmp1*vy-tmp2*(ux+wz))
         s23=mu0*(vz+wy)
         s33=mu0*(tmp1*wz-tmp2*(ux+vy))
!---(I-1/2,J,K)点上的速度--------------------------------------------   
         u1=(uu(i,j,k)+uu(i-1,j,k))*0.5d0
         v1=(v(i,j,k)+v(i-1,j,k))*0.5d0
         w1=(w(i,j,k)+w(i-1,j,k))*0.5d0
!--能量的粘性通量-------------------------------
         E1=u1*s11+v1*s12+w1*s13+k0*Tx
         E2=u1*s12+v1*s22+w1*s23+k0*Ty
         E3=u1*s13+v1*s23+w1*s33+k0*Tz
!-------通量=无粘通量+粘性通量---------------------------------------
! 粘性通量=Fv.n=Fv1*s1x+Fv2*s1y+Fv3*s1z---------------
         Flux(2,i,j,k)= (s11*s1x+s12*s1y+s13*s1z)* B%Si(i,j,k)  
         Flux(3,i,j,k)= (s12*s1x+s22*s1y+s23*s1z)* B%Si(i,j,k)
         Flux(4,i,j,k)= (s13*s1x+s23*s1y+s33*s1z)* B%Si(i,j,k)
         Flux(5,i,j,k)= (E1*s1x+E2*s1y+E3*s1z)* B%Si(i,j,k)
!----------------------------------------------------------
      enddo
    enddo
  enddo
!$OMP END DO  

!---------------------------------------------------------
!$OMP  DO  

     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
     do m=2,5
      B%Res(m,i,j,k)=B%Res(m,i,j,k)+Flux(m,i+1,j,k)-Flux(m,i,j,k)           
     enddo
     enddo
     enddo
     enddo
!$OMP END DO  

!---记录边界处的粘性力-------------------------------------
!$OMP  DO  
   do k=1,nz-1 
   do j=1,ny-1
     B%Surf1(j,k,1)=Flux(2,1,j,k)
     B%Surf1(j,k,2)=Flux(3,1,j,k)
     B%Surf1(j,k,3)=Flux(4,1,j,k)
     B%Surf4(j,k,1)=Flux(2,nx,j,k)
     B%Surf4(j,k,2)=Flux(3,nx,j,k)
     B%Surf4(j,k,3)=Flux(4,nx,j,k)
   enddo
   enddo
!$OMP END DO  
!$OMP  END PARALLEL

 end subroutine flux_viscous_i

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
! j方向的粘性通量   
   subroutine flux_viscous_j(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,i,j,k,m,nx,ny,nz,nMesh
   real(PRE_EC):: s1x,s1y,s1z
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz,ui,vi,wi,Ti,uj,vj,wj,Tj,uk,vk,wk,Tk
   real(PRE_EC):: ux,uy,uz,vx,vy,vz,wx,wy,wz,Tx,Ty,Tz
   real(PRE_EC):: s11,s12,s13,s22,s23,s33,u1,v1,w1,E1,E2,E3
   real(PRE_EC):: mu0,k0
   real(PRE_EC):: tmp1,tmp2
   Type (Block_TYPE),pointer:: B
   TYPE (Mesh_TYPE),pointer:: MP
   real(PRE_EC):: ui1,ui2,uj1,uj2,uk1,uk2,vi1,vi2,vj1,vj2,vk1,vk2, &
	               wi1,wi2,wj1,wj2,wk1,wk2,Ti1,Ti2,Tj1,Tj2,Tk1,Tk2
   
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz

   tmp1=4.d0/3.d0; tmp2=2.d0/3.d0

  if(B%IF_OverLimit == 1 ) then     ! 物理量超限， 采用1阶迎风计算无粘性， 不计算粘性项
 !$OMP Parallel DO
   do k=1,nz-1 
   do i=1,nx-1
     B%Surf2(i,k,1)=0.d0
     B%Surf2(i,k,2)=0.d0
     B%Surf2(i,k,3)=0.d0
     B%Surf5(i,k,1)=0.d0
     B%Surf5(i,k,2)=0.d0
     B%Surf5(i,k,3)=0.d0
   enddo
   enddo
!$OMP END Parallel DO
   return
   endif


!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(MP,B,nx,ny,nz,tmp1,tmp2,Cp,PrL,PrT,uu,v,w,T,flux)
!$OMP DO
   do k=1,nz-1 
   do j=1,ny
   do i=1,nx-1
    s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向

! 物理量对于计算坐标（下标）的导数;  (I,J-1/2,K)点
   uj=uu(i,j,k)-uu(i,j-1,k)
   vj=v(i,j,k)-v(i,j-1,k)
   wj=w(i,j,k)-w(i,j-1,k)
   Tj=T(i,j,k)-T(i,j-1,k)


	   if( (j==1 .or. j==B%ny) .and. (i==1 .or. i==B%nx-1) ) then
 
 	    ui1=0.d0; vi1=0.d0; wi1=0.d0; Ti1=0.d0
	    ui2=0.d0; vi2=0.d0; wi2=0.d0; Ti2=0.d0

	   else
		ui1=0.5d0*(uu(i-1,j,k)+uu(i-1,j-1,k))
		vi1=0.5d0*(v(i-1,j,k)+v(i-1,j-1,k))
		wi1=0.5d0*(w(i-1,j,k)+w(i-1,j-1,k))
		Ti1=0.5d0*(T(i-1,j,k)+T(i-1,j-1,k))
		ui2=0.5d0*(uu(i+1,j,k)+uu(i+1,j-1,k))
		vi2=0.5d0*(v(i+1,j,k)+v(i+1,j-1,k))
		wi2=0.5d0*(w(i+1,j,k)+w(i+1,j-1,k))
		Ti2=0.5d0*(T(i+1,j,k)+T(i+1,j-1,k))
       endif
	
	   if( (j==1 .or. j==B%ny) .and. (k==1 .or. k==B%nz-1) ) then
   		   uk1=0.d0; vk1=0.d0; wk1=0.d0; Tk1=0.d0
		   uk2=0.d0; vk2=0.d0; wk2=0.d0; Tk2=0.d0

	   else	 
		 uk1=0.5d0*(uu(i,j,k-1)+uu(i,j-1,k-1))
		 vk1=0.5d0*(v(i,j,k-1)+v(i,j-1,k-1))
		 wk1=0.5d0*(w(i,j,k-1)+w(i,j-1,k-1))
		 Tk1=0.5d0*(T(i,j,k-1)+T(i,j-1,k-1))
		 uk2=0.5d0*(uu(i,j,k+1)+uu(i,j-1,k+1))
		 vk2=0.5d0*(v(i,j,k+1)+v(i,j-1,k+1))
		 wk2=0.5d0*(w(i,j,k+1)+w(i,j-1,k+1))
		 Tk2=0.5d0*(T(i,j,k+1)+T(i,j-1,k+1))
       endif
	    ui=0.5d0*(ui2-ui1)
	    vi=0.5d0*(vi2-vi1)
	    wi=0.5d0*(wi2-wi1)
	    Ti=0.5d0*(Ti2-Ti1)
	    uk=0.5d0*(uk2-uk1)
	    vk=0.5d0*(vk2-vk1)
	    wk=0.5d0*(wk2-wk1)
	    Tk=0.5d0*(Tk2-Tk1)


      ix=B%ix2(i,j,k); iy=B%iy2(i,j,k); iz=B%iz2(i,j,k)
      jx=B%jx2(i,j,k); jy=B%jy2(i,j,k); jz=B%jz2(i,j,k)
      kx=B%kx2(i,j,k); ky=B%ky2(i,j,k); kz=B%kz2(i,j,k)

!----对物理坐标的偏导数----------------------------------------------
   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
   Tx=Ti*ix+Tj*jx+Tk*kx

   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky
   Ty=Ti*iy+Tj*jy+Tk*ky

   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz
   Tz=Ti*iz+Tj*jz+Tk*kz

!---粘性应力----------------------------------------------------------
   mu0=0.5d0*(B%mu(i,j,k)+B%mu_t(i,j,k) + B%mu(i,j-1,k)+B%mu_t(i,j-1,k))  
   k0=0.5d0*Cp*(B%mu(i,j,k)/PrL + B%mu_t(i,j,k)/PrT + B%mu(i,j-1,k)/PrL +B%mu_t(i,j-1,k)/PrT)   ! (I-1/2,J,K) 点的热传导系数
   
   s11=mu0*(tmp1*ux-tmp2*(vy+wz))    ! tmp1=4.d0/3.d0; tmp2=2.d0/3.d0
   s12=mu0*(uy+vx)
   s13=mu0*(uz+wx)
   s22=mu0*(tmp1*vy-tmp2*(ux+wz))
   s23=mu0*(vz+wy)
   s33=mu0*(tmp1*wz-tmp2*(ux+vy))
!---(I,J-1/2,K)点上的速度--------------------------------------------   
   u1=(uu(i,j,k)+uu(i,j-1,k))*0.5d0
   v1=(v(i,j,k)+v(i,j-1,k))*0.5d0
   w1=(w(i,j,k)+w(i,j-1,k))*0.5d0
!--能量的粘性通量-------------------------------
   E1=u1*s11+v1*s12+w1*s13+k0*Tx
   E2=u1*s12+v1*s22+w1*s23+k0*Ty
   E3=u1*s13+v1*s23+w1*s33+k0*Tz

!-------通量=无粘通量+粘性通量---------------------------------------
! 粘性通量=Fv.n=Fv1*s1x+Fv2*s1y+Fv3*s1z---------------
    Flux(2,i,j,k)=(s11*s1x+s12*s1y+s13*s1z)* B%Sj(i,j,k)  
    Flux(3,i,j,k)=(s12*s1x+s22*s1y+s23*s1z)* B%Sj(i,j,k)
    Flux(4,i,j,k)=(s13*s1x+s23*s1y+s33*s1z)* B%Sj(i,j,k)
    Flux(5,i,j,k)=(E1*s1x+E2*s1y+E3*s1z)* B%Sj(i,j,k)

   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
     do m=2,5
      B%Res(m,i,j,k)=B%Res(m,i,j,k)+Flux(m,i,j+1,k)-Flux(m,i,j,k)           
     enddo
     enddo
     enddo
     enddo
!$OMP END DO

!---记录边界处的粘性力-------------------------------------

!$OMP  DO
   do k=1,nz-1 
   do i=1,nx-1
     B%Surf2(i,k,1)=Flux(2,i,1,k)
     B%Surf2(i,k,2)=Flux(3,i,1,k)
     B%Surf2(i,k,3)=Flux(4,i,1,k)
     B%Surf5(i,k,1)=Flux(2,i,ny,k)
     B%Surf5(i,k,2)=Flux(3,i,ny,k)
     B%Surf5(i,k,3)=Flux(4,i,ny,k)
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


 end subroutine flux_viscous_j
 
!------------------------------------------------------------------------------------------------------------- 
!------------------------------------------------------------------------------------------------------------- 
! k方向的粘性通量   
   subroutine flux_viscous_k(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,i,j,k,m,nx,ny,nz,nMesh
   real(PRE_EC):: s1x,s1y,s1z
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz,ui,vi,wi,Ti,uj,vj,wj,Tj,uk,vk,wk,Tk
   real(PRE_EC):: ux,uy,uz,vx,vy,vz,wx,wy,wz,Tx,Ty,Tz
   real(PRE_EC):: s11,s12,s13,s22,s23,s33,u1,v1,w1,E1,E2,E3
   real(PRE_EC):: mu0,k0
   real(PRE_EC):: tmp1,tmp2
   Type (Block_TYPE),pointer:: B
   TYPE (Mesh_TYPE),pointer:: MP
   real(PRE_EC):: ui1,ui2,uj1,uj2,uk1,uk2,vi1,vi2,vj1,vj2,vk1,vk2, &
	               wi1,wi2,wj1,wj2,wk1,wk2,Ti1,Ti2,Tj1,Tj2,Tk1,Tk2

   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz

   tmp1=4.d0/3.d0; tmp2=2.d0/3.d0

 if(B%IF_OverLimit == 1 ) then     ! 物理量超限， 采用1阶迎风计算无粘性， 不计算粘性项
!$OMP  Parallel DO
   do j=1,ny-1 
   do i=1,nx-1
     B%Surf3(i,j,1)=0.d0
     B%Surf3(i,j,2)=0.d0
     B%Surf3(i,j,3)=0.d0
     B%Surf6(i,j,1)=0.d0
     B%Surf6(i,j,2)=0.d0
     B%Surf6(i,j,3)=0.d0
   enddo
   enddo
!$OMP END Parallel DO
 return
  endif



!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(MP,B,nx,ny,nz,tmp1,tmp2,Cp,PrL,PrT,uu,v,w,T,flux)

!$OMP  DO
   do k=1,nz 
   do j=1,ny-1
   do i=1,nx-1
    s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向

! 物理量对于计算坐标（下标）的导数;  (I,J,K-1/2)点

   uk=uu(i,j,k)-uu(i,j,k-1)
   vk=v(i,j,k)-v(i,j,k-1)
   wk=w(i,j,k)-w(i,j,k-1)
   Tk=T(i,j,k)-T(i,j,k-1)

 	 if( (k==1 .or. k==B%nz) .and. (i==1 .or. i==B%nx-1) ) then
 	    ui1=0.d0; vi1=0.d0; wi1=0.d0; Ti1=0.d0
	    ui2=0.d0; vi2=0.d0; wi2=0.d0; Ti2=0.d0
    
	 else
	  ui1=0.5d0*(uu(i-1,j,k)+uu(i-1,j,k-1))
	  vi1=0.5d0*(v(i-1,j,k)+v(i-1,j,k-1))
	  wi1=0.5d0*(w(i-1,j,k)+w(i-1,j,k-1))
	  Ti1=0.5d0*(T(i-1,j,k)+T(i-1,j,k-1))
	  ui2=0.5d0*(uu(i+1,j,k)+uu(i+1,j,k-1))
	  vi2=0.5d0*(v(i+1,j,k)+v(i+1,j,k-1))
	  wi2=0.5d0*(w(i+1,j,k)+w(i+1,j,k-1))
	  Ti2=0.5d0*(T(i+1,j,k)+T(i+1,j,k-1))
     endif

  	 if( (k==1 .or. k==B%nz) .and. (j==1 .or. j==B%ny-1) ) then
	   uj1=0.d0; vj1=0.d0; wj1=0.d0; Tj1=0.d0
	   uj2=0.d0; vj2=0.d0; wj2=0.d0; Tj2=0.d0

	 else
      uj1=0.5d0*(uu(i,j-1,k)+uu(i,j-1,k-1))
      vj1=0.5d0*(v(i,j-1,k)+v(i,j-1,k-1))
      wj1=0.5d0*(w(i,j-1,k)+w(i,j-1,k-1))
      Tj1=0.5d0*(T(i,j-1,k)+T(i,j-1,k-1))
      uj2=0.5d0*(uu(i,j+1,k)+uu(i,j+1,k-1))
      vj2=0.5d0*(v(i,j+1,k)+v(i,j+1,k-1))
      wj2=0.5d0*(w(i,j+1,k)+w(i,j+1,k-1))
      Tj2=0.5d0*(T(i,j+1,k)+T(i,j+1,k-1))
	 endif
	  ui=0.5d0*(ui2-ui1)
	  vi=0.5d0*(vi2-vi1)
	  wi=0.5d0*(wi2-wi1)
	  Ti=0.5d0*(Ti2-Ti1)

	  uj=0.5d0*(uj2-uj1)
	  vj=0.5d0*(vj2-vj1)
	  wj=0.5d0*(wj2-wj1)
	  Tj=0.5d0*(Tj2-Tj1)


      ix=B%ix3(i,j,k); iy=B%iy3(i,j,k); iz=B%iz3(i,j,k)
      jx=B%jx3(i,j,k); jy=B%jy3(i,j,k); jz=B%jz3(i,j,k)
      kx=B%kx3(i,j,k); ky=B%ky3(i,j,k); kz=B%kz3(i,j,k)

!----对物理坐标的偏导数----------------------------------------------
   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
   Tx=Ti*ix+Tj*jx+Tk*kx

   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky
   Ty=Ti*iy+Tj*jy+Tk*ky

   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz
   Tz=Ti*iz+Tj*jz+Tk*kz

!---粘性应力----------------------------------------------------------
   mu0=0.5d0*(B%mu(i,j,k)+B%mu_t(i,j,k) + B%mu(i,j,k-1)+B%mu_t(i,j,k-1))          ! (I,J,k-1/2)点的粘性系数 （层流+湍流）  
   k0=0.5d0*Cp*(B%mu(i,j,k)/PrL + B%mu_t(i,j,k)/PrT + B%mu(i,j,k-1)/PrL +B%mu_t(i,j,k-1)/PrT)   ! (I,J,K-1/2) 点的热传导系数

   s11=mu0*(tmp1*ux-tmp2*(vy+wz))    ! tmp1=4.d0/3.d0; tmp2=2.d0/3.d0
   s12=mu0*(uy+vx)
   s13=mu0*(uz+wx)
   s22=mu0*(tmp1*vy-tmp2*(ux+wz))
   s23=mu0*(vz+wy)
   s33=mu0*(tmp1*wz-tmp2*(ux+vy))

!---(I,J,K-1/2)点上的速度--------------------------------------------   
   u1=(uu(i,j,k)+uu(i,j,k-1))*0.5d0
   v1=(v(i,j,k)+v(i,j,k-1))*0.5d0
   w1=(w(i,j,k)+w(i,j,k-1))*0.5d0
!--能量的粘性通量-------------------------------
   E1=u1*s11+v1*s12+w1*s13+k0*Tx
   E2=u1*s12+v1*s22+w1*s23+k0*Ty
   E3=u1*s13+v1*s23+w1*s33+k0*Tz
!-------通量=无粘通量+粘性通量---------------------------------------
! 粘性通量=Fv.n=Fv1*s1x+Fv2*s1y+Fv3*s1z---------------
    Flux(2,i,j,k)=(s11*s1x+s12*s1y+s13*s1z)* B%Sk(i,j,k)  
    Flux(3,i,j,k)=(s12*s1x+s22*s1y+s23*s1z)* B%Sk(i,j,k)
    Flux(4,i,j,k)=(s13*s1x+s23*s1y+s33*s1z)* B%Sk(i,j,k)
    Flux(5,i,j,k)=(E1*s1x+E2*s1y+E3*s1z)* B%Sk(i,j,k)
!--------------------------------------------------------------------  
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP  DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
     do m=2,5
      B%Res(m,i,j,k)=B%Res(m,i,j,k)+Flux(m,i,j,k+1)-Flux(m,i,j,k)           
     enddo
	 enddo
     enddo
     enddo
!$OMP END DO

!---记录边界处的粘性力-------------------------------------
!$OMP  DO
   do j=1,ny-1 
   do i=1,nx-1
     B%Surf3(i,j,1)=Flux(2,i,j,1)
     B%Surf3(i,j,2)=Flux(3,i,j,1)
     B%Surf3(i,j,3)=Flux(4,i,j,1)
     B%Surf6(i,j,1)=Flux(2,i,j,nz)
     B%Surf6(i,j,2)=Flux(3,i,j,nz)
     B%Surf6(i,j,3)=Flux(4,i,j,nz)
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


 end subroutine flux_viscous_k








!--------------------------------------------------------------------------
! 原始变量重构  
   subroutine Reconstuction_original(U0,UL,UR,gamma,Iflag_Scheme)
     use  const_var
     use mod_struct_scheme, only: scheme_fP, scheme_fm
  	 implicit none
     real(PRE_EC):: U0(1-LAP:LAP,5) ,UL(5),UR(5),gamma
     integer:: Iflag_Scheme,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1    ;   m=1,5 for d,u,v,w,p
     do m=1,5
       call scheme_fP(UL(m),U0(:,m),Iflag_Scheme) 
       call scheme_fm(UR(m),U0(:,m),Iflag_Scheme) 
     enddo
   end subroutine Reconstuction_original

! 守恒变量重构
   subroutine Reconstuction_conservative(U0,UL,UR,gamma,Iflag_Scheme)
     use  const_var
     use mod_struct_scheme, only: scheme_fP, scheme_fm
   	 implicit none
     real(PRE_EC):: U0(1-LAP:LAP,5),UL(5),UR(5),QL(5),QR(5),gamma
     integer:: Iflag_Scheme,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1   ; m for the conservative variables U0(1,m)=d, U0(2,m)=d*u, ....
     do m=1,5
        call scheme_fP(QL(m),U0(:,m),Iflag_Scheme) 
        call scheme_fm(QR(m),U0(:,m),Iflag_Scheme) 
      enddo
!          find a bug  UL(4)=(QL(4)-(UL(2)*QL(2)+ .... 
       UL(1)=QL(1); UL(2)=QL(2)/UL(1); UL(3)=QL(3)/UL(1); UL(4)=QL(4)/UL(1)    ! density and velocities
       UL(5)=(QL(5)-(UL(2)*QL(2)+UL(3)*QL(3)+UL(4)*QL(4))*0.5d0)*(gamma-1.d0)  ! pressure
       UR(1)=QR(1); UR(2)=QR(2)/UR(1); UR(3)=QR(3)/UR(1); UR(4)=QR(4)/UR(1) 
       UR(5)=(QR(5)-(UR(2)*QR(2)+UR(3)*QR(3)+UR(4)*QR(4))*0.5d0)*(gamma-1.d0)  ! pressure
   end subroutine Reconstuction_conservative


!------------------------------ 特征变量重构 --------------------------------------------------
! 代码由冷岩开发
   subroutine Reconstuction_Characteristic(U0,UL,UR,gamma,Iflag_Scheme)
   use  const_var
   use mod_struct_scheme, only: scheme_fP, scheme_fm
   implicit none
   real(PRE_EC):: U0(1-LAP:LAP,5),V0(1-LAP:LAP,5),UL(5),UR(5),gamma
   real(PRE_EC):: Uh(5),S(5,5),S1(5,5),VL(5),VR(5),QL(5),QR(5)
   real(PRE_EC):: v2,d1,u1,v1,p1,c1,w1,tmp0,tmp1,tmp3,tmp5
   integer:: Iflag_Scheme,i,j,k,m
! U0(k,m) : k=1-LAP,LAP    ; m for the conservative variables U0(1,m)=d, U0(2,m)=d*u, ....
   Uh(:)=0.5d0*(U0(0,:)+U0(1,:))           ! conservative variables in the point I-1/2  (or i)
   d1=Uh(1); u1=Uh(2)/d1; v1=Uh(3)/d1; w1=Uh(4)/d1; p1=(Uh(5)-(Uh(2)*u1+Uh(3)*v1+Uh(4)*w1)*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
   c1=sqrt(gamma*p1/d1)
   v2=(u1*u1+v1*v1+w1*w1)*0.5d0
   tmp1=(gamma-1.d0)/c1
   tmp3=(gamma-1.d0)/(c1*c1)
   tmp5=1.d0/(2.d0*c1)
   tmp0=1.d0/tmp3

! A=S(-1)*LAMDA*S    see 《计算空气动力学》 158-159页   (with alfa1=1, alfa2=0, alfa3=0)
   S(1,1)=V2-tmp0;       S(1,2)=-u1 ;   S(1,3)=-v1 ;        S(1,4)=-w1;                  S(1,5)=1.d0
   S(2,1)=-v1 ;          S(2,2)=0.d0 ;  S(2,3)=1.d0 ;       S(2,4)=0.d0;                 S(2,5)=0.d0 
   S(3,1)=-w1 ;          S(3,2)=0.d0 ;  S(3,3)=0.d0 ;       S(3,4)=1.d0 ;                S(3,5)=0.d0 
   S(4,1)=-u1-V2*tmp1;   S(4,2)=1.d0+tmp1*u1;     S(4,3)=tmp1*v1;     S(4,4)=tmp1*w1;    S(4,5)=-tmp1
   S(5,1)=-u1+V2*tmp1;   S(5,2)=1.d0-tmp1*u1;     S(5,3)=-tmp1*v1;    S(5,4)=-tmp1*w1;   S(5,5)=tmp1 
   
   S1(1,1)=-tmp3;    S1(1,2)=0.d0;   S1(1,3)=0.d0;    S1(1,4)=-tmp5 ;                   S1(1,5)=tmp5
   S1(2,1)=-tmp3*u1; S1(2,2)=0.d0;   S1(2,3)=0.d0;    S1(2,4)=0.5d0-u1*tmp5 ;           S1(2,5)=0.5d0+u1*tmp5
   S1(3,1)=-tmp3*v1; S1(3,2)=1.d0;   S1(3,3)=0.d0;    S1(3,4)=-v1*tmp5;                 S1(3,5)=v1*tmp5
   S1(4,1)=-tmp3*w1; S1(4,2)=0.d0;   S1(4,3)=1.d0;    S1(4,4)=-w1*tmp5;                 S1(4,5)=w1*tmp5
   S1(5,1)=-tmp3*V2; S1(5,2)=v1;     S1(5,3)=w1;      S1(5,4)=(c1*u1-V2-tmp0)*tmp5;     S1(5,5)=(c1*u1+V2+tmp0)*tmp5
! V=SU      V(k)=S*U(k)
   do k=1-LAP,LAP     
     do m=1,5
       V0(k,m)=0.d0
       do j=1,5
         V0(k,m)=V0(k,m)+S(m,j)*U0(k,j)
       enddo
     enddo
   enddo  
   do m=1,5
     call scheme_fP(VL(m),V0(:,m),Iflag_Scheme) 
     call scheme_fm(VR(m),V0(:,m),Iflag_Scheme) 
   enddo
   do m=1,5
     QL(m)=0.d0; QR(m)=0.d0
     do j=1,5
       QL(m)=QL(m)+S1(m,j)*VL(j)
       QR(m)=QR(m)+S1(m,j)*VR(j)
     enddo
   enddo
   UL(1)=QL(1); UL(2)=QL(2)/UL(1)
   UL(3)=QL(3)/UL(1); UL(4)=QL(4)/UL(1)
   UL(5)=(QL(5)-(UL(2)*QL(2)+UL(3)*QL(3)+UL(4)*QL(4))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
  
   UR(1)=QR(1); UR(2)=QR(2)/UR(1)
   UR(3)=QR(3)/UR(1); UR(4)=QR(4)/UR(1)
   UR(5)=(QR(5)-(UR(2)*QR(2)+UR(3)*QR(3)+UR(4)*QR(4))*0.5d0)*(gamma-1.d0)  
   end subroutine Reconstuction_Characteristic
!-----------------------------------------------------------------------------------------------------------





!----------计算粘性系数 (Surthland公式)------------------------
   subroutine get_viscous(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   real(PRE_EC):: Tsb
   integer:: mBlock,i,j,k,nMesh,nx,ny,nz
  Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)
   Tsb=110.4d0/T_inf
    nx=B%nx ; ny=B%ny ; nz=B%nz

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,nx,ny,nz,Tsb,T,Re)
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       B%mu(i,j,k)=1.d0/Re*(1.d0+Tsb)*sqrt(T(i,j,k)**3)/(Tsb+T(i,j,k))
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO

!  粘性系数虚网格上的值 （复制临近面的值）
    B%mu(0,1:ny-1,1:nz-1)=B%mu(1,1:ny-1,1:nz-1)
	B%mu(nx,1:ny-1,1:nz-1)=B%mu(nx-1,1:ny-1,1:nz-1)
	B%mu(1:nx-1,0,1:nz-1)=B%mu(1:nx-1,1,1:nz-1)
	B%mu(1:nx-1,ny,1:nz-1)=B%mu(1:nx-1,ny-1,1:nz-1)
 	B%mu(1:nx-1,1:ny-1,0)=B%mu(1:nx-1,1:ny-1,1)
  	B%mu(1:nx-1,1:ny-1,nz)=B%mu(1:nx-1,1:ny-1,nz-1)
  end subroutine get_viscous

!----------------------------------------------------
!  对湍流粘性系数进行限制 (不能出现负值，不能超过层流粘性系数的MUT_MAX倍）
  subroutine limit_mut(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
!   real(PRE_EC),parameter:: MUT_MAX=200.d0
   integer:: mBlock,i,j,k,nMesh,nx,ny,nz
  Type (Block_TYPE),pointer:: B
  B => Mesh(nMesh)%Block(mBlock)
  nx=B%nx ; ny=B%ny ; nz=B%nz

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(B,nx,ny,nz,MUT_MAX)
     do k=0,nz
     do j=0,ny
     do i=0,nx
       if(B%mu_t(i,j,k) .lt. 0.)  B%mu_t(i,j,k)=0.
       if( MUT_MAX >=0.d0 .and. B%mu_t(i,j,k) > MUT_MAX*B%mu(i,j,k)) B%mu_t(i,j,k)=MUT_MAX*B%mu(i,j,k)
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO
   call Amut_boundary(nMesh,mBlock)
   end subroutine limit_mut




!---------------------------------------------------------
!  利用守恒变量，计算基本量 (d,u,v,T,p,c) 
!  处理计算异常 (例如 温度为负值)
!----------------------------------------------------------
  subroutine comput_duvtpc(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   Type (Block_TYPE),pointer:: B
   integer nMesh,mBlock,nx,ny,nz,i,j,k
   real(PRE_EC) p00  
   
   p00=1.d0/(gamma*Ma*Ma)
   B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
   nx=B%nx; ny=B%ny; nz=B%nz


!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(p00,nx,ny,nz,B,d,uu,v,w,T,p,cc,Ma,Cv)
   do k=1-LAP,nz+LAP-1
     do j=1-LAP,ny+LAP-1
       do i=1-LAP,nx+LAP-1
         d(i,j,k)= B%U(1,i,j,k)
         uu(i,j,k)=B%U(2,i,j,k)/d(i,j,k)
         v(i,j,k)= B%U(3,i,j,k)/d(i,j,k)
         w(i,j,k)= B%U(4,i,j,k)/d(i,j,k)
         T(i,j,k)=(B%U(5,i,j,k)-0.5d0*d(i,j,k)*(uu(i,j,k)*uu(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)))/(Cv*d(i,j,k))
         p(i,j,k)= p00*d(i,j,k)*T(i,j,k)
		 cc(i,j,k)=sqrt(T(i,j,k))/Ma                   ! 声速
       enddo
     enddo
   enddo
!$OMP END PARALLEL  DO 

! Debug Debug message
   if(IF_Debug == 2) then   ! show values at one point 
     if(B%Block_no == Pdebug(1)) then
	  i=Pdebug(2); j=Pdebug(3); k=Pdebug(4)
	  print*, "-----Debug , d,u,v,w,T,p= ----"
	  print*, d(i,j,k),uu(i,j,k),v(i,j,k),w(i,j,k),T(i,j,k),p(i,j,k)
	endif
  endif


  end subroutine comput_duvtpc

!----------------------------------------------------------------------
! 在给定的网格上求解N-S方程 （推进1个时间步）
! 对于单重网格，nMesh=1;  对于多重网格，nMesh=1,2,3, ... 分别对应用细网格、粗网格、更粗网格 ...
! 2015-11-26: A bug in Line 299 is removed   (KRK should be a shared data)
 
  subroutine NS_Time_advance(nMesh)
   use Global_var
   implicit none
   integer:: nMesh
   if(Time_Method .eq. Time_Euler1) then
     call NS_Time_advance_1Euler(nMesh)                    ! 1阶Euler
   else if (Time_Method .eq. Time_LU_SGS ) then            !  LU_SGS
     call NS_Time_advance_LU_SGS(nMesh)
   else if (Time_Method .eq. Time_Dual_LU_SGS) then        ! Dual_LU_SGS
     call  NS_Time_Dual_LU_SGS(nMesh)
   else if (Time_Method .eq. Time_RK3) then
     call NS_Time_advance_RK3(nMesh)                       ! 3阶RK
   else
     print*, "This time advance method is not supported!!!"
   endif
    call force_vt_kw(nMesh)     ! 强制 vt, k,w 非负

  end subroutine NS_Time_advance

!---------------------------------------------------------------------------------------------
! 强制vt, k,w非负   
   subroutine force_vt_kw(nMesh)
    use Global_var
    implicit none
    integer:: nMesh,mBlock,nx,ny,nz,i,j,k
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
  
    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny ; nz=B%nz
    if(MP%NVAR == 6) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
      do k=1,nz-1
 	  do j=1,ny-1
      do i=1,nx-1
       if(B%U(6,i,j,k) < 0)  B%U(6,i,j,k)=1.d-10
      enddo
      enddo
	  enddo
!$OMP END PARALLEL DO
    else if (MP%NVAR == 7) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
 	  do k=1,nz-1
	  do j=1,ny-1
      do i=1,nx-1
       if(B%U(6,i,j,k) < 0)  B%U(6,i,j,k)=1.d-10
       if(B%U(7,i,j,k) < 0)  B%U(7,i,j,k)=1.d-10
	  enddo
      enddo
	  enddo
!$OMP END PARALLEL DO
    endif
   enddo    
  end
!--------------------------------------------------------------------------------------











! 采用 LU_SGS方法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
  subroutine NS_Time_advance_LU_SGS(nMesh)
   use Global_var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: du
   call Set_Un(nMesh)
   call Comput_Residual_one_mesh(nMesh)              ! 单重网格上计算残差 (以及Du)
   if(nMesh .ne. 1) call Add_force_function(nMesh)   !  添加强迫函数（多重网格的粗网格使用）
  
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz
!--------------------------------------------------------------------------------------
!   时间推进 

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
             B%U(m,i,j,k)=B%Un(m,i,j,k)+B%dU(m,i,j,k)           ! LU_SGS方法
            enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO       
  
  enddo

!----------------------------------------------------------------   
    if( IFLAG_LIMIT_Flow == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 

!---------------------------------------------------------------------------------------  
   call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
   
   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global            ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1              ! 计算步数

  end subroutine NS_Time_advance_LU_SGS
!--------------------------------------------------------------------------------------



!  采用双时间步长法 LU_SGS方法进行时间推进一个时间步 
!  目前Dual LU_SGS 方法尚不支持多重网格,因而nMesh只能为1

  subroutine NS_Time_Dual_LU_SGS(nMesh)
    use Global_var
    use mod_struct_bc, only: Boundary_condition_onemesh
    use mod_struct_mpi, only: update_buffer_onemesh
    implicit none
    integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz,Kt_in
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
    real(PRE_EC):: max_res
    
	MP=>Mesh(nMesh)
    NVAR1=MP%NVAR
 do kt_in=1, step_inner_Limit                      ! 内循环迭代

   call Comput_Residual_one_mesh(nMesh)              ! 单重网格上计算残差及Du
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
             B%U(m,i,j,k)=B%U(m,i,j,k)+B%dU(m,i,j,k)           ! LU_SGS方法
            enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO       
   enddo

!----------------------------------------------------------------   
    if( IFLAG_LIMIT_FLOW == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 


  call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
  call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
  call comput_max_Res_onemesh(nMesh)                 ! 计算最大残差及均方根残差

     max_res=MP%Res_rms(1)       ! 最大均方根残差 (作为内迭代标准)
     do m=1,NVAR1
 	  max_res=max(max_res,MP%Res_rms(m))
	 enddo
     if( max_res .le. Res_Inner_Limit) exit   ! 达到残差标准，跳出内迭代
 enddo
   
   if(my_id .eq. 0) then 	
	 print*, "Inner step ... ", kt_in
	 print*, "rms residual eq =", MP%Res_rms(1:NVAR1)
   endif

    
 
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
			 B%Un1(m,i,j,k)=B%Un(m,i,j,k)        
             B%Un(m,i,j,k)=B%U(m,i,j,k)  
            enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO       
   
  enddo


   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global            ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1              ! 计算步数

  end subroutine NS_Time_Dual_LU_SGS
!--------------------------------------------------------------------------------------












! 采用1阶Euler法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
  subroutine NS_Time_advance_1Euler(nMesh)
   use Global_var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: du
   call Set_Un(nMesh)
   call Comput_Residual_one_mesh(nMesh)              ! 单重网格上计算残差
   if(nMesh .ne. 1) call Add_force_function(nMesh)   !  添加强迫函数（多重网格的粗网格使用）

    NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz
!--------------------------------------------------------------------------------------
!   时间推进 
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)
     do k=1,nz-1
       do j=1,ny-1
         do i=1,nx-1
           do m=1,NVAR1
             du=B%Res(m,i,j,k)/B%vol(i,j,k)   
             B%U(m,i,j,k)=B%Un(m,i,j,k)+B%dt(i,j,k)*du
           enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
  
  enddo    


!----------------------------------------------------------------   
    if( IFLAG_LIMIT_FLOW == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 

!---------------------------------------------------------------------------------------  
   call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global            ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1              ! 计算步数

  end subroutine NS_Time_advance_1Euler
!----------------------------------------------------------------------------------------


! 采用3阶RK方法推进1个时间步 （第nMesh重网格 的单重网格）
  subroutine NS_Time_advance_RK3(nMesh)
   use Global_var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m,nx,ny,nz
   Type (Block_TYPE),pointer:: B
   real(PRE_EC):: du
 
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)

!$OMP PARALLEL DO PRIVATE(i,j,k,m) SHARED(NVAR1,B)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%Un(m,i,j,k)=B%U(m,i,j,k)
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
 
   enddo

   do KRK=1,3                                          ! 3-step Runge-Kutta Method
	 call Comput_Residual_one_mesh(nMesh)              ! 计算残差
     if(nMesh .ne. 1) call Add_force_function(nMesh)   ! 添加强迫函数（多重网格的粗网格使用）
	 do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)                  ! 第nMesh 重网格的第mBlock块
       nx=B%nx; ny=B%ny; nz=B%nz
!--------------------------------------------------------------------------------------
!    时间推进

!$OMP PARALLEL DO PRIVATE(i,j,k,m,du) SHARED(NVAR1,nx,ny,nz,Ralfa,Rbeta,Rgamma,B,KRK)
       do k=1,nz-1 
         do j=1,ny-1
           do i=1,nx-1
             do m=1,NVAR1
		       du=B%Res(m,i,j,k)/B%Vol(i,j,k)  
               B%U(m,i,j,k)=Ralfa(KRK)*B%Un(m,i,j,k)+Rgamma(KRK)*B%U(m,i,j,k)+B%dt(i,j,k)*Rbeta(KRK)*du        ! 3阶RK
             enddo
           enddo
         enddo
	   enddo
 !$OMP END PARALLEL DO 
   enddo    

!---------------------------------------------------------------------------------------

    if( IFLAG_LIMIT_FLOW == 1) then                      ! 对压力、密度进行限制
	  call limit_flow(nMesh)
	endif 

 
     call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
     call update_buffer_onemesh(nMesh)              ! 同步各块的交界区
   enddo   
   Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global          ! 时间 （使用全局时间步长法时有意义）
   Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1            ! 计算步数

  end subroutine NS_Time_advance_RK3


! 计算最大残差和均方根残差（整个网格）
  subroutine comput_max_Res_onemesh(nMesh)
   use Global_var
   use, intrinsic :: ieee_arithmetic
   implicit none
	integer:: nMesh,mBlock,i,j,k,m,ierr
	logical F_NaN
 	real(PRE_EC):: Res,Res_max(7),Res_rms(7)
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
 
     MP=> Mesh(nMesh)
   
      Res_max(:)=0.d0
	  Res_rms(:)=0.d0
 
   do mBlock=1,MP%NUM_BLOCK 
!     call comput_max_Res_oneblock(nMesh,mBlock)
   	  B => MP%Block(mBlock)                 !第nMesh 重网格的第mBlock块

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(MP,B) REDUCTION(MAX: Res_max) REDUCTION(+: Res_rms)   
	 do k=1,B%nz-1
       do j=1,B%ny-1
         do i=1,B%nx-1
! -------------------------------------------------------------------------------------------
!    时间推进
           do m=1,MP%NVAR
             Res=B%Res(m,i,j,k)
!--------------------------------------------------------------------------------------------------
! detech "NaN", Since Ver 0.72, which is useful for debug
             F_NaN=ieee_is_nan(Res)
             if(F_NaN) then
 		       print*, "NaN in Residual is found !, In block",B%block_no
			   print*, "location i,j,k,m=",i,j,k,m
!               B%Res(m,i,j,k)=0.d0    ! 强制为0
			    print*, "Stop"
			   stop
		     endif  
              Res_max(m)=max(Res_max(m),abs(Res))    ! 最大残差
			  Res_rms(m)=Res_rms(m)+Res*Res           ! 均方根残差
!--------------------------------------------------------------------------------------------------       
	       enddo
         enddo
       enddo
	 enddo
!$OMP END PARALLEL DO
  enddo
 
    call MPI_ALLREDUCE(Res_max(1),MP%Res_max(1),MP%NVAR,OCFD_DATA_TYPE,MPI_MAX,Struct_Comm,ierr)
    call MPI_ALLREDUCE(Res_rms(1),MP%Res_rms(1),MP%NVAR,OCFD_DATA_TYPE,MPI_SUM,Struct_Comm,ierr)
    MP%Res_rms(:)=sqrt(MP%Res_rms(:)/(MP%Num_Cell))   !均方根残差

  end  subroutine comput_max_Res_onemesh

!-------------------------------------------------------------

!--------------------------------------------------------------
! 打印残差（最大残差和均方根残差）
  subroutine output_Res(nMesh)
   use Global_var
   implicit none
	integer:: nMesh
    call   comput_max_Res_onemesh(nMesh)
!-----------------------------------
   if(my_id .eq. 0) then
    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
    print*, "----------The Max Residuals are-------- ", " ---Mesh---",nMesh
    write(*, "(7E20.10)") Mesh(nMesh)%Res_max(:)
    print*, "  The R.M.S Residuals are "
    write(*, "(7E20.10)") Mesh(nMesh)%Res_rms(:)
    open(99,file="Residual.dat",position="append")
    write(99,"(I8,15E20.10)") Mesh(nMesh)%Kstep, Mesh(nMesh)%Res_max(:),Mesh(nMesh)%Res_rms(:)
    close(99) 
   endif

  end  subroutine output_Res

!-------------------------------------------------------------


!----------------------------------------------------------
! 对SA,SST方程的物理量进行限制
  subroutine limit_vt(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   Type (Block_TYPE),pointer:: B
   integer nMesh,mBlock,NVAR1,nx,ny,nz,i,j,k
   
   B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
   nx=B%nx; ny=B%ny; nz=B%nz
   NVAR1=Mesh(nMesh)%NVAR
   if(NVAR1 .eq. 6) then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B)
    do k=0,nz
    do j=0,ny
	do i=0,nx
	 if(B%U(6,i,j,k) .lt. 0.d0) B%U(6,i,j,k)=0.d0
	enddo
	enddo
	enddo
!$OMP END PARALLEL DO
   else if (NVAR1 .eq. 7) then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B)
    do k=0,nz
	do j=0,ny
	do i=0,nx
	 if(B%U(6,i,j,k) .lt. 0.d0) B%U(6,i,j,k)=0.d0
	 if(B%U(7,i,j,k) .lt. 0.d0) B%U(7,i,j,k)=0.d0
    enddo
	enddo
	enddo
!$OMP END PARALLEL DO

  endif
  end	 







!----------------------------------------------------------------------
! 多重网格求解N-S方程 （推进1个时间步）
! nMesh=1,2,3 分别对应用细网格、粗网格、更粗网格
! 包括2重网格和3重网格两个子程序；
! Code by Li Xinliang & Leng Yan
!---------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! 两重网格上推进1个时间步 (3阶RK or 1th Euler)
  subroutine NS_2stge_multigrid
   use Global_var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   implicit none
   integer::nMesh,m
   Type (Block_TYPE),pointer:: B
   integer,parameter:: Time_step_coarse_mesh=3       ! 粗网格迭代步数
!---------------------------------------------------
! -------------------------  网格1 -----------------
   if(Time_Method .eq. Time_Euler1) then
	 call  NS_Time_advance_1Euler(1)                 ! 细网格，1阶Euler方法推进1步 -> U(n+1)
   else
	 call  NS_Time_advance_RK3(1)                    ! 细网格，RK方法推进1步 -> U(n+1)
   endif 
   call  Comput_Residual_one_mesh(1)                 ! 计算网格1的残差 R(n+1)  
   call  interpolation2h(1,2,2)                      ! 把残差插值到网格2 (储存在QF里面)
   call  interpolation2h(1,2,1)                      ! 把守恒变量从网格1插值到网格2   （flag=1 插值守恒变量，=2 插值残差）
!------------------------------
   call  Boundary_condition_onemesh(2)               ! 物理边界条件
   call  update_buffer_onemesh(2)                    ! 内边界条件
   call  Comput_Residual_one_mesh(2)                 ! 计算网格2的残差
   call  comput_force_function(2)                    ! 计算强迫函数QF
   if(Time_Method .eq. Time_Euler1) then
	 call Set_Un(2)                                  ! 记录初始值  （RK方法中已经包含了该步） 
     do m=1, Time_step_coarse_mesh
	   call  NS_Time_advance_1Euler(2)               ! 1阶Euler迭代若干步
     enddo
   else 
	 call  NS_Time_advance_RK3(2)                    ! RK方法推进1步 （网格2）
   endif
   call  comput_delt_U(2)                            ! 计算修正量deltU （储存在Un里面）
   call  prolong_U(2,1,2)                            ! 把修正量插值到细网格 (储存在Un里面); flag=2 插值deltU (储存在Un里)
!------------------------------------	 
   call  comput_new_U(1)                             ! 计算新的U  (U=U+deltU)
   call  Boundary_condition_onemesh(1)               ! 物理边界条件
   call  update_buffer_onemesh(1)                    ! 内边界条件

  end subroutine NS_2stge_multigrid
!------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! 三重网格上迭代1个时间步 （V-型迭代） 3阶RK or 1阶Euler
  subroutine NS_3stge_multigrid
   use Global_var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   implicit none
   integer::nMesh,m
   integer,parameter:: Time_step_coarse_mesh=3       ! 粗网格迭代步数 (对1阶Euler有效)
!---------------------------------------------------
! ---- ---------------------------------- 网格1 -----------------
   if(Time_Method .eq. Time_Euler1) then
	 call  NS_Time_advance_1Euler(1)                 ! 细网格，1阶Euler方法推进1步 -> U(n+1)
   else
	 call  NS_Time_advance_RK3(1)                    ! 细网格，RK方法推进1步 -> U(n+1)
   endif
   call  Comput_Residual_one_mesh(1)                 ! 计算网格1的残差 R(n+1)   ! ????? 该步似乎可以省略 ?????  
! -----------------------------  
   call  interpolation2h(1,2,2)                      ! 把残差插值到网格2 (储存在网格2的QF里面)
   call  interpolation2h(1,2,1)                      ! 把守恒变量从网格1插值到网格2 （储存到U里面）  （flag=1 插值守恒变量，=2 插值残差）
!-------网格2 --------------------
   call  Boundary_condition_onemesh(2)               ! 物理边界条件
   call  update_buffer_onemesh(2)                    ! 内边界条件
   call  Comput_Residual_one_mesh(2)                 ! 计算网格2的残差         Res_2h(0)
   call  comput_force_function(2)                    ! 计算强迫函数QF （网格2）QF_2h=QF_2h-Res_2h(0) 
   if(Time_Method .eq. Time_Euler1) then
	 call Set_Un(2)                                  ! 记录初始值  （RK方法中已经包含了该步） 
     do m=1, Time_step_coarse_mesh
	   call  NS_Time_advance_1Euler(2)               ! 1阶Euler迭代若干步
     enddo
   else 
	 call  NS_Time_advance_RK3(2)                    ! RK方法推进1步 （网格2）
   endif
   call  Comput_Residual_one_mesh(2)                 ! 计算网格2的残差 R_2h(n+1) 
   call  Add_force_function(2)                       ! 添加上强迫残差(储存在Res里面)  RF_2h(n+1)=R_2h(n+1)+QF_2h  ;  目的：插值到网格3上
   call  interpolation2h(2,3,2)                      ! 把残差插值到网格3 (储存在网格3的QF里面)
   call  interpolation2h(2,3,1)                      ! 把守恒变量从网格2插值到网格3 （储存到U里面）  （flag=1 插值守恒变量，=2 插值残差）
!------网格3----------------------	  
   call  Boundary_condition_onemesh(3)               ! 边界条件: 物理边界 
   call  update_buffer_onemesh(3)                    ! 内边界
   call  Comput_Residual_one_mesh(3)                 ! 计算网格3的残差
   call  comput_force_function(3)                    ! 计算强迫函数QF （网格3）
   if(Time_Method .eq. Time_Euler1) then
	 call Set_Un(3)
     do m=1, Time_step_coarse_mesh
	   call  NS_Time_advance_1Euler(3)               ! 1阶Euler迭代若干步
     enddo
   else 
	 call  NS_Time_advance_RK3(3)                    ! RK方法推进1步 （网格3）
   endif
   call  comput_delt_U(3)                            ! 计算修正量deltU (=U-Un)
   call  prolong_U(3,2,2)                            ! 把修正量插值到网格2 (储存在deltU里面); flag=2 插值deltU 
!------网格2------------------------      
   call  comput_new_U(2)                             ! 网格2计算新的U  (U=U+deltU)
   call  comput_delt_U(2)                            ! 计算修正量deltU =U-Un
   call  prolong_U(2,1,2)                            ! 把修正量插值到细网格 (储存在deltU里面); flag=2 插值deltU 
!------网格1------------------------------
   call  comput_new_U(1)                             ! 计算新的U  (U=U+deltU)
   call Boundary_condition_onemesh(1)                ! 物理边界条件
   call update_buffer_onemesh(1)                     ! 内边界条件

  end subroutine NS_3stge_multigrid

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
  
!  计算强迫函数 QF=Ih_to_2h Res(n-1) - Res(n)        ! QF中储存着细网格插值过来的残差
  subroutine comput_force_function(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m
   Type (Block_TYPE),pointer:: B
    NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,NVAR1)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%QF(m,i,j,k)=B%QF(m,i,j,k)-B%Res(m,i,j,k)            ! QF原先储存着从细网格插值过来的残差
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo

  end  subroutine comput_force_function

!------------------------------------------------------------
!  把强迫函数添加到残差中 RF=R+QF        
  subroutine Add_force_function(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m
   Type (Block_TYPE),pointer:: B
    NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,NVAR1)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%Res(m,i,j,k)=B%Res(m,i,j,k)+B%QF(m,i,j,k)            ! 添加强迫函数后的残差仍储存在B%Res里面 （节省内存）
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo
 
  end  subroutine Add_force_function

!----------------------------------------------------------------------  
!  计算修正量 deltU=U-Un
  subroutine comput_delt_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,5
	         B%deltU(m,i,j,k)=B%U(m,i,j,k)-B%Un(m,i,j,k)
           enddo
	     enddo
       enddo
	 enddo
 !$OMP END PARALLEL DO 
  enddo
  end  subroutine comput_delt_U
!-----------------------------------------------------------------------
! 设定Un=U
  subroutine Set_Un(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,NVAR1,i,j,k,m
   Type (Block_TYPE),pointer:: B
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B,NVAR1)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,NVAR1
	         B%Un(m,i,j,k)=B%U(m,i,j,k)
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo
  
  end  subroutine Set_Un
!-------------------------------修正U --------------------------------------
  subroutine comput_new_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B)
	 do k=-1,B%nz+1
       do j=-1,B%ny+1
	     do i=-1,B%nx+1
	       do m=1,5                                             ! 第6个量是湍流粘性系数vt (SA模型使用),不需要修正
	         B%U(m,i,j,k)=B%U(m,i,j,k)+B%deltU(m,i,j,k)         ! Un里面储存的是U的修正量 （从粗网格插值而来）
           enddo
	     enddo
       enddo
	 enddo
!$OMP END PARALLEL DO 
   enddo

  end  subroutine comput_new_U
!---------------------------------------------------------------------------
! 粗网格向细网格的插值(Prolong) 及 细网格向粗网格上插值 (interpolation)
!----------------------------------------------------------------------
! 将网格m1的守恒变量(U) 或U的差插值到网格m2 (上一级细网格)
! flag=1时，将U插值到上一级网格；  (准备初值时使用)
! flag=2时，将deltU插值到上一级网格 （deltU储存着本时间步与上个时间步U的差） 

  Subroutine prolong_U(m1,m2,flag)
   use Global_Var
   implicit none
   integer:: m1,m2,mb,flag
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   if(m1 .le. 1 .or. m1-m2 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2)
     do mb=1,MP1%Num_Block
       B1=>MP1%Block(mb)
	   B2=>Mp2%Block(mb)
	   if(flag .eq. 1) then
!	   call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%U(1,-1,-1,-1),B2%U(1,-1,-1,-1))   ! 旧的程序接口，与新版Fortran不兼容
	    call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%U,B2%U)
  	   else
!	    call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%deltU(1,-1,-1,-1),B2%deltU(1,-1,-1,-1))
	    call prolongation(B1%nx,B1%ny,B1%nz,B2%nx,B2%ny,B2%nz,B1%deltU,B2%deltU)
     endif
   enddo

  end Subroutine prolong_U
!-------------------------------------------------------------
!   粗网格向细网格上的插值     
!   U1是粗网格上的值； U2是细网格上的值
  subroutine prolongation(nx1,ny1,nz1,nx2,ny2,nz2,U1,U2)
   use precision_EC
   implicit none
    integer:: i,j,k,m,nx1,ny1,nz1,nx2,ny2,nz2,NV
    real(PRE_EC),dimension(:,:,:,:),pointer:: U1,U2
 !   real(PRE_EC):: U1(NVAR,-1:nx1+1,-1:ny1+1,-1:nz1+1),U2(NVAR,-1:nx2+1,-1:ny2+1,-1:nz2+1)
    integer:: ia(2,0:nx2),ja(2,0:ny2),ka(2,0:nz2),U_bound(4)
 !   integer,parameter::NVAR=6
    real(PRE_EC),parameter:: a1=27.d0/64.d0,a2=9.d0/64.d0,a3=3.d0/64.d0,a4=1.d0/64.d0   ! 插值系数
    
!  寻找插值基架点的下标 
!  ia(1,i) 是距离i点最近的粗网格点的下标；ia(2,i)是次近点的下标	 
    U_bound=UBOUND(U1)   ! 第1维的上界 （NVAR)
    NV=U_bound(1)   ! NVAR= 5, 6 or 7 
    
   do i=0,nx2
	 if(mod(i,2).eq.0) then
	   ia(1,i)=i/2                    !最近点
	   ia(2,i)=i/2+1                  !次近点
	 else  
	   ia(1,i)=i/2+1                  !最近点
	   ia(2,i)=i/2                    !次近点
	 endif
   enddo
   do j=0,ny2
	 if( mod(j,2).eq. 0) then
	   ja(1,j)=j/2
	   ja(2,j)=j/2+1
	 else
	   ja(1,j)=j/2+1
	   ja(2,j)=j/2
	 endif
   enddo
   do k=0,nz2
	 if(mod(k,2).eq.0) then
	   ka(1,k)=k/2                    !最近点
	   ka(2,k)=k/2+1                  !次近点
	 else  
	   ka(1,k)=k/2+1                  !最近点
	   ka(2,k)=k/2                    !次近点
	 endif
   enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,m)

   do k=0,nz2
	 do j=0,ny2
	   do i=0,nx2
	     do m=1,NV
!               插值，最近点的权重a1, 次近点的权重a2, 最远点的权重a3	 
	       U2(m,i,j,k)=a1*U1(m,ia(1,i),ja(1,j),ka(1,k))+a2*(U1(m,ia(2,i),ja(1,j),ka(1,k))+U1(m,ia(1,i),ja(2,j),ka(1,k)) &
	                   +U1(m,ia(1,i),ja(1,j),ka(2,k)))+a3*(U1(m,ia(2,i),ja(1,j),ka(2,k))+U1(m,ia(1,i),ja(2,j),ka(2,k)) &
	                   +U1(m,ia(2,i),ja(2,j),ka(1,k)))+a4*U1(m,ia(2,i),ja(2,j),ka(2,k))
         enddo
	   enddo
     enddo
   enddo
!$OMP END PARALLEL DO 
  end subroutine prolongation
!---------------------------------------------------------------------------------
! 将网格m1的守恒变量U插值到网格m2 (细网格->粗网格) 
  Subroutine interpolation2h(m1,m2,flag)
   use Global_Var
   implicit none
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   real(PRE_EC),dimension(:,:,:,:),pointer:: P1,P2
   integer:: NVAR1,flag,m1,m2,mb,i,j,k,m,i1,i2,j1,j2,k1,k2
!  flag==1 插值守恒变量； flag==2 插值残差
   if( m2-m1 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2) 
	 NVAR1=5              ! 只插值5个守恒变量  
     do mb=1,MP1%Num_Block
       B1=>MP1%Block(mb)
  	   B2=>Mp2%Block(mb)
       if(flag .eq. 1) then  ! 插值守恒变量
	     P1=>B1%U
	     P2=>B2%U
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B1,B2,NVAR1,P1,P2)
	     do k=1,B2%nz-1
           do j=1,B2%ny-1
	         do i=1,B2%nx-1
	           i1=2*i-1 ; i2=2*i
	           j1=2*j-1 ; j2=2*j
	           k1=2*k-1 ; k2=2*k
               do m=1,NVAR1
!     以控制体体积为权重的加权平均 	 
	             P2(m,i,j,k)=(P1(m,i1,j1,k1)*B1%Vol(i1,j1,k1)+P1(m,i1,j2,k1)*B1%Vol(i1,j2,k1)   &
		                      +P1(m,i1,j2,k2)*B1%Vol(i1,j2,k2)+P1(m,i2,j1,k1)*B1%Vol(i2,j1,k1)   &
	                          +P1(m,i2,j1,k2)*B1%Vol(i2,j1,k2)+P1(m,i2,j2,k1)*B1%Vol(i2,j2,k1)   &
					          +P1(m,i2,j2,k2)*B1%Vol(i2,j2,k2)+P1(m,i1,j1,k2)*B1%Vol(i1,j1,k2))/B2%Vol(i,j,k)
  	           enddo
	         enddo
	       enddo
	     enddo
!$OMP END PARALLEL DO 

       else     ! 插值残差  （把m1网格上的残差B%Res 插值到m2网格上B%QF (然后减去本m2网格上的残差，形成强迫函数)）
   	     P1=>B1%Res
	     P2=>B2%QF

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(B2,P1,P2,NVAR1)
	     do k=0,B2%nz
           do j=0,B2%ny
	         do i=0,B2%nx
	           i1=2*i-1 ; i2=2*i
	           j1=2*j-1 ; j2=2*j
	           k1=2*k-1 ; k2=2*k
               do m=1,NVAR1
	             P2(m,i,j,k)=P1(m,i1,j1,k1)+P1(m,i1,j2,k1)+P1(m,i1,j2,k2)+P1(m,i2,j1,k1)   &    ! 残差的插值： 简单相加
		                     +P1(m,i2,j1,k2)+P1(m,i2,j2,k1)+P1(m,i2,j2,k2)+P1(m,i1,j1,k2)
  	           enddo
	         enddo
	       enddo
	     enddo
 !$OMP END PARALLEL DO 

 	   endif
     enddo

  end Subroutine interpolation2h

!------------------------------------------------------------------------------
! Source term of SA model  (See: J. Blazek's book, P240-243) 
! Do not consider transition (full turbulence)
! Ver 0.81a, using Eq. (7.4.2) (J. Blazek's book)
! Ver 0.81b, limit for sorce term  ( not less than 0.3*Omega)
! Ver 0.98c, Nondimensional (See CFL3D manual)

  subroutine  Turbulence_model_SA(nMesh,mBlock)
   use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,nMesh
   real(PRE_EC):: ui,vi,wi,uj,vj,wj,uk,vk,wk,ux,vx,wx,uy,vy,wy,uz,vz,wz
   real(PRE_EC):: s1x,s1y,s1z,v0,vti,vtj,vtk,vtx,vty,vtz
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz
   real(PRE_EC):: S, S1,X,fv1,fv2,fv3,ft1,ft2,r,g,fw,Q_SA,vn1,vn2,vfi,Q3
   
   real(PRE_EC),parameter:: SA_sigma=2.d0/3.d0,Cv1=7.1d0,Cv2=5.d0,Cb1=0.1355d0,Cb2=0.622d0,SA_k=0.41d0
   real(PRE_EC),parameter:: Cw1=Cb1/(SA_k*SA_k)+(1.d0+Cb2)/SA_sigma, Cw2=0.3d0, Cw3=2.d0, &
                            Ct1=1.d0, Ct2=2.d0, Ct3=1.2d0, Ct4=0.5d0  !Ct3=1.3d0
   real(PRE_EC),Pointer,dimension(:,:,:):: vt,Fluxv,fluxv2
   Type (Block_TYPE),pointer:: B
 
  
   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz



 ! 计算湍流粘性系数
   allocate(vt(0:nx,0:ny,0:nz),Fluxv(nx,ny,nz),fluxv2(nx,ny,nz))
	
! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量



!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Re)  
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)*Re             ! 量纲转换
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO   



!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,vt,d)
 
    do k=0,nz
    do j=0,ny
    do i=0,nx
     vt(i,j,k)=B%U(6,i,j,k)
	 X=d(i,j,k)*vt(i,j,k)/B%mu(i,j,k)   ! 湍流粘性系数与层流粘性系数之比
     fv1=X**3/(X**3+Cv1**3)
     B%mu_t(i,j,k)=fv1*d(i,j,k)*vt(i,j,k)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
   
   ! 限定湍流粘性系数
   call limit_mut(nMesh,mBlock)

! 设定湍流粘性系数虚网格的值
 

! 计算vt方程的残差 B%Res(6,:,:,:)

!$OMP PARALLEL  DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,Re,vt,d,uu,v,w,Fluxv,fluxv2)
!------i- direcion ---------------------------
!$OMP DO
    do k=1,nz-1 
     do j=1,ny-1
      do i=1,nx
          s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)  ! 归一化的法方向
          vti=vt(i,j,k)-vt(i-1,j,k)                               ! SA模型中的vt
          vtj=0.25d0*(vt(i,j+1,k)-vt(i,j-1,k)+vt(i-1,j+1,k)-vt(i-1,j-1,k))
          vtk=0.25d0*(vt(i,j,k+1)-vt(i,j,k-1)+vt(i-1,j,k+1)-vt(i-1,j,k-1))
          ix=B%ix1(i,j,k); iy=B%iy1(i,j,k); iz=B%iz1(i,j,k)
          jx=B%jx1(i,j,k); jy=B%jy1(i,j,k); jz=B%jz1(i,j,k)
          kx=B%kx1(i,j,k); ky=B%ky1(i,j,k); kz=B%kz1(i,j,k)
          vtx=vti*ix+vtj*jx+vtk*kx
          vty=vti*iy+vtj*jy+vtk*ky
          vtz=vti*iz+vtj*jz+vtk*kz
          v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i-1,j,k)+B%mu(i-1,j,k)/d(i-1,j,k))  ! (I-1/2,J,K) 点的动力学粘性系数
 
          vn1=uu(i-1,j,k)*s1x+v(i-1,j,k)*s1y+w(i-1,j,k)*s1z   ! 法向速度
          vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
          vfi=0.5d0*((vn1+abs(vn1))*vt(i-1,j,k)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
          Fluxv(i,j,k)= (-vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))* B%Si(i,j,k)    !!! Re
          Fluxv2(i,j,k)=(vtx*s1x+vty*s1y+vtz*s1z)*B%Si(i,j,k)
	  enddo
     enddo
    enddo
!$OMP END DO
!$OMP DO
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
	     v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
         B%Res(6,i,j,k)=Fluxv(i+1,j,k)-Fluxv(i,j,k) -v0/Re*(Fluxv2(i+1,j,k)-Fluxv2(i,j,k))    !!! Re     
      enddo
      enddo
      enddo
!$OMP END DO
!-----j- direction ------------------------
!$OMP DO
    do k=1,nz-1 
    do j=1,ny
    do i=1,nx-1
      s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向
      vti=0.25d0*(vt(i+1,j,k)-vt(i-1,j,k)+vt(i+1,j-1,k)-vt(i-1,j-1,k))
      vtj=vt(i,j,k)-vt(i,j-1,k)
      vtk=0.25d0*(vt(i,j,k+1)-vt(i,j,k-1)+vt(i,j-1,k+1)-vt(i,j-1,k-1))
      ix=B%ix2(i,j,k); iy=B%iy2(i,j,k); iz=B%iz2(i,j,k)
      jx=B%jx2(i,j,k); jy=B%jy2(i,j,k); jz=B%jz2(i,j,k)
      kx=B%kx2(i,j,k); ky=B%ky2(i,j,k); kz=B%kz2(i,j,k)
      vtx=vti*ix+vtj*jx+vtk*kx
      vty=vti*iy+vtj*jy+vtk*ky
      vtz=vti*iz+vtj*jz+vtk*kz
      v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i,j-1,k)+B%mu(i,j-1,k)/d(i,j-1,k))  ! (I-1/2,J,K) 点的动力学粘性系数   s11=Amu1*(tmp1*ux-tmp2*(vy+wz))    ! tmp1=4.d0/3.d0; tmp2=2.d0/3.d0
      vn1=uu(i,j-1,k)*s1x+v(i,j-1,k)*s1y+w(i,j-1,k)*s1z   ! 法向速度
      vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
      vfi=0.5d0*((vn1+abs(vn1))*vt(i,j-1,k)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
      Fluxv(i,j,k)= (-vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))* B%Sj(i,j,k)   !!! Re
      Fluxv2(i,j,k)= (vtx*s1x+vty*s1y+vtz*s1z)* B%Sj(i,j,k)
    enddo
    enddo
    enddo
!$OMP END DO
!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
	  v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
      B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxv(i,j+1,k)-Fluxv(i,j,k)- v0/Re*(Fluxv2(i,j+1,k)-Fluxv2(i,j,k))   !!! Re        
     enddo
     enddo
     enddo
!$OMP END DO
!-----k- direction ------------------------
!$OMP DO
   do k=1,nz 
   do j=1,ny-1
   do i=1,nx-1
    s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向
    vti=0.25d0*(vt(i+1,j,k)-vt(i-1,j,k)+vt(i+1,j,k-1)-vt(i-1,j,k-1))
    vtj=0.25d0*(vt(i,j+1,k)-vt(i,j-1,k)+vt(i,j+1,k-1)-vt(i,j-1,k-1))
    vtk=vt(i,j,k)-vt(i,j,k-1)

      ix=B%ix3(i,j,k); iy=B%iy3(i,j,k); iz=B%iz3(i,j,k)
      jx=B%jx3(i,j,k); jy=B%jy3(i,j,k); jz=B%jz3(i,j,k)
      kx=B%kx3(i,j,k); ky=B%ky3(i,j,k); kz=B%kz3(i,j,k)
      vtx=vti*ix+vtj*jx+vtk*kx
      vty=vti*iy+vtj*jy+vtk*ky
      vtz=vti*iz+vtj*jz+vtk*kz
      v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i,j,k-1)+B%mu(i,j,k-1)/d(i,j,k-1))  ! (I,J,K-1/2) 点的动力学粘性系数
     vn1=uu(i,j,k-1)*s1x+v(i,j,k-1)*s1y+w(i,j,k-1)*s1z   ! 法向速度
     vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
     vfi=0.5d0*((vn1+abs(vn1))*vt(i,j,k-1)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
     Fluxv(i,j,k)= ( -vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
     Fluxv2(i,j,k)= (vtx*s1x+vty*s1y+vtz*s1z)*B%Sk(i,j,k)   
    enddo
    enddo
    enddo
!$OMP END DO
!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
 	   v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
       B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxv(i,j,k+1)-Fluxv(i,j,k)-v0/Re*(Fluxv2(i,j,k+1)-Fluxv2(i,j,k))           
	 enddo
     enddo
     enddo
!$OMP END DO
!--------------源项---------------------------------------------------------
!$OMP DO
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
!----- get S (normal of vorticity)  S=sqrt(0.5*Omiga_ij*Omiga_ij) at the cell's center ------
!  计算涡量；计算湍流粘性系数的梯度
      

   ui=uu(i+1,j,k)-uu(i-1,j,k)            
   vi=v(i+1,j,k)-v(i-1,j,k)  
   wi=w(i+1,j,k)-w(i-1,j,k)  

   uj=uu(i,j+1,k)-uu(i,j-1,k)   
   vj=v(i,j+1,k)-v(i,j-1,k)
   wj=w(i,j+1,k)-w(i,j-1,k) 

   uk=uu(i,j,k+1)-uu(i,j,k-1)  
   vk=v(i,j,k+1)-v(i,j,k-1)
   wk=w(i,j,k+1)-w(i,j,k-1)  
   
   ix=B%ix0(i,j,k); iy=B%iy0(i,j,k); iz=B%iz0(i,j,k)
   jx=B%jx0(i,j,k); jy=B%jy0(i,j,k); jz=B%jz0(i,j,k)
   kx=B%kx0(i,j,k); ky=B%ky0(i,j,k); kz=B%kz0(i,j,k)

   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
 
   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky

   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz

! 涡量
   S=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
   X=d(i,j,k)*vt(i,j,k)/B%mu(i,j,k)   ! 湍流粘性系数与层流粘性系数之比
!--------------------------------------------------------------------------   

! 源项 采用 Blazek's book,  p241 (7.38), (7.39)


!   Blazek's book 的公式稳定性不好 (容易算出负的湍流粘性系数)
!   Source term, Blazek's Book section (7.2.1),  modified from original form
!   fv1=X**3/(X**3+Cv1**3)
!   fv2=1.d0/(1.d0+X/Cv2)**3
!   fv3=(1.d0+X*fv1)*(1.d0-fv2)/max(X,0.001d0)
!   S1=fv3*S+vt(i,j,k)*fv2/(SA_k*B%dw(i,j,k))**2           
!   r=vt(i,j,k)/(S1*SA_K*SA_K*B%dw(i,j,k)*B%dw(i,j,k))

!------------------------------------------------------------------------
! Source term, original form; See: http://turbmodels.larc.nasa.gov/spalart.html
   fv1=X**3/(X**3+Cv1**3)
   fv2=1.d0-X/(1.d0+X*fv1)
   ft2=Ct3*exp(-Ct4*X*X)            
   S1=max(S+fv2*vt(i,j,k)/(Re*(SA_k*B%dw(i,j,k))**2),0.d0)    !!! Re
   r=min(vt(i,j,k)/(Re*S1*SA_K*SA_K*B%dw(i,j,k)*B%dw(i,j,k)),10.0_PRE_EC)    !!! Re
   g=r+Cw2*(r**6-r)
   fw=g*((1.d0+Cw3**6)/(g**6+Cw3**6))**(1.d0/6.d0)
!   Q_SA=Cb1*(1.d0-ft2)*S1*vt(i,j,k) -(Cw1*fw-Cb1*ft2/(SA_k*SA_k))*(vt(i,j,k)/B%dw(i,j,k))**2/Re   !!! Re
!---------Source term, see: CFL3D manual ---------------------
    Q3=Cb1*((1-ft2)*fv2+ft2)/(SA_k*SA_k)-Cw1*fw
	Q_SA=Cb1*(1.d0-ft2)*S1*vt(i,j,k)   &
	    +Q3*(vt(i,j,k)/B%dw(i,j,k))**2/Re   !!! Re


!---------------------------------------------------------------------------
   B%Res(6,i,j,k)=B%Res(6,i,j,k)+Q_SA*B%vol(i,j,k)     ! Bug removed

   enddo
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Re)  
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)/Re             ! 量纲转换
    B%mu_t(i,j,k)=B%mu_t(i,j,k)/Re
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO   



  deallocate(vt,Fluxv,fluxv2)

end  subroutine  Turbulence_model_SA

!-----------------------------------------------------------------------------------
!-----根据vt (v~) 计算出SA模型中的湍流粘性系数mut-----------------------------------
! Blazek's Book p241 (7.37)



! 粘性系数虚网格上的值 （固壁边界采用反值，以保证固壁上的平均湍流粘性系数为0）
subroutine Amut_boundary(nMesh,mBlock)
   use Global_Var
   use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,m,ksub,nMesh
   integer:: ib,ie,jb,je,kb,ke
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc

   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz

! mut in Ghost Cell of the boundary
! 采用临近点的值

!   B%Amu_t(0,:,:)=B%Amu_t(1,:,:)
!   B%Amu_t(nx,:,:)=B%Amu_t(nx-1,:,:)
!   B%Amu_t(:,0,:)=B%Amu_t(:,1,:)
!   B%Amu_t(:,ny,:)=B%Amu_t(:,ny-1,:)
!   B%Amu_t(:,:,0)=B%Amu_t(:,:,1)
!   B%Amu_t(:,:,nz)=B%Amu_t(:,:,nz-1)

! Ghost Cell 点的 mut值为 内点mut值*（-1） (这样可以使壁面上mut=0)

  do  ksub=1,B%subface
    Bc=> B%bc_msg(ksub)
      if( Bc%bc .eq. BC_Wall  ) then   ! (粘性) 壁面边界条件 
      Bc => B%bc_msg(ksub)
      ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)  
     if(Bc%face .eq. 1 ) then   ! i- 
!$OMP DO
       do k=kb,ke-1
       do j=jb,je-1
         B%mu_t(0,j,k)= -B%mu_t(1,j,k)       ! mut
       enddo
       enddo
!$OMP ENDDO
     else if (Bc%face .eq. 4 ) then   ! i+
!$OMP DO
       do k=kb,ke-1
       do j=jb,je-1
         B%mu_t(ie,j,k)= -B%mu_t(ie-1,j,k)       ! mut
       enddo
       enddo
!$OMP ENDDO
     else if(Bc%face .eq. 2 ) then   !j-
!$OMP DO
       do k=kb,ke-1
       do i=ib,ie-1
         B%mu_t(i,0,k)= -B%mu_t(i,1,k)       ! mut
      enddo
      enddo
!$OMP ENDDO
     else if(Bc%face .eq. 5 ) then   !j+
!$OMP DO
       do k=kb,ke-1
       do i=ib,ie-1
         B%mu_t(i,je,k)= -B%mu_t(i,je-1,k)       ! mut
      enddo
      enddo
!$OMP ENDDO
     else if(Bc%face .eq. 3 ) then   !k-
!$OMP DO
       do j=jb,je-1
       do i=ib,ie-1
         B%mu_t(i,j,0)= -B%mu_t(i,j,1) 
       enddo
       enddo
!$OMP ENDDO
     else if(Bc%face .eq. 6 ) then   !k+
!$OMP DO
       do j=jb,je-1
       do i=ib,ie-1
         B%mu_t(i,j,ke)= -B%mu_t(i,j,ke-1)
       enddo
       enddo
!$OMP ENDDO

      endif
!$OMP END PARALLEL

    endif
   enddo

end  subroutine Amut_boundary

!------------------------------------------------------------------------------
! Source term of SA model  (See: J. Blazek's book, P240-243) 
! Do not consider transition (full turbulence)
! Ver 0.81a, using Eq. (7.4.2) (J. Blazek's book)
! Ver 0.81b, limit for sorce term  ( not less than 0.3*Omega)
! Ver 0.98c, Nondimensional (See CFL3D manual)

  subroutine  Turbulence_model_NewSA(nMesh,mBlock)
   use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,nMesh
   real(PRE_EC):: ui,vi,wi,uj,vj,wj,uk,vk,wk,ux,vx,wx,uy,vy,wy,uz,vz,wz
   real(PRE_EC):: s1x,s1y,s1z,v0,vti,vtj,vtk,vtx,vty,vtz
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz
   real(PRE_EC):: S, S1,X,fv1,fv2,fv3,ft1,ft2,r,g,fw,Q_SA,vn1,vn2,vfi,Q3
   real(PRE_EC):: ppi,pj,pk,px,py,pz,p_plus,Cb1,Cw1
  
   real(PRE_EC),parameter:: SA_sigma=2.d0/3.d0,Cv1=7.1d0,Cv2=5.d0,Cb2=0.622d0,SA_k=0.41d0
   real(PRE_EC),parameter:: Cw2=0.3d0, Cw3=2.d0, Ct1=1.d0, Ct2=2.d0, Ct3=1.2d0, Ct4=0.5d0  !Ct3=1.3d0
!   real(PRE_EC),parameter:: Cb1=0.1355d0, Cw1=Cb1/(SA_k*SA_k)+(1.d0+Cb2)/SA_sigma
   real(PRE_EC),parameter:: ep=1.d-6
   real(PRE_EC),Pointer,dimension(:,:,:):: vt,Fluxv,fluxv2
   Type (Block_TYPE),pointer:: B
 
  
   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz



 ! 计算湍流粘性系数
   allocate(vt(0:nx,0:ny,0:nz),Fluxv(nx,ny,nz),fluxv2(nx,ny,nz))
	
! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量



!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Re)  
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)*Re             ! 量纲转换
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO   



!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,vt,d)
 
    do k=0,nz
    do j=0,ny
    do i=0,nx
     vt(i,j,k)=B%U(6,i,j,k)
	 X=d(i,j,k)*vt(i,j,k)/B%mu(i,j,k)   ! 湍流粘性系数与层流粘性系数之比
     fv1=X**3/(X**3+Cv1**3)
     B%mu_t(i,j,k)=fv1*d(i,j,k)*vt(i,j,k)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
   
   ! 限定湍流粘性系数
   call limit_mut(nMesh,mBlock)

! 设定湍流粘性系数虚网格的值
 

! 计算vt方程的残差 B%Res(6,:,:,:)

!$OMP PARALLEL  DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,Re,vt,d,uu,v,w,Fluxv,fluxv2,CP1_NSA,CP2_NSA)
!------i- direcion ---------------------------
!$OMP DO
    do k=1,nz-1 
     do j=1,ny-1
      do i=1,nx
          s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)  ! 归一化的法方向
          vti=vt(i,j,k)-vt(i-1,j,k)                               ! SA模型中的vt
          vtj=0.25d0*(vt(i,j+1,k)-vt(i,j-1,k)+vt(i-1,j+1,k)-vt(i-1,j-1,k))
          vtk=0.25d0*(vt(i,j,k+1)-vt(i,j,k-1)+vt(i-1,j,k+1)-vt(i-1,j,k-1))
          ix=B%ix1(i,j,k); iy=B%iy1(i,j,k); iz=B%iz1(i,j,k)
          jx=B%jx1(i,j,k); jy=B%jy1(i,j,k); jz=B%jz1(i,j,k)
          kx=B%kx1(i,j,k); ky=B%ky1(i,j,k); kz=B%kz1(i,j,k)
          vtx=vti*ix+vtj*jx+vtk*kx
          vty=vti*iy+vtj*jy+vtk*ky
          vtz=vti*iz+vtj*jz+vtk*kz
          v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i-1,j,k)+B%mu(i-1,j,k)/d(i-1,j,k))  ! (I-1/2,J,K) 点的动力学粘性系数
 
          vn1=uu(i-1,j,k)*s1x+v(i-1,j,k)*s1y+w(i-1,j,k)*s1z   ! 法向速度
          vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
          vfi=0.5d0*((vn1+abs(vn1))*vt(i-1,j,k)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
          Fluxv(i,j,k)= (-vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))* B%Si(i,j,k)    !!! Re
          Fluxv2(i,j,k)=(vtx*s1x+vty*s1y+vtz*s1z)*B%Si(i,j,k)
	  enddo
     enddo
    enddo
!$OMP END DO
!$OMP DO
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
	     v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
         B%Res(6,i,j,k)=Fluxv(i+1,j,k)-Fluxv(i,j,k) -v0/Re*(Fluxv2(i+1,j,k)-Fluxv2(i,j,k))    !!! Re     
      enddo
      enddo
      enddo
!$OMP END DO
!-----j- direction ------------------------
!$OMP DO
    do k=1,nz-1 
    do j=1,ny
    do i=1,nx-1
      s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向
      vti=0.25d0*(vt(i+1,j,k)-vt(i-1,j,k)+vt(i+1,j-1,k)-vt(i-1,j-1,k))
      vtj=vt(i,j,k)-vt(i,j-1,k)
      vtk=0.25d0*(vt(i,j,k+1)-vt(i,j,k-1)+vt(i,j-1,k+1)-vt(i,j-1,k-1))
      ix=B%ix2(i,j,k); iy=B%iy2(i,j,k); iz=B%iz2(i,j,k)
      jx=B%jx2(i,j,k); jy=B%jy2(i,j,k); jz=B%jz2(i,j,k)
      kx=B%kx2(i,j,k); ky=B%ky2(i,j,k); kz=B%kz2(i,j,k)
      vtx=vti*ix+vtj*jx+vtk*kx
      vty=vti*iy+vtj*jy+vtk*ky
      vtz=vti*iz+vtj*jz+vtk*kz
      v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i,j-1,k)+B%mu(i,j-1,k)/d(i,j-1,k))  ! (I-1/2,J,K) 点的动力学粘性系数   s11=Amu1*(tmp1*ux-tmp2*(vy+wz))    ! tmp1=4.d0/3.d0; tmp2=2.d0/3.d0
      vn1=uu(i,j-1,k)*s1x+v(i,j-1,k)*s1y+w(i,j-1,k)*s1z   ! 法向速度
      vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
      vfi=0.5d0*((vn1+abs(vn1))*vt(i,j-1,k)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
      Fluxv(i,j,k)= (-vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))* B%Sj(i,j,k)   !!! Re
      Fluxv2(i,j,k)= (vtx*s1x+vty*s1y+vtz*s1z)* B%Sj(i,j,k)
    enddo
    enddo
    enddo
!$OMP END DO
!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
	  v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
      B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxv(i,j+1,k)-Fluxv(i,j,k)- v0/Re*(Fluxv2(i,j+1,k)-Fluxv2(i,j,k))   !!! Re        
     enddo
     enddo
     enddo
!$OMP END DO
!-----k- direction ------------------------
!$OMP DO
   do k=1,nz 
   do j=1,ny-1
   do i=1,nx-1
    s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向
    vti=0.25d0*(vt(i+1,j,k)-vt(i-1,j,k)+vt(i+1,j,k-1)-vt(i-1,j,k-1))
    vtj=0.25d0*(vt(i,j+1,k)-vt(i,j-1,k)+vt(i,j+1,k-1)-vt(i,j-1,k-1))
    vtk=vt(i,j,k)-vt(i,j,k-1)

      ix=B%ix3(i,j,k); iy=B%iy3(i,j,k); iz=B%iz3(i,j,k)
      jx=B%jx3(i,j,k); jy=B%jy3(i,j,k); jz=B%jz3(i,j,k)
      kx=B%kx3(i,j,k); ky=B%ky3(i,j,k); kz=B%kz3(i,j,k)
      vtx=vti*ix+vtj*jx+vtk*kx
      vty=vti*iy+vtj*jy+vtk*ky
      vtz=vti*iz+vtj*jz+vtk*kz
      v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k) + vt(i,j,k-1)+B%mu(i,j,k-1)/d(i,j,k-1))  ! (I,J,K-1/2) 点的动力学粘性系数
     vn1=uu(i,j,k-1)*s1x+v(i,j,k-1)*s1y+w(i,j,k-1)*s1z   ! 法向速度
     vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
     vfi=0.5d0*((vn1+abs(vn1))*vt(i,j,k-1)+(vn2-abs(vn2))*vt(i,j,k))  ! 一阶 L-F格式
     Fluxv(i,j,k)= ( -vfi+v0/Re*(vtx*s1x+vty*s1y+vtz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
     Fluxv2(i,j,k)= (vtx*s1x+vty*s1y+vtz*s1z)*B%Sk(i,j,k)   
    enddo
    enddo
    enddo
!$OMP END DO
!$OMP DO
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
 	   v0=(vt(i,j,k)+B%mu(i,j,k)/d(i,j,k))*Cb2/SA_sigma
       B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxv(i,j,k+1)-Fluxv(i,j,k)-v0/Re*(Fluxv2(i,j,k+1)-Fluxv2(i,j,k))           
	 enddo
     enddo
     enddo
!$OMP END DO
!--------------源项---------------------------------------------------------
!$OMP DO
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1
!----- get S (normal of vorticity)  S=sqrt(0.5*Omiga_ij*Omiga_ij) at the cell's center ------
!  计算涡量；计算湍流粘性系数的梯度
      

   ui=uu(i+1,j,k)-uu(i-1,j,k)            
   vi=v(i+1,j,k)-v(i-1,j,k)  
   wi=w(i+1,j,k)-w(i-1,j,k)  

   uj=uu(i,j+1,k)-uu(i,j-1,k)   
   vj=v(i,j+1,k)-v(i,j-1,k)
   wj=w(i,j+1,k)-w(i,j-1,k) 

   uk=uu(i,j,k+1)-uu(i,j,k-1)  
   vk=v(i,j,k+1)-v(i,j,k-1)
   wk=w(i,j,k+1)-w(i,j,k-1)  
  

  
   
   ix=B%ix0(i,j,k); iy=B%iy0(i,j,k); iz=B%iz0(i,j,k)
   jx=B%jx0(i,j,k); jy=B%jy0(i,j,k); jz=B%jz0(i,j,k)
   kx=B%kx0(i,j,k); ky=B%ky0(i,j,k); kz=B%kz0(i,j,k)

   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
 
   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky

   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz


   ppi=p(i+1,j,k)-p(i-1,j,k)  
   pj=p(i,j+1,k)-p(i,j-1,k) 
   pk=w(i,j,k+1)-p(i,j,k-1)  
   px=ppi*ix+pj*jx+pk*kx
   py=ppi*iy+pj*jy+pk*ky
   pz=ppi*iz+pj*jz+pk*kz


! 涡量
   S=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
   X=d(i,j,k)*vt(i,j,k)/B%mu(i,j,k)   ! 湍流粘性系数与层流粘性系数之比
!--------------------------------------------------------------------------   

! 源项 采用 Blazek's book,  p241 (7.38), (7.39)


!   Blazek's book 的公式稳定性不好 (容易算出负的湍流粘性系数)
!   Source term, Blazek's Book section (7.2.1),  modified from original form
!   fv1=X**3/(X**3+Cv1**3)
!   fv2=1.d0/(1.d0+X/Cv2)**3
!   fv3=(1.d0+X*fv1)*(1.d0-fv2)/max(X,0.001d0)
!   S1=fv3*S+vt(i,j,k)*fv2/(SA_k*B%dw(i,j,k))**2           
!   r=vt(i,j,k)/(S1*SA_K*SA_K*B%dw(i,j,k)*B%dw(i,j,k))

!------------------------------------------------------------------------
! Source term, original form; See: http://turbmodels.larc.nasa.gov/spalart.html
   fv1=X**3/(X**3+Cv1**3)
   fv2=1.d0-X/(1.d0+X*fv1)
   ft2=Ct3*exp(-Ct4*X*X)            
   S1=max(S+fv2*vt(i,j,k)/(Re*(SA_k*B%dw(i,j,k))**2),0.d0)    !!! Re
   r=min(vt(i,j,k)/(Re*S1*SA_K*SA_K*B%dw(i,j,k)*B%dw(i,j,k)),10.0_PRE_EC)    !!! Re
   g=r+Cw2*(r**6-r)
   fw=g*((1.d0+Cw3**6)/(g**6+Cw3**6))**(1.d0/6.d0)
!---------Source term, see: CFL3D manual ---------------------
   p_plus=B%mu(i,j,k)/(d(i,j,k)**3*(B%mu(i,j,k)/d(i,j,k)*S)**1.5d0+ep)*sqrt(px*px+py*py+pz*pz)
 
    Cb1=0.1355d0*(1.d0 +Cp1_NSA*(1.d0-exp(-p_plus/Cp2_NSA)) )
    Cw1=0.1355d0/(SA_k*SA_k)+(1.d0+Cb2)/SA_sigma
!   Cw1=Cb1/(SA_k*SA_k)*exp(-P_plus/Cp2)+(1.d0+Cb2)/SA_sigma



    Q3=Cb1*((1-ft2)*fv2+ft2)/(SA_k*SA_k)-Cw1*fw
	Q_SA=Cb1*(1.d0-ft2)*S1*vt(i,j,k)   &
	    +Q3*(vt(i,j,k)/B%dw(i,j,k))**2/Re   !!! Re


!---------------------------------------------------------------------------
   B%Res(6,i,j,k)=B%Res(6,i,j,k)+Q_SA*B%vol(i,j,k)     ! Bug removed

   enddo
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Re)  
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)/Re             ! 量纲转换
    B%mu_t(i,j,k)=B%mu_t(i,j,k)/Re
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO   



  deallocate(vt,Fluxv,fluxv2)

end  subroutine  Turbulence_model_NewSA

!-----------------------------------------------------------------------------------
!-----根据vt (v~) 计算出SA模型中的湍流粘性系数mut-----------------------------------
! Blazek's Book p241 (7.37)

!------------------------------------------------------------------------------
! SST model  (See: J. Blazek's book, section 7.2.3) 
! Do not consider transition (full turbulence)
! See:  CFL3D 5.0 manual
! 2018-9-29:  A bug is removed;

  subroutine  Turbulence_model_SST(nMesh,mBlock)
   use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,nMesh
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz
   real(PRE_EC):: s1x,s1y,s1z
   real(PRE_EC):: ui,vi,wi,Kti,Wti,uj,vj,wj,Ktj,Wtj,uk,vk,wk,Ktk,Wtk,  &
                  ux,vx,wx,Ktx,Wtx,uy,vy,wy,Kty,Wty,uz,vz,wz,Ktz,Wtz
   real(PRE_EC):: omega,arg1,arg2,arg3,f2,Kws,CD_kw,Pk,Pk1,Pk0,Qk,Qw
   real(PRE_EC):: t11,t22,t33,t12,t13,t23,muk,muw,vn1,vn2,kfi,wfi
   real(PRE_EC):: Cw_SST,beta_SST,sigma_k_SST,sigma_w_SST
   
   real(PRE_EC),parameter:: sigma_k1_SST=0.85d0,sigma_w1_SST=0.5d0,beta1_SST=0.075d0,Cw1_SST=0.533d0, &
                        sigma_k2_SST=1.d0,  sigma_w2_SST=0.856d0,beta2_SST=0.0828d0,Cw2_SST=0.440d0
   real(PRE_EC),parameter::a1_SST=0.31d0, betas_SST=0.09d0
   real(PRE_EC),Pointer,dimension(:,:,:):: Kt,Wt,Fluxk,Fluxw  ! 湍能、湍能比耗散率；源项；通量项
   real(PRE_EC),Pointer,dimension(:,:,:):: f1

   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
 
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz
 ! 计算湍流粘性系数
   allocate(Kt(0:nx,0:ny,0:nz),Wt(0:nx,0:ny,0:nz),Fluxk(nx,ny,nz),Fluxw(nx,ny,nz))
   allocate(f1(nx,ny,nz))

! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Kt,Wt,Fluxk,Fluxw,d,uu,v,w,f1)
    
!$OMP DO   
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)*Re             ! 量纲转换
   enddo
   enddo
   enddo
!$OMP END DO   

!$OMP DO   
    do k=0,nz
	do j=0,ny
	do i=0,nx
	 Kt(i,j,k)=B%U(6,i,j,k)/B%U(1,i,j,k)
	 Wt(i,j,k)=B%U(7,i,j,k)/B%U(1,i,j,k)
    enddo
	enddo
	enddo
!$OMP END DO   
   

!$OMP DO   
    do k=1,nz-1
    do j=1,ny-1
    do i=1,nx-1
! 计算涡量, 交叉对流项
    ui=uu(i+1,j,k)-uu(i-1,j,k)            
    vi=v(i+1,j,k)-v(i-1,j,k)  
    wi=w(i+1,j,k)-w(i-1,j,k)  
    Kti=Kt(i+1,j,k)-Kt(i-1,j,k)  
    wti=wt(i+1,j,k)-wt(i-1,j,k)  
 
    uj=uu(i,j+1,k)-uu(i,j-1,k)   
    vj=v(i,j+1,k)-v(i,j-1,k)
    wj=w(i,j+1,k)-w(i,j-1,k) 
    Ktj=Kt(i,j+1,k)-Kt(i,j-1,k) 
    wtj=wt(i,j+1,k)-wt(i,j-1,k) 
    
	uk=uu(i,j,k+1)-uu(i,j,k-1)  
    vk=v(i,j,k+1)-v(i,j,k-1)
    wk=w(i,j,k+1)-w(i,j,k-1)  
    Ktk=Kt(i,j,k+1)-Kt(i,j,k-1)  
    wtk=wt(i,j,k+1)-wt(i,j,k-1) 
	 
    ix=B%ix0(i,j,k); iy=B%iy0(i,j,k); iz=B%iz0(i,j,k)
    jx=B%jx0(i,j,k); jy=B%jy0(i,j,k); jz=B%jz0(i,j,k)
    kx=B%kx0(i,j,k); ky=B%ky0(i,j,k); kz=B%kz0(i,j,k)
   
    ux=ui*ix+uj*jx+uk*kx
    vx=vi*ix+vj*jx+vk*kx
    wx=wi*ix+wj*jx+wk*kx
    Ktx=kti*ix+ktj*jx+ktk*kx
    wtx=wti*ix+wtj*jx+wtk*kx

    uy=ui*iy+uj*jy+uk*ky
    vy=vi*iy+vj*jy+vk*ky
    wy=wi*iy+wj*jy+wk*ky
    Kty=kti*iy+ktj*jy+ktk*ky
    wty=wti*iy+wtj*jy+wtk*ky

    uz=ui*iz+uj*jz+uk*kz
    vz=vi*iz+vj*jz+vk*kz
    wz=wi*iz+wj*jz+wk*kz
    Ktz=kti*iz+ktj*jz+ktk*kz
    wtz=wti*iz+wtj*jz+wtk*kz

! 涡量
     omega=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
	 arg2=max( 2.d0* sqrt(abs(Kt(i,j,k)))/(0.09*Wt(i,j,k)*B%dw(i,j,k)*Re) , &
	          500.d0*B%mu(i,j,k)/(d(i,j,k)*Wt(i,j,k)*B%dw(i,j,k)**2 *Re*Re) )
     f2=tanh(arg2**2)
     B%mu_t(i,j,k)=a1_SST*d(i,j,k)*Kt(i,j,k)/max(a1_SST*Wt(i,j,k),f2*abs(omega)/Re)


 ! 计算f1 (识别是否为近壁区，近壁区趋近于1）      
     
 !    Kws=2.d0*(ktx*wtx+kty*wty+ktz*ktz)*d(i,j,k)*sigma_w2_SST/(Wt(i,j,k)+1.d-20)      ! Bug
      Kws=2.d0*(ktx*wtx+kty*wty+ktz*wtz)*d(i,j,k)*sigma_w2_SST/(Wt(i,j,k)+1.d-20)      ! 交叉输运项
    
     CD_kw=max(Kws,1.d-20)
     arg3=max(sqrt(abs(Kt(i,j,k)))/(0.09*Wt(i,j,k)*B%dw(i,j,k) *Re)  , &
	          500.d0*B%mu(i,j,k)/(d(i,j,k)*Wt(i,j,k)*B%dw(i,j,k)**2 *Re*Re) )
	 arg1=min(arg3,4.d0*d(i,j,k)*sigma_w2_SST*Kt(i,j,k)/(CD_kw*B%dw(i,j,k)**2 ))
     f1(i,j,k)=tanh(arg1**4)             ! 开关函数，近壁区趋近于1，远壁区趋近于0  （用来切换k-w及k-epsl方程)
     
     
!    湍应力 （使用了涡粘模型）     ! Blazek's Book, Eq. (7.25)
        
!	     t11=(4.d0/3.d0)*ux-(2.d0/3.d0)*(vy+wz)  
!         t22=(4.d0/3.d0)*vy-(2.d0/3.d0)*(ux+wz) 
!         t33=(4.d0/3.d0)*wz-(2.d0/3.d0)*(ux+vy) 
!         t12=uy+vx
!         t13=uz+wx
!         t23=vz+wy


!    湍能方程的源项（生成-耗散)     
!     Pk1=t11*ux+t22*vy+t33*wz+t12*(uy+vx)+t13*(uz+wx)+t23*(vz+wy)     
!	  Pk=B%mu_t(i,j,k)*Pk1                                                   ! 湍能生成项 （湍应力乘以应变率）       

      Pk=B%mu_t(i,j,k)*omega*omega
!     Pk0=min(Pk,20.d0*betas_SST*Kt(i,j,k)*Wt(i,j,k)*Re*Re)                        ! 对湍能生成项进行限制，防止湍能过大
   
     Pk0=Pk        ! 不进行限制

	 Qk=Pk0/Re-Re*betas_SST*d(i,j,k)*Wt(i,j,k)*Kt(i,j,k)    ! k方程的源项  （生成项-耗散项）

     Cw_SST=f1(i,j,k)*Cw1_SST+(1.d0-f1(i,j,k))*Cw2_SST    ! 模型系数，利用f1函数进行切换
     beta_SST=f1(i,j,k)*beta1_SST+(1.d0-f1(i,j,k))*beta2_SST    ! 模型系数，利用f1函数进行切换
!    Qw= Cw_SST*d(i,j,k)*Pk1          &
!	           -beta_SST*d(i,j,k)*Wt(i,j,k)**2+(1.d0-f1(i,j,k))*Kws     ! W方程的源项    
     Qw= Cw_SST*d(i,j,k)*omega*omega/Re          &
	           -Re*beta_SST*d(i,j,k)*Wt(i,j,k)**2+(1.d0-f1(i,j,k))*Kws/Re     ! W方程的源项    

!-------------------------------------------

	B%Res(6,i,j,k)=QK*B%vol(i,j,k)
	B%Res(7,i,j,k)=Qw*B%vol(i,j,k)
 
	enddo
	enddo
    enddo
!$OMP END DO   
!$OMP END PARALLEL
	

! 设定湍流粘性系数虚网格的值
! mut in Ghost Cell of the boundary
! 采用临近点的值

   B%mu_t(0,:,:)=B%mu_t(1,:,:)
   B%mu_t(nx,:,:)=B%mu_t(nx-1,:,:)
   B%mu_t(:,0,:)=B%mu_t(:,1,:)
   B%mu_t(:,ny,:)=B%mu_t(:,ny-1,:)
   B%mu_t(:,:,0)=B%mu_t(:,:,1)
   B%mu_t(:,:,nz)=B%mu_t(:,:,nz-1)
!  固壁上
   call  Amut_boundary(nMesh,mBlock)
 
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,B,Kt,Wt,Fluxk,Fluxw,d,uu,v,w,f1)
!$OMP DO   
!------i- direcion ---------------------------
     do k=1,nz-1 
     do j=1,ny-1
     do i=1,nx
! 扩散项
! 扩散系数，界面上的值=两侧值的平均, 边界上的扩散系数=内侧的值   
       sigma_K_SST=f1(i,j,k)*sigma_k1_SST+(1.d0-f1(i,j,k))*sigma_k2_SST
       sigma_W_SST=f1(i,j,k)*sigma_w1_SST+(1.d0-f1(i,j,k))*sigma_w2_SST
       muk=(B%mu(i-1,j,k)+B%mu(i,j,k) + sigma_K_SST*(B%mu_t(i-1,j,k)+B%mu_t(i,j,k)) )*0.5d0 /Re        ! 扩散系数 (k方程), 界面上的值=两侧值的平均
       muw=(B%mu(i-1,j,k)+B%mu(i,j,k) + sigma_W_SST*(B%mu_t(i-1,j,k)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (w方程)
       s1x=B%ni1(i,j,k); s1y=B%ni2(i,j,k) ; s1z= B%ni3(i,j,k)  ! 归一化的法方向
          Kti=Kt(i,j,k)-Kt(i-1,j,k)                               ! K
          Wti=wt(i,j,k)-wt(i-1,j,k)                               ! W
          Ktj=0.25d0*(Kt(i,j+1,k)-Kt(i,j-1,k)+Kt(i-1,j+1,k)-Kt(i-1,j-1,k))
          Wtj=0.25d0*(Wt(i,j+1,k)-Wt(i,j-1,k)+Wt(i-1,j+1,k)-Wt(i-1,j-1,k))
          Ktk=0.25d0*(Kt(i,j,k+1)-Kt(i,j,k-1)+Kt(i-1,j,k+1)-Kt(i-1,j,k-1))
          Wtk=0.25d0*(Wt(i,j,k+1)-Wt(i,j,k-1)+Wt(i-1,j,k+1)-wt(i-1,j,k-1))

          ix=B%ix1(i,j,k); iy=B%iy1(i,j,k); iz=B%iz1(i,j,k)
          jx=B%jx1(i,j,k); jy=B%jy1(i,j,k); jz=B%jz1(i,j,k)
          kx=B%kx1(i,j,k); ky=B%ky1(i,j,k); kz=B%kz1(i,j,k)
          ktx=kti*ix+ktj*jx+ktk*kx                                ! Kt对坐标的导数
          wtx=wti*ix+wtj*jx+wtk*kx
          kty=kti*iy+ktj*jy+ktk*ky
          wty=wti*iy+wtj*jy+wtk*ky
          ktz=kti*iz+ktj*jz+ktk*kz
          wtz=wti*iz+wtj*jz+wtk*kz
!          对流项
          vn1=uu(i-1,j,k)*s1x+v(i-1,j,k)*s1y+w(i-1,j,k)*s1z   ! 法向速度
          vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
          kfi=0.5d0*((vn1+abs(vn1))*kt(i-1,j,k)+(vn2-abs(vn2))*kt(i,j,k))  ! 对流项，一阶 L-F格式
          wfi=0.5d0*((vn1+abs(vn1))*wt(i-1,j,k)+(vn2-abs(vn2))*wt(i,j,k))  ! 对流项，一阶 L-F格式
		  Fluxk(i,j,k)= (-kfi+muk*(ktx*s1x+kty*s1y+ktz*s1z))* B%Si(i,j,k)
		  Fluxw(i,j,k)= (-wfi+muw*(wtx*s1x+wty*s1y+wtz*s1z))* B%Si(i,j,k)
      enddo
     enddo
    enddo
!$OMP END DO   

!$OMP DO   
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
         B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxk(i+1,j,k)-Fluxk(i,j,k)        
         B%Res(7,i,j,k)=B%Res(7,i,j,k)+Fluxw(i+1,j,k)-Fluxw(i,j,k)        
	  enddo
      enddo
      enddo
!$OMP END DO   

!-----j- direction ------------------------

!$OMP  DO   
    do k=1,nz-1 
    do j=1,ny
    do i=1,nx-1
       sigma_K_SST=f1(i,j,k)*sigma_k1_SST+(1.d0-f1(i,j,k))*sigma_k2_SST
       sigma_W_SST=f1(i,j,k)*sigma_w1_SST+(1.d0-f1(i,j,k))*sigma_w2_SST
       muk=(B%mu(i,j-1,k)+B%mu(i,j,k) + sigma_K_SST*(B%mu_t(i,j-1,k)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (k方程), 界面上的值=两侧值的平均
       muw=(B%mu(i,j-1,k)+B%mu(i,j,k) + sigma_W_SST*(B%mu_t(i,j-1,k)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (w方程)

      s1x=B%nj1(i,j,k); s1y=B%nj2(i,j,k) ; s1z= B%nj3(i,j,k)  ! 归一化的法方向
      kti=0.25d0*(kt(i+1,j,k)-kt(i-1,j,k)+kt(i+1,j-1,k)-kt(i-1,j-1,k))
      ktj=kt(i,j,k)-kt(i,j-1,k)
      ktk=0.25d0*(kt(i,j,k+1)-kt(i,j,k-1)+kt(i,j-1,k+1)-kt(i,j-1,k-1))
      wti=0.25d0*(wt(i+1,j,k)-wt(i-1,j,k)+wt(i+1,j-1,k)-wt(i-1,j-1,k))
      wtj=wt(i,j,k)-wt(i,j-1,k)
      wtk=0.25d0*(wt(i,j,k+1)-wt(i,j,k-1)+wt(i,j-1,k+1)-wt(i,j-1,k-1))
     
	  ix=B%ix2(i,j,k); iy=B%iy2(i,j,k); iz=B%iz2(i,j,k)
      jx=B%jx2(i,j,k); jy=B%jy2(i,j,k); jz=B%jz2(i,j,k)
      kx=B%kx2(i,j,k); ky=B%ky2(i,j,k); kz=B%kz2(i,j,k)
      ktx=kti*ix+ktj*jx+ktk*kx
      kty=kti*iy+ktj*jy+ktk*ky
      ktz=kti*iz+ktj*jz+ktk*kz
      wtx=wti*ix+wtj*jx+wtk*kx
      wty=wti*iy+wtj*jy+wtk*ky
      wtz=wti*iz+wtj*jz+wtk*kz
      vn1=uu(i,j-1,k)*s1x+v(i,j-1,k)*s1y+w(i,j-1,k)*s1z   ! 法向速度
      vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
      kfi=0.5d0*((vn1+abs(vn1))*kt(i,j-1,k)+(vn2-abs(vn2))*kt(i,j,k))  ! 一阶 L-F格式
      wfi=0.5d0*((vn1+abs(vn1))*wt(i,j-1,k)+(vn2-abs(vn2))*wt(i,j,k))  ! 一阶 L-F格式
      
	  Fluxk(i,j,k)= (-kfi+muk*(ktx*s1x+kty*s1y+ktz*s1z))* B%Sj(i,j,k)
	  Fluxw(i,j,k)= (-wfi+muw*(wtx*s1x+wty*s1y+wtz*s1z))* B%Sj(i,j,k)

    enddo
    enddo
    enddo
!$OMP END DO 
!$OMP  DO 
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
      B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxk(i,j+1,k)-Fluxk(i,j,k)           
      B%Res(7,i,j,k)=B%Res(7,i,j,k)+Fluxw(i,j+1,k)-Fluxw(i,j,k)           
	 enddo
     enddo
     enddo
!$OMP END DO 
!-----k- direction ------------------------
!$OMP  DO 
   do k=1,nz 
   do j=1,ny-1
   do i=1,nx-1
      sigma_K_SST=f1(i,j,k)*sigma_k1_SST+(1.d0-f1(i,j,k))*sigma_k2_SST
      sigma_W_SST=f1(i,j,k)*sigma_w1_SST+(1.d0-f1(i,j,k))*sigma_w2_SST
      muk=(B%mu(i,j,k-1)+B%mu(i,j,k) + sigma_K_SST*(B%mu_t(i,j,k-1)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (k方程), 界面上的值=两侧值的平均
      muw=(B%mu(i,j,k-1)+B%mu(i,j,k) + sigma_W_SST*(B%mu_t(i,j,k-1)+B%mu_t(i,j,k)) )*0.5d0 /Re       ! 扩散系数 (w方程)

      s1x=B%nk1(i,j,k); s1y=B%nk2(i,j,k) ; s1z= B%nk3(i,j,k)  ! 归一化的法方向
      kti=0.25d0*(kt(i+1,j,k)-kt(i-1,j,k)+kt(i+1,j,k-1)-kt(i-1,j,k-1))
      ktj=0.25d0*(kt(i,j+1,k)-kt(i,j-1,k)+kt(i,j+1,k-1)-kt(i,j-1,k-1))
      ktk=kt(i,j,k)-kt(i,j,k-1)
      wti=0.25d0*(wt(i+1,j,k)-wt(i-1,j,k)+wt(i+1,j,k-1)-wt(i-1,j,k-1))
      wtj=0.25d0*(wt(i,j+1,k)-wt(i,j-1,k)+wt(i,j+1,k-1)-wt(i,j-1,k-1))
      wtk=wt(i,j,k)-wt(i,j,k-1)
      ix=B%ix3(i,j,k); iy=B%iy3(i,j,k); iz=B%iz3(i,j,k)
      jx=B%jx3(i,j,k); jy=B%jy3(i,j,k); jz=B%jz3(i,j,k)
      kx=B%kx3(i,j,k); ky=B%ky3(i,j,k); kz=B%kz3(i,j,k)
      ktx=kti*ix+ktj*jx+ktk*kx
      kty=kti*iy+ktj*jy+ktk*ky
      ktz=kti*iz+ktj*jz+ktk*kz
      wtx=wti*ix+wtj*jx+wtk*kx
      wty=wti*iy+wtj*jy+wtk*ky
      wtz=wti*iz+wtj*jz+wtk*kz

     vn1=uu(i,j,k-1)*s1x+v(i,j,k-1)*s1y+w(i,j,k-1)*s1z   ! 法向速度
     vn2=uu(i,j,k)*s1x+v(i,j,k)*s1y+w(i,j,k)*s1z
     kfi=0.5d0*((vn1+abs(vn1))*kt(i,j,k-1)+(vn2-abs(vn2))*kt(i,j,k))  ! 一阶 L-F格式
     wfi=0.5d0*((vn1+abs(vn1))*wt(i,j,k-1)+(vn2-abs(vn2))*wt(i,j,k))  ! 一阶 L-F格式
    
	 Fluxk(i,j,k)= ( -kfi+muk*(ktx*s1x+kty*s1y+ktz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
 	 Fluxw(i,j,k)= ( -wfi+muw*(wtx*s1x+wty*s1y+wtz*s1z))*B%Sk(i,j,k)    ! 无粘+粘性通量
   
	enddo
    enddo
    enddo
!$OMP END DO
 
!$OMP  DO 
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       B%Res(6,i,j,k)=B%Res(6,i,j,k)+Fluxk(i,j,k+1)-Fluxk(i,j,k)           
       B%Res(7,i,j,k)=B%Res(7,i,j,k)+Fluxw(i,j,k+1)-Fluxw(i,j,k)           
	 enddo
     enddo
     enddo
!$OMP END DO 

 !$OMP DO   
   do k=0,nz
   do j=0,ny
   do i=0,nx
    B%mu(i,j,k)=B%mu(i,j,k)/Re
    B%mu_t(i,j,k)=B%mu_t(i,j,k)/Re
   enddo
   enddo
   enddo
!$OMP END DO   


!$OMP END PARALLEL 
!--------------源项---------------------------------------------------------
  deallocate(f1,Kt,Wt,Fluxk,Fluxw)

end  subroutine  Turbulence_model_SST

!------------------------------------------------------------------------------
   subroutine turbulence_model_BL(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,nz,ksub,i,j,k,kflag,i1,j1,k1,i2,j2,k2,i0,j0,k0,nMesh
   real(PRE_EC):: ui,vi,wi,uj,vj,wj,uk,vk,wk,ux,vx,wx,uy,vy,wy,uz,vz,wz
   real(PRE_EC):: ix,iy,iz,jx,jy,jz,kx,ky,kz,x0,y0,z0
   real(PRE_EC),allocatable,dimension(:,:,:):: omiga
   real(PRE_EC),allocatable,dimension(:):: Amu1d,Amut1d,d1d,u1d,yy,omiga1d
   integer,allocatable:: flag1(:,:,:) 
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   
   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny; nz=B%nz
   B%mu_t(:,:,:)=0.d0

! test if the block contains wall    
   kflag=0
   do ksub=1, B%subface
     if(B%bc_msg(ksub)%bc .eq. BC_WALL) kflag=1
   enddo  
   if(Kflag .eq. 0) return    ! No wall in this block

!  This Block Contains Wall
   allocate(omiga(0:nx,0:ny,0:nz),flag1(0:nx,0:ny,0:nz))
   omiga=0.d0
   flag1=0

!----- get Omiga (vorticity)  omiga=sqrt(omigax**2+omigay**2+omigaz**2) at the cell's center ------
! 采用中心差分求解
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (nx,ny,nz,uu,v,w,omiga,B)
   do k=1,nz-1
   do j=1,ny-1
   do i=1,nx-1

! 物理量对于计算坐标（下标）的导数
 
   ui=uu(i+1,j,k)-uu(i-1,j,k)            
   vi=v(i+1,j,k)-v(i-1,j,k)  
   wi=w(i+1,j,k)-w(i-1,j,k)  
   uj=uu(i,j+1,k)-uu(i,j-1,k)   
   vj=v(i,j+1,k)-v(i,j-1,k)
   wj=w(i,j+1,k)-w(i,j-1,k) 
   uk=uu(i,j,k+1)-uu(i,j,k-1)  
   vk=v(i,j,k+1)-v(i,j,k-1)
   wk=w(i,j,k+1)-w(i,j,k-1)  
 
   ix=B%ix0(i,j,k); iy=B%iy0(i,j,k); iz=B%iz0(i,j,k)
   jx=B%jx0(i,j,k); jy=B%jy0(i,j,k); jz=B%jz0(i,j,k)
   kx=B%kx0(i,j,k); ky=B%ky0(i,j,k); kz=B%kz0(i,j,k)

!----对物理坐标的偏导数----------------------------------------------
   ux=ui*ix+uj*jx+uk*kx
   vx=vi*ix+vj*jx+vk*kx
   wx=wi*ix+wj*jx+wk*kx
   uy=ui*iy+uj*jy+uk*ky
   vy=vi*iy+vj*jy+vk*ky
   wy=wi*iy+wj*jy+wk*ky
   uz=ui*iz+uj*jz+uk*kz
   vz=vi*iz+vj*jz+vk*kz
   wz=wi*iz+wj*jz+wk*kz
!-------------------------------------------------
    omiga(i,j,k)=sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
do ksub=1, B%subface
  Bc=> B%bc_msg(ksub)

 if(Bc%bc .eq. BC_WALL) then   ! Wall boundary
  ! 沿网格线一维处理（假设网格线垂直壁面）---------------------------------------
  if(Bc%face .eq. 1 .or. Bc%face .eq. 4) then   ! i+ or i- 
    allocate(yy(nx),Amu1d(nx),Amut1d(nx),d1d(nx),u1d(nx),omiga1d(nx))
    Amut1d=0.d0; Amu1d=0.d0  ! 初始化
   
    
    if(Bc%face .eq. 1) then
      i1=1; i2=0 
    else
      i1=nx-1; i2=nx 
    endif  
   
    do k=Bc%kb,Bc%ke-1
    do j=Bc%jb,Bc%je-1
       x0=(B%xc(i1,j,k)+B%xc(i2,j,k))*0.5d0 
       y0=(B%yc(i1,j,k)+B%yc(i2,j,k))*0.5d0
       z0=(B%zc(i1,j,k)+B%zc(i2,j,k))*0.5d0

     do i=1,nx-1
      if(Bc%face .eq. 1) then
       i0=i
      else
       i0=nx-i
      endif

       yy(i0)=sqrt((B%xc(i,j,k)-x0)**2+(B%yc(i,j,k)-y0)**2+(B%zc(i,j,k)-z0)**2)
       u1d(i0)=sqrt(uu(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
       d1d(i0)=d(i,j,k)  
       omiga1d(i0)=omiga(i,j,k)
       Amu1d(i0)=B%mu(i,j,k)
     enddo
 !  BL模型（一维）  
     call BL_model_1d(nx-1,yy,Amu1d,Amut1d,d1d,u1d,omiga1d)
   
     do i=1,nx-1
! 如果一个点处于多条线上，取粘性系数最小的值
      if(Bc%face .eq. 1) then
       i0=i
      else
       i0=nx-i
      endif
   
     if(flag1(i,j,k) .eq. 0) then
      flag1(i,j,k)=1
      B%mu_t(i,j,k)=Amut1d(i0)
     else
      B%mu_t(i,j,k)=min(B%mu_t(i,j,k),Amut1d(i0))
     endif
     enddo
   
    enddo
    enddo
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,omiga1d)

! !!! To set Amu_t=0 in the wall      设置虚网格点上的mut(-1)=-mut(1), 以保证壁面上mut=0 (mut=0.5*(mut(-1)+mut(1))

!  设置壁面第1层网格上的湍流粘性系数为0
    if(Bc%face .eq. 1) then
     B%mu_t(0,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(1,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0   
    else
     B%mu_t(nx,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(nx-1,Bc%jb:Bc%je-1,Bc%kb:Bc%ke-1)=0.d0
    endif


 else if(Bc%face .eq. 2 .or. Bc%face .eq. 4) then   ! face of j- or j+ 
    allocate(yy(ny),Amu1d(ny),Amut1d(ny),d1d(ny),u1d(ny),omiga1d(ny))
     Amut1d=0.d0; Amu1d=0.d0  ! 初始化

    
    if(Bc%face .eq. 2) then
      j1=1; j2=0 
    else
      j1=ny-1; j2=ny
    endif  
    do k=Bc%kb,Bc%ke-1
    do i=Bc%ib,Bc%ie-1
       x0=(B%xc(i,j1,k)+B%xc(i,j2,k))*0.5d0 
       y0=(B%yc(i,j1,k)+B%yc(i,j2,k))*0.5d0
       z0=(B%zc(i,j1,k)+B%zc(i,j2,k))*0.5d0

     do j=1,ny-1
      if(Bc%face .eq. 2) then
       j0=j
      else
       j0=ny-j
      endif

       yy(j0)=sqrt((B%xc(i,j,k)-x0)**2+(B%yc(i,j,k)-y0)**2+(B%zc(i,j,k)-z0)**2)
       u1d(j0)=sqrt(uu(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
       d1d(j0)=d(i,j,k)  
       omiga1d(j0)=omiga(i,j,k)
       Amu1d(j0)=B%mu(i,j,k)
     enddo
 !  BL模型（一维）  
     call BL_model_1d(ny-1,yy,Amu1d,Amut1d,d1d,u1d,omiga1d)
   
     do j=1,ny-1
! 如果一个点处于多条线上，取粘性系数最小的值
      if(Bc%face .eq. 2) then
       j0=j
      else
       j0=ny-j
      endif
   
     if(flag1(i,j,k) .eq. 0) then
      flag1(i,j,k)=1
      B%mu_t(i,j,k)=Amut1d(j0)
     else
      B%mu_t(i,j,k)=min(B%mu_t(i,j,k),Amut1d(j0))
     endif
     enddo
    enddo
    enddo
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,omiga1d)

    if(Bc%face .eq. 2) then
     B%mu_t(Bc%ib:Bc%ie-1, 0,  Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, 1,  Bc%kb:Bc%ke-1)=0.d0   
    else
     B%mu_t(Bc%ib:Bc%ie-1, ny,  Bc%kb:Bc%ke-1)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, ny-1, Bc%kb:Bc%ke-1)=0.d0
    endif




  else if(Bc%face .eq. 3 .or. Bc%face .eq. 6) then   ! face of k- or k+ 
    allocate(yy(nz),Amu1d(nz),Amut1d(nz),d1d(nz),u1d(nz),omiga1d(nz))
    Amut1d=0.d0; Amu1d=0.d0  ! 初始化

    
    if(Bc%face .eq. 3) then
      k1=1; k2=0 
    else
      k1=nz-1; j2=nz
    endif  
    do j=Bc%jb,Bc%je-1
    do i=Bc%ib,Bc%ie-1
       x0=(B%xc(i,j,k1)+B%xc(i,j,k2))*0.5d0 
       y0=(B%yc(i,j,k1)+B%yc(i,j,k2))*0.5d0
       z0=(B%zc(i,j,k1)+B%zc(i,j,k2))*0.5d0

     do k=1,nz-1
      if(Bc%face .eq. 3) then
        k0=k
      else
        k0=nz-k
      endif

       yy(k0)=sqrt((B%xc(i,j,k)-x0)**2+(B%yc(i,j,k)-y0)**2+(B%zc(i,j,k)-z0)**2)
       u1d(k0)=sqrt(uu(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) 
       d1d(k0)=d(i,j,k)  
       omiga1d(k0)=omiga(i,j,k)
       Amu1d(k0)=B%mu(i,j,k)
     enddo
 !  BL模型（一维）  
     call BL_model_1d(nz-1,yy,Amu1d,Amut1d,d1d,u1d,omiga1d)
   
     do k=1,nz-1
! 如果一个点处于多条线上，取粘性系数最小的值
      if(Bc%face .eq. 3) then
       k0=k
      else
       k0=nz-k
      endif
   
     if(flag1(i,j,k) .eq. 0) then
      flag1(i,j,k)=1
      B%mu_t(i,j,k)=Amut1d(k0)
     else
      B%mu_t(i,j,k)=min(B%mu_t(i,j,k),Amut1d(k0))
     endif
     enddo
   
    enddo
    enddo
  
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,omiga1d)

    if(Bc%face .eq. 3) then
     B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, 0)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, 1)=0.d0   
    else
     B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, nz)=0.d0
	 B%mu_t(Bc%ib:Bc%ie-1, Bc%jb:Bc%je-1, nz-1)=0.d0
    endif
   
  endif
 
 endif
 enddo


   call Amut_boundary(nMesh,mBlock)
   deallocate(omiga,flag1)

end


!c------------------------------------------------------------------------
! B-L model of turbulence
! Ref:  Wilox DC. Turbulence Modeling for CFD (2nd Edition), p77
   subroutine BL_model_1d(ny,yy,Amu,Amu_t,d,u,omiga)
      use precision_EC
	  implicit none
      integer ny,j,Iflag
      real(PRE_EC),dimension(ny) :: yy, Amu,Amu_t,d,u,omiga
      real(PRE_EC),parameter::  AP=26.,Ccp=1.6,Ckleb=0.3,Cwk=0.25d0,AKT=0.4,AK=0.0168  ! Cwk= 1.d0
      real(PRE_EC):: Tw,Ret,Fmax,etamax,Udif,etap,FF,Fwak,bl,Fkleb,Visti,Visto
      
           TW=abs(Amu(1)*omiga(1))
            do j=1,Ny
             if(abs(Amu(j)*omiga(j)) .gt. TW) Tw=abs(Amu(j)*omiga(j))  
            enddo
             Ret=sqrt(d(1)*TW)/Amu(1)

           Fmax=0.d0 ;   etamax=0.d0 ;     Udif=0.d0
  

 !----------------------------------
           do j=1,Ny
            if(u(j) .gt. Udif) Udif=u(j)
             etap=yy(j)*Ret
             FF=yy(j)*abs(omiga(j))*(1.d0-exp(-etap/AP))

            if(FF.gt.Fmax) then
             Fmax=FF
             etamax=yy(j)   
            endif

!            if(FF.gt.Fmax) then
!             Fmax=FF
!             etamax=yy(j)    
!            else
!             goto 100    ! Find the first peak of F(y)   ! 只要第1个峰值         
!            endif

           enddo
100        continue

           IFlag=0
           Fwak=min(etamax*Fmax,Cwk*etamax*Udif*Udif/Fmax)
           do j=1,Ny
            etap=Ret*yy(j)
            bl=AKT*yy(j)*(1.d0-exp(-etap/AP))
            visti=d(j)*bl*bl*abs(omiga(j))
            Fkleb=1.d0/(1.d0+5.5d0*(Ckleb*yy(j)/etamax)**6)
            visto=AK*Ccp*d(j)*Fwak*Fkleb
            if(abs(visto).lt.abs(visti)) IFlag=1
            if(Iflag.eq.0) then
             Amu_t(j)=visti
            else
             Amu_t(j)=visto
            endif
           enddo 

  end

!  2017-7-11  对流场进行限制 （防止物理量超界）
!-------------------------------------------
  subroutine limit_flow(nMesh)
   use Global_var
   implicit none
   integer::nMesh
   call   limit_flow_U(nMesh)             ! Limit d,u,v,w,p
   if(Mesh(nMesh)%NVAR == 6) then 
    call limit_flow_SA (nMesh)       ! limit U(6) (usuall for SA)
   endif

   if(Mesh(nMesh)%NVAR == 7) then 
!    call limit_flow_U7       ! limit U(6) (usuall for SA)
     ! you can add your code here  !!!!!
   endif
  
  
  end

!----------------------------------------
!  限定压力及密度的变化幅度 (see CFL3D User's manual: p236, Time Advancement)
!  变化幅度超过阈值(如, -20%) 则进行特殊处理

   subroutine limit_flow_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,nx,ny,nz,i,j,k
   real(PRE_EC):: dn,un,vn,wn,pn,dd,dp
   real(PRE_EC),parameter:: alfac=-0.2d0, phic=2.d0     ! (Phic > 1)
   real(PRE_EC),dimension(:,:,:),pointer:: d1,u1,v1,w1,p1
   integer:: ia,ja,ka,i1,j1,k1,kn,Nneg,Kneg(3,10)
   real(PRE_EC):: d0,u0,v0,w0,p0, sn1
   ! alfac 限定值； phic 目标值 (p1 -> pn/phic)
   
   Type (Block_TYPE),pointer:: B

   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz

     allocate(d1(nx-1,ny-1,nz-1),u1(nx-1,ny-1,nz-1),v1(nx-1,ny-1,nz-1), &
	          w1(nx-1,ny-1,nz-1),p1(nx-1,ny-1,nz-1))

!--------------------------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,gamma,d1,u1,v1,w1,p1)
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
 !        dn=B%Un(1,i,j,k)
 !        un=B%Un(2,i,j,k)/dn
 !        vn=B%Un(3,i,j,k)/dn
 !        wn=B%Un(4,i,j,k)/dn
 !        pn=(B%Un(5,i,j,k)-0.5d0*dn*(un*un+vn*vn+wn*wn)) *(gamma-1.d0)

         d1(i,j,k)=B%U(1,i,j,k)
		 u1(i,j,k)=B%U(2,i,j,k)/d1(i,j,k)
         v1(i,j,k)=B%U(3,i,j,k)/d1(i,j,k)
         w1(i,j,k)=B%U(4,i,j,k)/d1(i,j,k)
         p1(i,j,k)=(B%U(5,i,j,k)-0.5d0*d1(i,j,k)    &
		           *(u1(i,j,k)*u1(i,j,k)+v1(i,j,k)*v1(i,j,k)+w1(i,j,k)*w1(i,j,k))) *(gamma-1.d0)

!        dd=d1(i,j,k)-dn
!        dp=p1(i,j,k)-pn
     
!        if( dd/dn < alfac .or. dp/pn < alfac) then
!	      if(dd/dn < alfac)   d1(i,j,k)=dn+dd/(1.d0+phic*(alfac+abs(dd/dn)))
!		  if(dp/pn < alfac)   p1(i,j,k)=pn+dp/(1.d0+phic*(alfac+abs(dp/pn)))
! 		   B%U(1,i,j,k)=d1(i,j,k)
!		   B%U(2,i,j,k)=d1(i,j,k)*u1(i,j,k)
!		   B%U(3,i,j,k)=d1(i,j,k)*v1(i,j,k)
!		   B%U(4,i,j,k)=d1(i,j,k)*w1(i,j,k)
!		   B%U(5,i,j,k)=p1(i,j,k)/(gamma-1.d0)+   &
!		          0.5d0*d1(i,j,k)*(u1(i,j,k)**2+v1(i,j,k)**2+w1(i,j,k)**2)
!  	     endif
	 
	  enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO   
	  
	   Nneg=0  		   
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,gamma,d1,u1,v1,w1,p1,Nneg,Kneg,Ldmin,Lpmin,Ldmax,Lpmax,Lumax)
       do k=1,nz-1
       do j=1,ny-1
       do i=1,nx-1
        if(d1(i,j,k) <Ldmin .or. p1(i,j,k) < Lpmin .or. d1(i,j,k) > Ldmax .or. p1(i,j,k) > Lpmax    &
		  .or. abs(u1(i,j,k)) > Lumax .or. abs(v1(i,j,k)) > Lumax  .or. abs(w1(i,j,k)) > Lumax   ) then                 ! 物理量超限
	
	
		  Nneg=Nneg+1
          if(Nneg <=10) then        ! The location of Limited point
		   Kneg(1,Nneg)=i; Kneg(2,Nneg)=j; Kneg(3,Nneg)=k
		  endif
		  
		  d0=0.d0; u0=0.d0; v0=0.d0 ; w0=0.d0; p0=0.d0; kn=0		   		   
          do ka=-1,1
		  do ja=-1,1
		  do ia=-1,1
		    k1=k+ka ; j1=j+ja ; i1=i+ia
			if(i1 > 0 .and. i1< nx .and. j1 > 0 .and. j1< ny .and. k1>0 .and. k1<nz ) then
            if( .not. (d1(i1,j1,k1) <Ldmin .or. p1(i1,j1,k1) < Lpmin .or. d1(i1,j1,k1) > Ldmax .or. p1(i1,j1,k1) > Lpmax   &
		       .or. abs(u1(i1,j1,k1)) > Lumax .or. abs(v1(i1,j1,k1)) > Lumax  .or. abs(w1(i1,j1,k1)) > Lumax )   ) then

			d0=d0+d1(i1,j1,k1)
			u0=u0+u1(i1,j1,k1)
			v0=v0+v1(i1,j1,k1)
			w0=w0+w1(i1,j1,k1)
			p0=p0+p1(i1,j1,k1)
			kn=kn+1
		    endif
			endif
		   enddo
		   enddo
		   enddo

		   if( kn ==0 ) then         ! 周围全部为“坏点”
			 d0=d1(i,j,k)
			 u0=u1(i,j,k)
			 v0=v1(i,j,k)
			 w0=w1(i,j,k)
			 p0=p1(i,j,k)
			 
			 if( d0 < Ldmin)  d0=Ldmin
			 if( d0 > Ldmax)  d0=Ldmax
			 if( p0 < Lpmin)  p0=Lpmin
			 if( p0 > Lpmax)  p0=Lpmax
			 if(abs(u0) > Lumax ) u0=sign(Lumax,u0)
			 if(abs(v0) > Lumax ) v0=sign(Lumax,v0)
			 if(abs(w0) > Lumax ) w0=sign(Lumax,w0)
		 
		   else
             sn1=1.d0/(kn)
             d0=d0*sn1
			 u0=u0*sn1
			 v0=v0*sn1
			 w0=w0*sn1
			 p0=p0*sn1
            endif

		   B%U(1,i,j,k)=d0
		   B%U(2,i,j,k)=d0*u0
		   B%U(3,i,j,k)=d0*v0
		   B%U(4,i,j,k)=d0*w0
		   B%U(5,i,j,k)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0+w0*w0)
       endif
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO   
    
   if(Nneg > 0) then
     print*, "------------------------------------------------------"
	 print*, "Limters ..."
	 print*, Ldmin,Ldmax,Lpmin,Lpmax,Lumax
	 print*, "Warning !!! Mesh, Block=", nMesh, B%block_no ,  "has", Nneg, " Limitted points !!!!" 
	 do k=1,Min(Nneg,10)
	 print*, (Kneg(i,k),i=1,3)
     i1=Kneg(1,k); j1=Kneg(2,k); k1=Kneg(3,k)
	 print*, d1(i1,j1,k1),u1(i1,j1,k1),v1(i1,j1,k1),w1(i1,j1,k1),p1(i1,j1,k1)
	 enddo

     open(101,file="error.log",position="append" )
     write(101,*) "------------------------------------------------------"
	 write(101,*) "Warning !!! Mesh, Block=", nMesh, B%block_no ,  "has", Nneg, " Limitted points !!!!" 
	 do k=1,Min(Nneg,10)
	 write(101,*) (Kneg(i,k),i=1,3)
	 enddo
     close(101)

     B%IF_OverLimit=1                ! 设定物理量超限标志
   endif
    
   deallocate(d1,u1,v1,w1,p1)
  enddo
  end subroutine limit_flow_U





   subroutine limit_flow_SA(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,nx,ny,nz,i,j,k,kn,ia,ja,ka,i1,j1,k1
   real(PRE_EC)::  s0
   real(PRE_EC),parameter:: SAmin=1.d-8
    
   Type (Block_TYPE),pointer:: B

   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     nx=B%nx; ny=B%ny; nz=B%nz

!--------------------------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(nx,ny,nz,B,LSAmax)
	
     do k=1,nz-1
     do j=1,ny-1
     do i=1,nx-1
       if(B%U(6,i,j,k) < SAmin) B%U(6,i,j,k)=SAmin
       if(B%U(6,i,j,k) > LSAmax ) then
	  	 
   	    s0=0.d0; kn=0		   		   
         do ka=-1,1
	     do ja=-1,1
	     do ia=-1,1
	       k1=k+ka ; j1=j+ja ; i1=i+ia
	       if(i1 > 0 .and. i1< nx .and. j1 > 0 .and. j1< ny .and. k1>0 .and. k1<nz ) then
           if( B%U(6,i1,j1,k1) >=SAmin .and. B%U(6,i1,j1,k1) <= LSAmax  ) then
			s0=s0+B%U(6,i1,j1,k1)
			kn=kn+1
		    endif
			endif
		  enddo
		  enddo
		  enddo

		  if( kn ==0 ) then         ! 周围全部为“坏点”
            B%U(6,i,j,k)=LSAmax			 
		  else
		    B%U(6,i,j,k)=s0/kn
          endif
        endif
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO   
  enddo
  end subroutine limit_flow_SA

!  filtering data 
! Revised by Li Xinliang, 2013-10-4
   
! 高精度滤波（4阶精度）
  subroutine Filtering_oneMesh(nMesh)     
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1
   Type (Block_TYPE),pointer:: B

   if(my_id .eq. 0) print*, "filtering ......"
 
   NVAR1=Mesh(nMesh)%NVAR
   do mBlock=1,Mesh(nMesh)%Num_Block
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
    allocate(f(NVAR1,nx,ny,nz),f0(NVAR1,nx,ny,nz))
    
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
         f(1,i,j,k)= B%U(1,i,j,k)                 ! d
         f(2,i,j,k)= B%U(2,i,j,k)/B%U(1,i,j,k)    ! u
         f(3,i,j,k)= B%U(3,i,j,k)/B%U(1,i,j,k)    ! v
         f(4,i,j,k)= B%U(4,i,j,k)/B%U(1,i,j,k)    ! w
         f(5,i,j,k)=(B%U(5,i,j,k)-0.5d0*B%U(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2))*(gamma-1.d0)  !p
 !        f(5,i,j,k)=(B%U(5,i,j,k)-0.5d0*B%U(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2))/(Cv*f(1,i,j,k))  !T
 
        do m=6,NVAR1
         f(m,i,j,k)=B%U(m,i,j,k)
        enddo
	   
	   enddo
	   enddo
	   enddo
       


       call filter_x3d(nMesh,mBlock)
       call filter_y3d(nMesh,mBlock)
       call filter_z3d(nMesh,mBlock)

      do k=1,nz
      do j=1,ny
	  do i=1,nx
      B%U(1,i,j,k)=f(1,i,j,k)
      B%U(2,i,j,k)=f(1,i,j,k)*f(2,i,j,k)
      B%U(3,i,j,k)=f(1,i,j,k)*f(3,i,j,k)
      B%U(4,i,j,k)=f(1,i,j,k)*f(4,i,j,k)
      B%U(5,i,j,k)=f(5,i,j,k)/(gamma-1.d0)+0.5d0*f(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2)
!      B%U(5,i,j,k)=Cv*f(1,i,j,k)*f(5,i,j,k)+0.5d0*f(1,i,j,k)*(f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2)
       do m=6,NVAR1
	    B%U(m,i,j,k)=f(m,i,j,k)
	   enddo
	  enddo
      enddo
      enddo

      deallocate(f,f0)
    enddo

  end 

!---------------------------------------------------

  subroutine filter_x3d(nMesh,mBlock)      
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,i1,i2
   Type (Block_TYPE),pointer:: B
   integer,parameter:: KLP=4   ! 滤波的网格半宽度
   real(PRE_EC),parameter:: eps0=1.d-8,  a1= 1.d0/2.d0,  a2= 9.d0/32.d0,    a3=-1.d0/32.d0
   real(PRE_EC):: p1,p2,alpha
 
    NVAR1=Mesh(nMesh)%NVAR
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
       
	  i1=KLP
      i2=nx-KLP+1
       
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
       do m=1,NVAR1
	    f0(m,i,j,k)=f(m,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo

 !----Smooth indix 

        
		do k=1,nz
        do j=1,ny
        do i=i1,i2
            p1=dabs(f0(5,i+1,j,k)-f0(5,i  ,j,k))   ! f(5,:,:,:)=p()
            p2=dabs(f0(5,i  ,j,k)-f0(5,i-1,j,k))
		    
			alpha=0.1d0*min (dabs((p1-p2)/(p1+p2+eps0))**4 , 1.d0)
           do m=1,NVAR1
!		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
!			      alpha*(a1*f0(m,i,j,k)+a2*(f0(m,i+1,j,k)+f0(m,i-1,j,k))+a3*(f0(m,i+3,j,k)+f0(m,i-3,j,k)))
             f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
			      alpha*(29.d0/32.d0*f0(m,i,j,k)+1.d0/16.d0*(f0(m,i+1,j,k)+f0(m,i-1,j,k))-1.d0/64.d0*(f0(m,i+2,j,k)+f0(m,i-2,j,k)))

           enddo
		enddo
        enddo
        enddo
  
      end



!--------------------------------------------------------------------------

  subroutine filter_y3d(nMesh,mBlock)      
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,i1,i2
   Type (Block_TYPE),pointer:: B
   integer,parameter:: KLP=4   ! 滤波的网格半宽度
   real(PRE_EC),parameter:: eps0=1.d-8,  a1= 1.d0/2.d0,  a2= 9.d0/32.d0,    a3=-1.d0/32.d0
   real(PRE_EC):: p1,p2,alpha
 
    NVAR1=Mesh(nMesh)%NVAR
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
       
	  i1=KLP
      i2=ny-KLP+1
       
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
       do m=1,NVAR1
	    f0(m,i,j,k)=f(m,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo

 !----Smooth indix 

        
		do k=1,nz
        do j=i1,i2
        do i=1,nx
            p1=dabs(f0(5,i,j+1,k)-f0(5,i,  j,k))   ! f(5,:,:,:)=p()
            p2=dabs(f0(5,i,j,  k)-f0(5,i,  j-1,k))
		 
		    alpha=0.1d0*min (dabs((p1-p2)/(p1+p2+eps0))**4 , 1.d0)
           do m=1,NVAR1
!		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
!			      alpha*(a1*f0(m,i,j,k)+a2*(f0(m,i,j+1,k)+f0(m,i,j-1,k))+a3*(f0(m,i,j+3,k)+f0(m,i,j-3,k)))
             f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
			      alpha*(29.d0/32.d0*f0(m,i,j,k)+1.d0/16.d0*(f0(m,i,j+1,k)+f0(m,i,j-1,k))-1.d0/64.d0*(f0(m,i,j+2,k)+f0(m,i,j-2,k)))

		   enddo
		enddo
        enddo
        enddo
  
      end


!--------------------------------------------------------------------------

  subroutine filter_z3d(nMesh,mBlock)      
   use Global_Var
   use filting_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,i1,i2
   Type (Block_TYPE),pointer:: B
   integer,parameter:: KLP=4   ! 滤波的网格半宽度
   real(PRE_EC),parameter:: eps0=1.d-8,  a1= 1.d0/2.d0,  a2= 9.d0/32.d0,    a3=-1.d0/32.d0
   real(PRE_EC):: p1,p2,alpha
 
    NVAR1=Mesh(nMesh)%NVAR
    B=>Mesh(nMesh)%block(mBlock)
    nx=B%nx; ny=B%ny; nz=B%nz
       
	  i1=KLP
      i2=nz-KLP+1
       
       do k=1,nz
	   do j=1,ny
	   do i=1,nx
       do m=1,NVAR1
	    f0(m,i,j,k)=f(m,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo

 !----Smooth indix 

        
		do k=i1,i2    
        do j=1,ny
        do i=1,nx
            p1=dabs(f0(5,i,j,k+1)-f0(5,i,  j,k))   ! f(5,:,:,:)=p()
            p2=dabs(f0(5,i,j,  k)-f0(5,i,  j,k-1))
		    alpha=0.1d0* min (dabs((p1-p2)/(p1+p2+eps0))**4 , 1.d0)
           do m=1,NVAR1
!		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
!			      alpha*(a1*f0(m,i,j,k)+a2*(f0(m,i,j,k+1)+f0(m,i,j,k-1))+a3*(f0(m,i,j,k+3)+f0(m,i,j,k-3)))
		     f(m,i,j,k)=(1.d0-alpha)*f0(m,i,j,k)+ &  
			      alpha*(29.d0/32.d0*f0(m,i,j,k)+1.d0/16.d0*(f0(m,i,j,k+1)+f0(m,i,j,k-1))-1.d0/64.d0*(f0(m,i,j,k+2)+f0(m,i,j,k-2)))

           enddo
		enddo
        enddo
        enddo
  
      end

end module mod_struct_solver
