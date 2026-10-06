!===============================================================================
! mod_struct_time.f90 -- structured solver: time-stepping helpers
! Encapsulates sub_time_acceleraction.f90 + sub_LU_SGS.f90 routines.
!===============================================================================
  module mod_struct_time
  contains
 !  加速收敛技术
 !   1) 局部时间步长
 !   2) 残差光顺

 !--------------------------------------------------------------------------------------------------
 ! 计算谱半径 see: Blazek's book, p189-190    
     subroutine comput_Lijk(nMesh,mBlock) 
     use Global_Var
	 use Flow_Var
     implicit none
     integer:: mBlock,nx,ny,nz,i,j,k,nMesh
     real(PRE_EC) D0,un,S0,vol1
     Type (Block_TYPE),pointer:: B

     B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny; nz=B%nz

! $OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,nz,gamma,PrL,Prt,uu,v,w,cc,B,Lci,Lvi,Lcj,Lvj,Lck,Lvk)
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
       if(If_viscous .eq. 1) then
	    D0=max(gamma,4._PRE_EC/3._PRE_EC)/d(i,j,k)*(B%mu(i,j,k)/PrL+B%mu_t(i,j,k)/Prt)
       else
	    D0=0.d0
	   endif
	   vol1=1.d0/B%vol(i,j,k)
       
       un=0.5d0*( uu(i,j,k)*(B%ni1(i,j,k)+B%ni1(i+1,j,k))+ &
                   v(i,j,k)*(B%ni2(i,j,k)+B%ni2(i+1,j,k))+ &
                   w(i,j,k)*(B%ni3(i,j,k)+B%ni3(i+1,j,k)) )
       S0=(B%si(i,j,k)+B%si(i+1,j,k))*0.5
       Lci(i,j,k)=(abs(un)+cc(i,j,k))*S0
       Lvi(i,j,k)=D0*S0*S0*vol1

       un=0.5d0*( uu(i,j,k)*(B%nj1(i,j,k)+B%nj1(i,j+1,k))+ &
                   v(i,j,k)*(B%nj2(i,j,k)+B%nj2(i,j+1,k))+ &
                   w(i,j,k)*(B%nj3(i,j,k)+B%nj3(i,j+1,k)) )
       S0=(B%sj(i,j,k)+B%sj(i,j+1,k))*0.5
       Lcj(i,j,k)=(abs(un)+cc(i,j,k))*S0
       Lvj(i,j,k)=D0*S0*S0*vol1


       un=0.5d0*( uu(i,j,k)*(B%nk1(i,j,k)+B%nk1(i,j,k+1))+ &
                   v(i,j,k)*(B%nk2(i,j,k)+B%nk2(i,j,k+1))+ &
                   w(i,j,k)*(B%nk3(i,j,k)+B%nk3(i,j,k+1)) )
       S0=(B%sk(i,j,k)+B%sk(i,j,k+1))*0.5
       Lck(i,j,k)=(abs(un)+cc(i,j,k))*S0
       Lvk(i,j,k)=D0*S0*S0*vol1

       enddo
       enddo
       enddo
! $OMP END PARALLEL DO 

    end subroutine comput_Lijk


!---------------------------------------------------------------------------------
! 计算（当地）时间步长  ! J. Blazek, P.190
  subroutine comput_dt(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   integer  nMesh,mBlock,nx,ny,nz,i,j,k
   real(PRE_EC) D0,un,S0,vol1,dt_fac
   real(PRE_EC):: C

   Type (Block_TYPE),pointer:: B   
   C=1.d0
   B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
   nx=B%nx; ny=B%ny; nz=B%nz

   if( B%IF_OverLimit .eq. 0) then                ! 物理量超限
     dt_fac=1.0
   else
     dt_fac= 0.1d0              ! 时间步长降低10倍
	 print*, " ------- In Block No. ", B%Block_no,  "flow OverLimit, time step 1/10 ----"
   endif
 

   if(Iflag_local_dt .eq. 1)    then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
       do k=1,nz-1
       do j=1,ny-1
       do i=1,nx-1
         B%dt(i,j,k)=dt_fac*CFL*B%Vol(i,j,k)  &
		           /(Lci(i,j,k)+Lcj(i,j,k)+Lck(i,j,k)+C*(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k)))

         if(If_dtime_mesh .eq. 1) B%dt(i,j,k)=B%dt(i,j,k)*B%dtime_mesh(i,j,k)           ! 根据网格质量，修正时间步长

		 if(B%dt(i,j,k) .gt. dtmax) B%dt(i,j,k)=dtmax
         if(B%dt(i,j,k) .lt. dtmin) B%dt(i,j,k)=dtmin
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO 

   else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
      do k=1,nz-1
      do j=1,ny-1
      do i=1,nx-1
          B%dt(i,j,k)=dt_global ! 全局时间步长法
      enddo
	  enddo
      enddo
!$OMP END PARALLEL DO 
  endif


  end subroutine comput_dt

 
 
! 迎风隐式残差光顺
! Ver 0.76  Upwind Implicit Residual smoothing
! 见J. Blazek's Book "Computational Fluid Dynamic: Principle and Application" 
   Subroutine Residual_smoothing(nMesh,mBlock)
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: mBlock,nx,ny,nz,i,j,k,m,nMesh
   integer,parameter::Nmax=4000
   real(PRE_EC):: as(Nmax),bs(Nmax),cs(Nmax),R(Nmax),Rs(Nmax),ei,Mn
   real(PRE_EC),parameter:: epsl=1.d0
   Type (Block_TYPE),pointer:: B
   B => Mesh(nMesh)%Block(mBlock)                                         ! 指向其一块
   nx=B%nx;  ny= B%ny;  nz=B%nz

!--i-direction---------------------------
    do k=1,nz-1
    do j=1,ny-1
    do i=1,nx-1
      ei=epsl*min(1.0_PRE_EC,Lci(i,j,k)/Lcj(i,j,k),Lci(i,j,k)/Lck(i,j,k))
 !     法向Mach数     
	  Mn=0.5d0*( uu(i,j,k)*(B%ni1(i,j,k)+B%ni1(i+1,j,k))+ &
                 v(i,j,k)*(B%ni2(i,j,k)+B%ni2(i+1,j,k))+ &
                 w(i,j,k)*(B%ni3(i,j,k)+B%ni3(i+1,j,k)) )/cc(i,j,k)    
      if(Mn .gt. 1.0) then
       as(i)=-ei;  bs(i)=1.0+ei ;  cs(i)=0.0
      else if (Mn .lt. -1.0) then
	   as(i)=0.0; bs(i)=1.0+ei; cs(i)=-ei
	  else
	   as(i)=-ei; bs(i)=1.+2.*ei ; cs(i)=-ei
	  endif  
    enddo
     do m=1,5
       do i=1,nx-1
        R(i)=B%Res(m,i,j,k)
       enddo
       call tridiagonal(nx-1,as,bs,cs,R,Rs)    ! R残差； Rs光滑后的残差
       do i=1,nx-1
        B%Res(m,i,j,k)=Rs(i)
       enddo
     enddo
   enddo
   enddo

!---j- direction ------------------
   do k=1,nz-1
   do i=1,nx-1
   do j=1,ny-1
    ei=epsl*min(1.0_PRE_EC,Lcj(i,j,k)/Lci(i,j,k),Lcj(i,j,k)/Lck(i,j,k))

!     法向Mach数     
     Mn=0.5d0*( uu(i,j,k)*(B%nj1(i,j,k)+B%nj1(i,j+1,k))+ &
                v(i,j,k)*(B%nj2(i,j,k)+B%nj2(i,j+1,k))+ &
                w(i,j,k)*(B%nj3(i,j,k)+B%nj3(i,j+1,k)) )/cc(i,j,k)
     if(Mn .gt. 1.0) then
	  as(j)=-ei ; bs(j)=1.0+ei; cs(j)=0.0
     else if(Mn .lt. -1.0) then
	  as(j)=0.0; bs(j)=1.0+ei; cs(j)=-ei
	 else
	  as(j)=-ei; bs(j)=1.0+2.0*ei; cs(j)=-ei
	 endif
    enddo
      do m=1,5
       do j=1,ny-1
        R(j)=B%Res(m,i,j,k)
       enddo
       call tridiagonal(ny-1,as,bs,cs,R,Rs)
       do j=1,ny-1
        B%Res(m,i,j,k)=Rs(j)
       enddo
      enddo
    enddo
    enddo
!---k- direction ----------------
   do j=1,ny-1
   do i=1,nx-1
      do k=1,nz-1
       ei=epsl*min(1.0_PRE_EC,Lck(i,j,k)/Lci(i,j,k),Lck(i,j,k)/Lcj(i,j,k))
        Mn=0.5d0*( uu(i,j,k)*(B%nk1(i,j,k)+B%nk1(i,j,k+1))+ &
                   v(i,j,k)*(B%nk2(i,j,k)+B%nk2(i,j,k+1))+ &
                   w(i,j,k)*(B%nk3(i,j,k)+B%nk3(i,j,k+1)) )/cc(i,j,k)
      if(Mn .gt. 1.0) then
	   as(k)=-ei ; bs(k)=1.0+ei; cs(k)=0.0
      else if(Mn .lt. -1.0) then
	   as(k)=0.0; bs(k)=1.0+ei; cs(k)=-ei
	  else
	   as(k)=-ei; bs(k)=1.0+2.0*ei; cs(k)=-ei
	  endif
     enddo
      do m=1,5
       do k=1,nz-1
         R(k)=B%Res(m,i,j,k)
       enddo
       call tridiagonal(nz-1,as,bs,cs,R,Rs)
       do k=1,nz-1
        B%Res(m,i,j,k)=Rs(k)
       enddo
      enddo
    enddo
    enddo

   end subroutine Residual_smoothing
 

! 三对角方程组求解
! 求解 a(i)*Rs(i-1)+b(i)*Rs(i)+c(i)*Rs(i+1)=R(i)
   subroutine tridiagonal(n,a,b,c,R,Rs)
   use precision_EC
   implicit none
   integer:: n,i
   integer,parameter::Nmax=4000
   real(PRE_EC):: a(n),b(n),c(n),R(n),Rs(n),P(Nmax),Q(Nmax),tmp
   P(1)=0.d0; Q(1)=R(1)
   P(n)=0.d0; Q(n)=R(n)

   do i=2,n
    tmp=1.d0/(a(i)*P(i-1)+b(i))
    P(i)=-c(i)*tmp
    Q(i)=(R(i)-a(i)*Q(i-1))*tmp
   enddo
   Rs(n)=R(n)
   Rs(1)=R(1)
   do i=n-1,2,-1
   Rs(i)=P(i)*Rs(i+1)+Q(i)
   enddo
   end subroutine tridiagonal


!  采用LU-SGS方法，计算DU=U(n+1)-U(n)
!  Code by Li Xinliang, 2011-12-29
!-----------------------------------------------------------------------------------------
    subroutine  du_LU_SGS(nMesh,mBlock,Sfac1)                          ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
    use Global_Var
    use Flow_Var 
    implicit none
	integer:: nMesh,mBlock,NV,nx,ny,nz,plane,i,j,k,m
    real(PRE_EC),dimension(7)::alfa,dui,duj,duk,DF
    Type (Block_TYPE),pointer:: B
    real(PRE_EC):: Sfac1    ! 双时间步时使用

     NV=Mesh(nMesh)%NVAR
	 B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny; nz=B%nz

! LU-SGS的两次扫描
!----------------------------------
!   从i=1,j=1,k=1 到i=nx-1,j=ny-1,k=nz-1的扫描过程  (向上扫描过程)
!   扫描 i+j+k=plane 的平面
!   w_LU是松弛因子（1到2之间），增大w_LU会提高稳定性，但会降低收敛速度
   do plane=3,nx+ny+nz-3            

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(plane,nx,ny,nz,NV,B,Lci,Lcj,Lck,Lvi,Lvj,Lvk,gamma,w_LU,If_viscous)
     do k=1,nz-1
	 do j=1,ny-1
	 i=plane-k-j
	 if( i .lt. 1 .or. i .gt. nx-1) cycle    ! 超出了这个平面
       alfa(1:NV)=B%vol(i,j,k)/B%dt(i,j,k)+w_LU*(Lci(i,j,k)+Lcj(i,j,k)+Lck(i,j,k))     ! 对角线项
     if(If_viscous .eq. 1)  then
	   alfa(1:NV)=alfa(1:NV)+2.d0*(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k))          
       if(NV .eq. 6) then
	     alfa(6)=alfa(6)+(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k))      ! 再增大些对角线 （考虑到sigma_SA=2.d0/3.d0）
!	   else if(NV .eq. 7) then
!	      alfa(6)=alfa(6)+0.09*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)                          !处理源项刚性
!         alfa(7)=alfa(7)+2.d0*0.0828*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)
	   endif
     endif
	 	 alfa(1:NV)=alfa(1:NV)+Sfac1*B%vol(i,j,k)         ! 单时间步长Sfac1=0
		  		   
	 if(i.ne. 1) then
!                                              通量的差量，用来近似计算A*W (See Blazek's book, page 208)
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i-1,j,k),B%DU(1:NV,i-1,j,k),B%ni1(i,j,k),B%ni2(i,j,k),B%ni3(i,j,k),gamma)  
       dui(1:NV)=0.5d0*(DF(1:NV)*B%si(i,j,k)+w_LU*Lci(i-1,j,k)*B%DU(1:NV,i-1,j,k))
       if(If_viscous .eq. 1)    dui(1:NV)=dui(1:NV)+Lvi(i-1,j,k)*B%DU(1:NV,i-1,j,k)         
      else
	   dui(1:NV)=0.d0                             ! 左侧没有点
      endif
	 
	 if(j.ne.1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j-1,k),B%DU(1:NV,i,j-1,k),B%nj1(i,j,k),B%nj2(i,j,k),B%nj3(i,j,k),gamma)  
        duj(1:NV)=0.5d0*(DF(1:NV)*B%sj(i,j,k)+w_LU*Lcj(i,j-1,k)*B%DU(1:NV,i,j-1,k))
        if(If_viscous .eq. 1)    duj(1:NV)=duj(1:NV)+Lvj(i,j-1,k)*B%DU(1:NV,i,j-1,k)   ! 2012-2-29
	 else
	   duj(1:NV)=0.d0
	 endif

    if(k .ne. 1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j,k-1),B%DU(1:NV,i,j,k-1),B%nk1(i,j,k),B%nk2(i,j,k),B%nk3(i,j,k),gamma)  
        duk(1:NV)=0.5d0*(DF(1:NV)*B%sk(i,j,k)+w_LU*Lck(i,j,k-1)*B%DU(1:NV,i,j,k-1))
        if(If_viscous .eq. 1)    duk(1:NV)=duk(1:NV)+Lvk(i,j,k-1)*B%DU(1:NV,i,j,k-1)   ! 2012-2-29
	else
	   duk(1:NV)=0.d0
	endif
	do m=1,NV
    B%DU(m,i,j,k)=(B%Res(m,i,j,k)+dui(m)+duj(m)+duk(m))/alfa(m)
	enddo
    enddo
    enddo
!$OMP END PARALLEL DO 
   enddo
!----------------------------------------------------------
!  从 (nx-1,ny-1,nz-1)到(1,1,1)的扫描过程 （向下扫描过程）
!  plane=i+j+k
   do plane=nx+ny+nz-3,3,-1   

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(plane,nx,ny,nz,NV,B,Lci,Lcj,Lck,Lvi,Lvj,Lvk,gamma,w_LU,If_viscous)
     do k=nz-1,1,-1
	 do j=ny-1,1,-1
	 i=plane-k-j
	 if( i .lt. 1 .or. i .gt. nx-1) cycle            ! 超出了这个平面
      alfa(1:NV)=B%vol(i,j,k)/B%dt(i,j,k)+w_LU*(Lci(i,j,k)+Lcj(i,j,k)+Lck(i,j,k))
      if(If_viscous .eq. 1) then
	    alfa(1:NV)=alfa(1:NV)+2.d0*(Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k) )         
	   if(NV .eq. 6) then
	     alfa(6)=alfa(6)+Lvi(i,j,k)+Lvj(i,j,k)+Lvk(i,j,k)             ! Sigma_SA=2./3.
!	    elseif(NV .eq. 7) then
!	      alfa(6)=alfa(6)+0.09*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)
!         alfa(7)=alfa(7)+2.d0*0.0828*B%U(6,i,j,k)/B%U(1,i,j,k)*Re*B%vol(i,j,k)
	   endif
      endif
	 	 alfa(1:NV)=alfa(1:NV)+Sfac1*B%vol(i,j,k)         ! 单时间步长Sfac1=0

	 if(i.ne. nx-1) then
!                                              通量的差量，用来近似计算A*W (See Blazek's book, page 208)
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i+1,j,k),B%DU(1:NV,i+1,j,k),B%ni1(i,j,k),B%ni2(i,j,k),B%ni3(i,j,k),gamma)  
         dui(1:NV)=-0.5d0*(DF(1:NV)*B%si(i+1,j,k)-w_LU*Lci(i+1,j,k)*B%DU(1:NV,i+1,j,k))
         if(If_viscous .eq. 1)    dui(1:NV)=dui(1:NV)+Lvi(i+1,j,k)*B%DU(1:NV,i+1,j,k)
  	   
	   else
	   dui(1:NV)=0.d0
       endif
	 
	 if(j.ne. ny-1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j+1,k),B%DU(1:NV,i,j+1,k),B%nj1(i,j,k),B%nj2(i,j,k),B%nj3(i,j,k),gamma)  
       duj(1:NV)=-0.5d0*(DF(1:NV)*B%sj(i,j+1,k)-w_LU*Lcj(i,j+1,k)*B%DU(1:NV,i,j+1,k))
       if(If_viscous .eq. 1)    duj(1:NV)=duj(1:NV)+Lvj(i,j+1,k)*B%DU(1:NV,i,j+1,k)

	 else
	   duj(1:NV)=0.d0
	 endif

    if(k .ne. nz-1) then
       call comput_DFn(NV,DF(1:NV),B%U(1:NV,i,j,k+1),B%DU(1:NV,i,j,k+1),B%nk1(i,j,k),B%nk2(i,j,k),B%nk3(i,j,k),gamma)  
       duk(1:NV)=-0.5d0*(DF(1:NV)*B%sk(i,j,k+1)-w_LU*Lck(i,j,k+1)*B%DU(1:NV,i,j,k+1))
       if(If_viscous .eq. 1)    duk(1:NV)=duk(1:NV)+Lvk(i,j,k+1)*B%DU(1:NV,i,j,k+1)
	else
	   duk(1:NV)=0.d0
	endif
	do m=1,NV
      B%DU(m,i,j,k)=B%DU(m,i,j,k)+(dui(m)+duj(m)+duk(m))/alfa(m)
	enddo
   enddo
   enddo
!$OMP END PARALLEL DO 

   enddo

   end subroutine  du_LU_SGS

!-----------------------------------------------------------------------


!  计算通量的差量 DF=F(Unew)-F(Uold),  LU-SGS方法中使用，用来近似A*DU  
    subroutine comput_DFn(NVAR1,DF,U,DU,n1,n2,n3,gamma)
    use precision_EC
	implicit none
    integer:: NVAR1
    real(PRE_EC),dimension(NVAR1):: DF,U,DU,U2
	real(PRE_EC):: n1,n2,n3,un1,un2,gamma,p1,p2
    U2=U+DU
	un1=(U(2)*n1+U(3)*n2+U(4)*n3)/U(1)              !un
	p1=(gamma-1.d0)*(U(5)-0.5d0*(U(2)**2+U(3)**2+U(4)**2)/U(1))
	un2=(U2(2)*n1+U2(3)*n2+U2(4)*n3)/U2(1)          !un
	p2=(gamma-1.d0)*(U2(5)-0.5d0*(U2(2)**2+U2(3)**2+U2(4)**2)/U2(1))

    DF(1)=U2(1)*un2-U(1)*un1                  ! d*un
	DF(2)=(U2(2)*un2+p2*n1)-(U(2)*un1+p1*n1)
	DF(3)=(U2(3)*un2+p2*n2)-(U(3)*un1+p1*n2)
	DF(4)=(U2(4)*un2+p2*n3)-(U(4)*un1+p1*n3)
	DF(5)=(U2(5)+p2)*un2-(U(5)+p1)*un1
	if(NVAR1 .eq. 6) then
	  DF(6)=U2(6)*un2-U(6)*un1
	else if(NVAR1 .eq. 7) then
      DF(6)=U2(6)*un2-U(6)*un1   ! k方程的对流通量
	  DF(7)=U2(7)*un2-U(7)*un1   ! w方程的对流通量
    endif
	end subroutine comput_dFn
  end module mod_struct_time
