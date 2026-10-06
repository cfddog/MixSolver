!===============================================================================
! mod_struct_bc.f90 -- structured solver: boundary conditions
! Encapsulates sub_boundary.f90 + sub_boundary_user.f90 routines.
!===============================================================================
  module mod_struct_bc
  contains
! modules for Boundary layer condition
! Since Ver 1.1, 边界（非内部链接边界）只采用1层虚网格
! 周期性条件通过 内边界实现 （某些情况下需要特殊处理）
!---------------------------------------------------------------------
! 处理边界条件（非内边界） （处理一套网格）
    subroutine Boundary_condition_onemesh(nMesh)
     use Global_Var
     use mod_interface, only: BC_INTERFACE
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: nMesh,mBlock,ksub
     integer,parameter:: FLAG_OUTLET=1 , FLAG_INLET=2

 ! ------------------------------------------------------ 
     do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
       do  ksub=1,B%subface
        Bc=> B%bc_msg(ksub)
        if(Bc%bc <= 0 )  cycle              ! 非内边界


          if( Bc%bc .eq. BC_Wall  .and. If_viscous .eq. 1 ) then
             call boundary_wall(nMesh,mBlock,ksub)
          else if( Bc%bc .eq. BC_Farfield  ) then
             call boundary_Farfield(nMesh,mBlock,ksub,0)
          else if ( Bc%bc .eq. BC_Inflow  ) then
             call boundary_Farfield(nMesh,mBlock,ksub, FLAG_INLET)
          else if( Bc%bc .eq. BC_Outflow) then
             call boundary_Farfield(nMesh,mBlock,ksub,FLAG_OUTLET)
          else if( Bc%bc .eq. BC_Symmetry .or. (Bc%bc .eq. BC_Wall .and. If_viscous .eq. 0) ) then
             call boundary_Symmetry_or_SlideWall(nMesh,mBlock,ksub)
          else if ( Bc%bc .eq. BC_Extrapolate ) then
             call boundary_Extrapolate(nMesh,mBlock,ksub)
          else if ( Bc%bc >=900 ) then
             call boundary_USER(nMesh,mBlock,ksub)
          else if ( Bc%bc .eq. BC_INTERFACE ) then
             cycle
          else
             print*, "The boundary condition is not supported!!!"
             print*, "Block_no is ", B%block_no, "bc=",Bc%bc
             stop
          endif
      enddo
      enddo
      
	 end subroutine Boundary_condition_onemesh

!------------------------------------------------------------
!-------------------------------------------------------------------  
! Wall boundary 
! 设定两层虚网格(Ghost Cell) 

    subroutine boundary_wall(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     real(PRE_EC):: d1,u1,v1,w1,p1,T1,d2,u2,v2,w2,p2,T2
     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,nMesh,NVAR1
     real(PRE_EC),parameter:: beta1_SST=0.075d0
     
     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      


    if(Bc%face .eq. 1 ) then   ! i- 
       do k=kb,ke-1
       do j=jb,je-1
		                          !   inner-cell             ! Ghost Cell
	    call wall_bound(NVAR1,B%U(:,1,j,k),B%U(:,0,j,k),Ma,gamma,Twall,B%mu(1,j,k),B%dw(1,j,k),Re)
	   enddo
       enddo
     
	 else if (Bc%face .eq. 4 ) then   ! i+
       do k=kb,ke-1
       do j=jb,je-1
	   	call wall_bound(NVAR1,B%U(:,ie-1,j,k),B%U(:,ie,j,k),Ma,gamma,Twall,B%mu(ie-1,j,k),B%dw(ie-1,j,k),Re)
	   enddo
       enddo
     
	 else if(Bc%face .eq. 2 ) then   !j-
       do k=kb,ke-1
       do i=ib,ie-1
	    call wall_bound(NVAR1,B%U(:,i,1,k),B%U(:,i,0,k),Ma,gamma,Twall,B%mu(i,1,k),B%dw(i,1,k),Re)
	   enddo
       enddo
     
	 else if(Bc%face .eq. 5 ) then   !j+
       do k=kb,ke-1
       do i=ib,ie-1
  	    call wall_bound(NVAR1,B%U(:,i,je-1,k),B%U(:,i,je,k),Ma,gamma,Twall,B%mu(i,je-1,k),B%dw(i,je-1,k),Re)
       enddo
       enddo
     
	 else if(Bc%face .eq. 3 ) then   !k-
       do j=jb,je-1
       do i=ib,ie-1
  	   	call wall_bound(NVAR1,B%U(:,i,j,1),B%U(:,i,j,0),Ma,gamma,Twall,B%mu(i,j,1),B%dw(i,j,1),Re)
       enddo
       enddo
    
	 else if(Bc%face .eq. 6 ) then   !k+
       do j=jb,je-1
       do i=ib,ie-1
  	    call wall_bound(NVAR1,B%U(:,i,j,ke-1),B%U(:,i,j,ke),Ma,gamma,Twall,B%mu(i,j,ke-1),B%dw(i,j,ke-1),Re)
       enddo
       enddo
 	 endif

         
    end subroutine boundary_wall
!-----------------------------------------------------
!-----------------------------------------------------
! 远场边界条件 （区分亚、超声速及出口、入口） 
! Ref. J. Blazek et al. "CFD principles and applications", P281-283
    subroutine boundary_Farfield(nMesh,mBlock,ksub,Flag)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc

     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,nMesh,Flag,NVAR1
     integer:: i1,j1,k1,i2,j2,k2,m
	 real(PRE_EC):: d_inf,u_inf,v_inf,w_inf,p_inf,n1,n2,n3 
     real(PRE_EC):: d1,u1,v1,w1,p1,c1,d2,u2,v2,w2,p2,pb,db,ub,vb,wb,Ma_n
     integer,parameter:: FLAG_OUTLET=1 , FLAG_INLET=2

!  本软件目前用来计算内流，给定无穷远条件
!  A_alfa,A_beta  攻角及侧滑角 （根据坐标方向确定）；  A_alfa  (x-y)平面内的倾角； A_beta (x-z)平面内的倾角 
     d_inf=1.d0
     u_inf=cos(A_alfa)*cos(A_beta)
     v_inf=sin(A_alfa)*cos(A_beta) 
     w_inf=sin(A_beta)
     p_inf=1.d0/(gamma*Ma*Ma)
     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)

     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je  ; kb=Bc%kb; ke=Bc%ke    

       do k=kb,ke
         do j=jb,je
           do i=ib,ie
!-----------------------------------------------------------------

 !  (i1,j1,k1) 是靠近边界的内点， (i2,j2,k2)  是边界外的1层 Ghost Cell点    
 ! (n1,n2,n3)为外法线方向

             if(Bc%face .eq. 1) then                 ! i- 面
               i1=i; j1=j; k1=k; i2=i-1 ; j2=j ; k2=k 
               n1=-B%ni1(i,j,k) ; n2=-B%ni2(i,j,k); n3=-B%ni3(i,j,k)   ! 外法线  
		     else if(Bc%face .eq. 2) then
               i1=i; j1=j; k1=k; i2=i;  j2=j-1 ; k2=k 
               n1=-B%nj1(i,j,k) ; n2=-B%nj2(i,j,k); n3=-B%nj3(i,j,k)     
             else if(Bc%face .eq. 3) then             
               i1=i; j1=j; k1=k; i2=i;  j2=j ; k2=k-1 
               n1=-B%nk1(i,j,k) ; n2=-B%nk2(i,j,k); n3=-B%nk3(i,j,k)     
             else if(Bc%face .eq. 4) then             ! i+ 面 (i=ibegin=iend=nx), i1=i-1 是内点, i2=i=nx是Ghost Cell
               i1=i-1; j1=j; k1=k;  i2=i; j2=j ; k2=k 
               n1=B%ni1(i,j,k) ; n2=B%ni2(i,j,k); n3=B%ni3(i,j,k)   ! 外法线  
             else if(Bc%face .eq. 5) then
               i1=i; j1=j-1; k1=k;  i2=i; j2=j ; k2=k 
               n1=B%nj1(i,j,k) ; n2=B%nj2(i,j,k); n3=B%nj3(i,j,k)     
		     else if(Bc%face .eq. 6) then
               i1=i; j1=j; k1=k-1;  i2=i; j2=j ; k2=k 
               n1=B%nk1(i,j,k) ; n2=B%nk2(i,j,k); n3=B%nk3(i,j,k)     
		     endif

            d1=B%U(1,i1,j1,k1); u1=B%U(2,i1,j1,k1)/d1; v1=B%U(3,i1,j1,k1)/d1; w1=B%U(4,i1,j1,k1)/d1
            p1=(B%U(5,i1,j1,k1)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)              ! 内点处的值
            c1=sqrt(gamma*p1/d1) 
!                   
          if(Flag .eq. FLAG_OUTLET) then   ! 强制为(超声速)出口
              d2=d1 ; u2=u1 ; v2=v1 ; w2=w1; p2=p1
		  else if(Flag .eq.  FLAG_INLET) then   ! 强制为(超声速)入口
	          d2=d_inf; u2=u_inf; v2=v_inf; w2=w_inf; p2=p_inf  
	  
!------------------------------------------------------------------------------
		  else                ! 普通远场边界条件


!   修改2012-5-21： 以来流 (而不是当地) Mach数判断，计算外流效果好； 计算内流尚待研究
            if( P_OUTLET <= -1.d0 ) then            ! 强制按照来流定义 
			 Ma_n=(u_inf*n1+v_inf*n2+w_inf*n3)*Ma   ! 法向Mach数， 强制以来流方向定义
            else
	         Ma_n=(u1*n1+v1*n2+w1*n3)/c1    ! 法向Mach数， 以内点值定义 （在边界层出口处效果不好）
		    endif
			      
		   if(Ma_n > 1.d0) then   ! 超声速出口  
              d2=d1 ; u2=u1 ; v2=v1 ; w2=w1; p2=p1
 !           else if (Ma_n .gt. 0.d0) then  !  亚声速出口
           else if (Ma_n > -1.d-6) then  !  亚声速出口  ( -1.d-6 为小量； 考虑到平行于来流的远场边界，按照出口处理更好）
		     if(P_OUTLET > 0.d0) then 
              pb=P_OUTLET 
              db=d1+(pb-p1)/(c1*c1)
              ub=u1+(p1-pb)/(d1*c1)*n1
              vb=v1+(p1-pb)/(d1*c1)*n2
              wb=w1+(p1-pb)/(d1*c1)*n3
              p2=2.d0*pb-p1; d2=2.d0*db-d1; u2=2.d0*ub-u1; v2=2.d0*vb-v1; w2=2.d0*wb-w1  !ub界面值，u2为Ghost Cell值  ub=(u1+u2)/2 
             else
               d2=d1 ; u2=u1 ; v2=v1 ; w2=w1; p2=p1   ! 外推
			 endif 
			else if (Ma_n > -1.d0) then  ! 亚声速入口 
              pb=0.5d0*(p1+p_inf-d1*c1*((u_inf-u1)*n1+(v_inf-v1)*n2+(w_inf-w1)*n3 ))
              db=d_inf+(pb-p_inf)/(c1*c1)
              ub=u_inf-(p_inf-pb)/(d1*c1)*n1
              vb=v_inf-(p_inf-pb)/(d1*c1)*n2
              wb=w_inf-(p_inf-pb)/(d1*c1)*n3
              p2=2.d0*pb-p1 ; d2=2.d0*db-d1 ; u2=2.d0*ub-u1 ; v2=2.d0*vb-v1 ; w2=2.d0*wb-w1
            else   ! 超声速入口
              d2=d_inf; u2=u_inf; v2=v_inf; w2=w_inf; p2=p_inf  
            endif  
          endif
			 
            B%U(1,i2,j2,k2)=d2
            B%U(2,i2,j2,k2)=d2*u2
            B%U(3,i2,j2,k2)=d2*v2
            B%U(4,i2,j2,k2)=d2*w2
            B%U(5,i2,j2,k2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)
 
    
!    标量（SA, SST）的边界条件
            if(Ma_n .gt. 0.d0 .or. Flag .eq. FLAG_OUTLET) then
			  if(NVAR1 .eq. 6)  then
			    B%U(6,i2,j2,k2)=B%U(6,i1,j1,k1)           ! vt
              else if(NVAR1 .eq. 7) then
			    B%U(6,i2,j2,k2)=B%U(6,i1,j1,k1)           
			    B%U(7,i2,j2,k2)=B%U(7,i1,j1,k1)           
			  endif   
            else
			  if(NVAR1 .eq. 6)  then
! see:             http://turbmodels.larc.nasa.gov/spalart.html
! 		 	    B%U(6,i2,j2,k2)=0.1d0/Re           ! vt
!		 	    B%U(6,i2,j2,k2)=5.d0/Re            ! vt 设定为层流粘性系数的5倍
		 	    B%U(6,i2,j2,k2)=5.d0               ! vt 设定为层流粘性系数的5倍 （0.98c以后版本)
              
			  else if(NVAR1 .eq. 7) then
			    B%U(6,i2,j2,k2)=d_inf*Kt_inf          ! 湍动能来流值 
			    B%U(7,i2,j2,k2)=d_inf*Wt_inf          ! 湍能比耗散率 
			  endif   
			endif

	      enddo
        enddo 
      enddo

    end subroutine boundary_Farfield
!----------------------------------------------------------------------------------------

!-------------------------------------------------------------------  
! 对称或滑移壁面条件
! Symmetry boundary condition or slide wall boundary condition
! 仅适用1层虚网格

    subroutine boundary_Symmetry_or_SlideWall(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1
     integer:: i1,j1,k1,i2,j2,k2
     real(PRE_EC):: n1,n2,n3,Vn

     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      

!   i2,i3 : 2层Ghost Cell;  i1,i4 2层内点
     if(Bc%face .eq. 1 .or. Bc%face .eq. 4) then   ! i- or i+
       if(Bc%face .eq. 1) then
         i=ib; i1=ib; i2=ib-1
       else
         i=ie; i1=ie-1; i2=ie 
       endif

       do k=kb,ke-1
         do j=jb,je-1
          do m=1,NVAR1
		    B%U(m,i2,j,k)=B%U(m,i1,j,k)
          enddo
           
		   n1=B%ni1(i,j,k) ; n2=B%ni2(i,j,k); n3=B%ni3(i,j,k)   ! 归一化法方向  
		   Vn=B%U(2,i1,j,k)*n1+B%U(3,i1,j,k)*n2+B%U(4,i1,j,k)*n3   ! 法向动量
!  对称边界条件，标量保持不变；垂直壁面的速度分量变号；平行壁面的速度分量不变
           B%U(2,i2,j,k)= B%U(2,i2,j,k)-2.d0*Vn*n1 
           B%U(3,i2,j,k)= B%U(3,i2,j,k)-2.d0*Vn*n2       
           B%U(4,i2,j,k)= B%U(4,i2,j,k)-2.d0*Vn*n3       
		 enddo
       enddo

     else if(Bc%face .eq. 2  .or. Bc%face .eq. 5) then   ! j- or j+
       if(Bc%face .eq. 2) then
         j=jb; j1=jb; j2=jb-1 
       else
         j=je; j1=je-1; j2=je 
       endif

       do k=kb,ke-1
         do i=ib,ie-1
           do m=1,NVAR1
		    B%U(m,i,j2,k)=B%U(m,i,j1,k)
           enddo

           n1=B%nj1(i,j,k) ; n2=B%nj2(i,j,k); n3=B%nj3(i,j,k)     
           Vn=B%U(2,i,j1,k)*n1+B%U(3,i,j1,k)*n2+B%U(4,i,j1,k)*n3
           B%U(2,i,j2,k)=  B%U(2,i,j2,k)-2.d0*Vn*n1
           B%U(3,i,j2,k)=  B%U(3,i,j2,k)-2.d0*Vn*n2
           B%U(4,i,j2,k)=  B%U(4,i,j2,k)-2.d0*Vn*n3

         enddo
       enddo

     else if(Bc%face .eq. 3 .or. Bc%face .eq. 6) then   ! k- or k+
       if(Bc%face .eq. 3) then
         k=kb; k1=kb; k2=kb-1 
       else
         k=ke; k1=ke-1; k2=ke
       endif

       do j=jb,je-1
         do i=ib,ie-1
           do m=1,NVAR1
		    B%U(m,i,j,k2)=B%U(m,i,j,k1)
           enddo
          
		   n1=B%nk1(i,j,k) ; n2=B%nk2(i,j,k); n3=B%nk3(i,j,k)     
		   Vn=B%U(2,i,j,k1)*n1+B%U(3,i,j,k1)*n2+B%U(4,i,j,k1)*n3
           B%U(2,i,j,k2)=  B%U(2,i,j,k2)-2.d0*Vn*n1      
           B%U(3,i,j,k2)=  B%U(3,i,j,k2)-2.d0*Vn*n2       
           B%U(4,i,j,k2)=  B%U(4,i,j,k2)-2.d0*Vn*n3       
		 enddo
       enddo

     endif
        
    end subroutine boundary_Symmetry_or_SlideWall

!-------------------------------------------------------------
   subroutine comput_origin_var(U1,U2,U3,U4,U5,d,u,v,w,p,T,Ma,gamma)    ! 根据守恒变量，计算基本变量
   use   precision_EC
   implicit none
   real(PRE_EC):: U1,U2,U3,U4,U5,d,u,v,w,T,p,Ma,gamma
    d=U1
    u=U2/d
    v=U3/d
    w=U4/d
    p=(U5-0.5d0*d*(u*u+v*v+w*w))*(gamma-1.d0)             
    T=gamma*Ma*Ma*p/d 
   end
!--------------------------------------------------------------
   subroutine comput_conser_var(U1,U2,U3,U4,U5,d,u,v,w,p,Ma,gamma)    ! 根据基本变量,计算守恒变量
   use   precision_EC
   implicit none
   real(PRE_EC):: U1,U2,U3,U4,U5,d,u,v,w,p,Ma,gamma
    U1=d
    U2=d*u
    U3=d*v
    U4=d*w
    U5=p/(gamma-1.d0)+0.5d0*d*(u*u+v*v+w*w)
   end




!  绝热或等温边界条件
!  设定两层Ghost Cell
!  U1: 内点 (i=1); Ug1: Ghost点 (i=0)
	 subroutine wall_bound(NVAR,U1,Ug1,Ma,gamma,Twall,mu1,dw,Re)
	 use   precision_EC
	 implicit none
	 integer:: NVAR
	 real(PRE_EC):: U1(NVAR),Ug1(NVAR)
	 real(PRE_EC):: d1,uu1,v1,w1,T1,p1,d2,uu2,v2,w2,T2,p2,Ma,gamma,Twall,mu1,dw,wt,Re
     real(PRE_EC),parameter:: beta1_SST=0.075d0

!    U1 内点(i=1)； Ug1 Ghost边界点(i=0)
     if(Twall .gt. 0.d0) then
       d1=U1(1)
       uu1=U1(2)/d1
       v1=U1(3)/d1
       w1=U1(4)/d1
       p1=(U1(5)-0.5d0*d1*(uu1*uu1+v1*v1+w1*w1))*(gamma-1.d0)             
       T1=gamma*Ma*Ma*p1/d1 
       p2=p1               ! 边界层假设，壁面处法向压力梯度为0
       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
       if( T2 .lt. 0.5*T1) T2=0.5*T1         !!! 2012-2-1  防止虚网格上的温度过低
  	   uu2=-uu1              ! 无滑移壁    (u1+u2)*0.5=0
	   v2=-v1
	   w2=-w1
	   d2=gamma*Ma*Ma*p2/T2 

       Ug1(1)=d2
       Ug1(2)=d2*uu2
       Ug1(3)=d2*v2
       Ug1(4)=d2*w2
       Ug1(5)=p2/(gamma-1.d0)+0.5d0*d2*(uu2*uu2+v2*v2+w2*w2)
  	  else
	    Ug1(1)=U1(1)
		Ug1(2)=-U1(2)
		Ug1(3)=-U1(3)
		Ug1(4)=-U1(4)
		Ug1(5)=U1(5)
      endif


		if(NVAR .eq. 6) then
		 Ug1(6)=0.d0                 ! vt for SA model
		else if(NVAR .eq. 7) then
!		 Ug1(6)=0.d0                              ! K for SST model
!		 Ug1(7)=60.d0*mu1/(beta1_SST*dw**2)       ! W for SST model
         wt=60.d0*mu1/(U1(1)*beta1_SST*dw**2*Re)
		 Ug1(6)=-U1(6)
         Ug1(7)=(U1(1)+Ug1(1))*wt-U1(7)

		endif
  
	 end

!-----------------------------------------------------
! 外插边界条件 
    subroutine boundary_Extrapolate(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc

     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,nMesh,NVAR1
     integer:: i1,j1,k1,i2,j2,k2,m

!  本软件目前用来计算内流，给定无穷远条件

     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)

     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je  ; kb=Bc%kb; ke=Bc%ke    

       do k=kb,ke
         do j=jb,je
           do i=ib,ie
!-----------------------------------------------------------------

 !  (i1,j1,k1) 是靠近边界的内点， (i2,j2,k2)  是边界外的1层 Ghost Cell点    

             if(Bc%face .eq. 1) then                 ! i- 面
               i1=i; j1=j; k1=k; i2=i-1 ; j2=j ; k2=k 
             else if(Bc%face .eq. 2) then
               i1=i; j1=j; k1=k; i2=i;  j2=j-1 ; k2=k 
             else if(Bc%face .eq. 3) then             
               i1=i; j1=j; k1=k; i2=i;  j2=j ; k2=k-1 
             else if(Bc%face .eq. 4) then             ! i+ 面 (i=ibegin=iend=nx), i1=i-1 是内点, i2=i=nx是Ghost Cell
               i1=i-1; j1=j; k1=k;  i2=i; j2=j ; k2=k 
             else if(Bc%face .eq. 5) then
               i1=i; j1=j-1; k1=k;  i2=i; j2=j ; k2=k 
             else if(Bc%face .eq. 6) then
               i1=i; j1=j; k1=k-1;  i2=i; j2=j ; k2=k 
             endif

            do m=1,NVAR1
              B%U(m,i2,j2,k2)=B%U(m,i1,j1,k1)         ! 简单外推 (1阶)
			enddo

	      enddo
        enddo 
      enddo

    end subroutine boundary_Extrapolate

!-----------------------------------------------------
! User-defined BC dispatcher (901: fixed inlet, 902: inlet time series,
! 903: blow/suction wall)
    subroutine boundary_user(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: mBlock,ksub,nMesh
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     Select case  (Bc%bc)
     case ( BC_USER_FixedInlet )
       call boundary_user_Inlet(nMesh,mBlock,ksub)
     case ( BC_USER_Inlet_time )
       call boundary_user_Inlet_time(nMesh,mBlock,ksub)
     case ( BC_USER_Blow_Suction_Wall )
       call boundary_user_blow_suction_Wall(nMesh,mBlock,ksub)
     case default
       print*, "This USER defined boundary is not supported !!!"
       stop
     end Select
    end

!------------------user-defined BC: fixed inlet profile------------------
! Reads inlet distribution from file inlet.xxx (xxx = block number)
! Uses 1 ghost-cell layer only

    subroutine boundary_user_Inlet(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1
     integer:: n1,n2,n3
     integer,save:: read_flag=0
     real(PRE_EC),allocatable,save,dimension(:,:,:,:):: uc,fc
     character(len=30):: fname

     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke
     n1=max(ie-ib,1)
     n2=max(je-jb,1)
     n3=max(ke-kb,1)

 !-----读取数据--------------------------------------   
	 if(read_flag==0) then
	  read_flag=1
      allocate(uc(n1,n2,n3,5),fc(5,n1,n2,n3))
      write(fname,"('inlet.',I3.3)") B%Block_no
	  open(99,file=fname,form="unformatted")
	  print*, "read file : ", fname
	  print*, "n1,n2,n3=",n1,n2,n3
	  do m=1,5
	   read(99,*) (((uc(i,j,k,m),i=1,n1),j=1,n2),k=1,n3)
      enddo
	  close(99)
	 endif

        do k=1,n3
		do j=1,n2
		do i=1,n1
          fc(1,i,j,k)=uc(i,j,k,1)
          fc(2,i,j,k)=uc(i,j,k,1)*uc(i,j,k,2)
          fc(3,i,j,k)=uc(i,j,k,1)*uc(i,j,k,3)
          fc(4,i,j,k)=uc(i,j,k,1)*uc(i,j,k,4)
          fc(5,i,j,k)=uc(i,j,k,1)*(Cv*uc(i,j,k,5)+0.5d0*(uc(i,j,k,2)**2+uc(i,j,k,3)**2+uc(i,j,k,4)**2))
        enddo
	    enddo
	    enddo

!-----------------------------------------------------------------

 ! 仅使用1层Ghost Cell    
             if(Bc%face .eq. 1) then                 ! i- 面
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5 
				   B%U(m,0,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 2) then
 				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,0,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 3) then             
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,0)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
  
             else if(Bc%face .eq. 4) then 
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5
				   B%U(m,ie,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo

             else if(Bc%face .eq. 5) then
				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,je,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
   
             else if(Bc%face .eq. 6) then
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,ke)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
             endif

     end 
!----------------------------------------------------------------------------------------

!------------------用户自定义的边界条件-----------------------------------
! 给定入口时间序列，文件名 inlet-time.xxx  (xxx为块号)
! 仅使用1层虚网格

    subroutine boundary_user_Inlet_time(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc

     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1
     integer:: n1,n2,n3,Istep_inlet
	 integer,save:: read_flag=0
     real(PRE_EC),allocatable,save,dimension(:,:,:,:):: uc,fc       ! 边界值 
     real(PRE_EC):: tt_inlet
	 character(len=30):: fname

     NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je; kb=Bc%kb; ke=Bc%ke    
     n1=max(ie-ib,1)
	 n2=max(je-jb,1)
	 n3=max(ke-kb,1)

 !-----读取数据--------------------------------------   
	 if(read_flag==0) then
	  read_flag=1
      allocate(uc(n1,n2,n3,5),fc(5,n1,n2,n3))
      write(fname,"('inlet-time.',I3.3)") B%Block_no
	  open(999,file=fname,form="unformatted")
	  print*, "read file : ", fname
	  print*, "n1,n2,n3=",n1,n2,n3
    endif

! 初始时刻 (KRK=0) 以及 更新时间步  (KRK=3) 时， 读取 uc数据
! 如采用1步 Euler法， 则KRK保持0

   if(KRK ==0 .or. KRK==3) then
      read(999) Istep_inlet,tt_inlet
      print*, "Istep_inlet,tt_inlet=", Istep_inlet, tt_inlet
	  do m=1,5
	   read(999) (((uc(i,j,k,m),i=1,n1),j=1,n2),k=1,n3)
      enddo

        do k=1,n3
		do j=1,n2
		do i=1,n1
          fc(1,i,j,k)=uc(i,j,k,1)
          fc(2,i,j,k)=uc(i,j,k,1)*uc(i,j,k,2)
          fc(3,i,j,k)=uc(i,j,k,1)*uc(i,j,k,3)
          fc(4,i,j,k)=uc(i,j,k,1)*uc(i,j,k,4)
          fc(5,i,j,k)=uc(i,j,k,1)*(Cv*uc(i,j,k,5)+0.5d0*(uc(i,j,k,2)**2+uc(i,j,k,3)**2+uc(i,j,k,4)**2))
        enddo
	    enddo
	    enddo
    endif

!-----------------------------------------------------------------

 ! 仅使用1层Ghost Cell    
             if(Bc%face .eq. 1) then                 ! i- 面
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5 
				   B%U(m,0,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 2) then
 				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,0,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
             else if(Bc%face .eq. 3) then             
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,0)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
  
             else if(Bc%face .eq. 4) then 
 				 do k=kb,ke-1
				 do j=jb,je-1
                 do m=1,5
				   B%U(m,ie,j,k)=fc(m,1,j,k)
				 enddo
				 enddo
				 enddo

             else if(Bc%face .eq. 5) then
				 do k=kb,ke-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,je,k)=fc(m,i,1,k)
				 enddo
				 enddo
				 enddo
   
             else if(Bc%face .eq. 6) then
				 do j=jb,je-1
                 do i=ib,ie-1				   
                 do m=1,5 
				   B%U(m,i,j,ke)=fc(m,i,j,1)
				 enddo
				 enddo
				 enddo
             endif

     end 
!----------------------------------------------------------------------------------------


!-------------------------------------------------------------------  
! 包含吹吸扰动的壁面 （扰动范围x-z方向）
! 仅使用1层虚网格


    subroutine boundary_user_blow_suction_Wall(nMesh,mBlock,ksub)
     Use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: mBlock,ksub,ib,ie,jb,je,kb,ke,i,j,k,m,nMesh,NVAR1,i1,j1,k1,i2,j2,k2,nk,ierr
     real(PRE_EC):: n1,n2,n3,Vn
	 integer,save:: Iflag=0
     real(PRE_EC):: d1,u1,v1,w1,T1,p1,d2,u2,v2,w2,T2,p2,xb,xe,zb,ze,wx,wz,Aw
     real(PRE_EC),save:: rpara(100)
! --------- 仅在第1次调用本子程序时执行 ------------- 
!   读取辅助参数 (供自定义边界条件使用)
   if(Iflag==0) then
     Iflag=1
        open(99,file="bc_user.ec")
		read(99,*) 
		read(99,*)
		read(99,*) nk
		read(99,*)
		read(99,*) (rpara(k),k=1,nk)
        close(99)
      print*, "my_id=",  my_id, "rpara=", (rpara(k),k=1,nk)
   endif
! --------------------------------------------
	 
	 NVAR1=Mesh(nMesh)%NVAR
     B => Mesh(nMesh)%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      

    if(NVAR1 > 5) then
	  print*, "Warning, NVAR1 should =5 !!!"
	  print*, "if you using blow-and-suction wall, you should not use Turbulence model"
	endif

! 壁面吹吸扰动区域、波数及振幅    
	xb=rpara(1)
	xe=rpara(2)
	zb=rpara(3)
	ze=rpara(4)
	wx=rpara(5)
	wz=rpara(6)
	Aw=rpara(7)
!-------------------------

!   仅使用1层Ghost Cell
!   i2 Ghost Cell;  i1 内点
     if(Bc%face .eq. 1 .or. Bc%face .eq. 4) then   ! i- or i+
       if(Bc%face .eq. 1) then
         i=ib; i1=ib; i2=ib-1
       else
         i=ie; i1=ie-1; i2=ie 
       endif

       do k=kb,ke-1
         do j=jb,je-1
    

		    d1=B%U(1,i1,j,k); u1=B%U(2,i1,j,k)/d1; v1=B%U(3,i1,j,k)/d1; w1=B%U(4,i1,j,k)/d1
            p1=(B%U(5,i1,j,k)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)      
            T1=gamma*Ma*Ma*p1/d1 
			p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            
			if(Twall >0) then
		       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
		    else
			   T2=T1
		    endif
	        d2=gamma*Ma*Ma*p2/T2 
		
		   if( B%xc(i1,j,k) >= xb  .and. B%xc(i1,j,k) <=xe .and. B%zc(i1,j,k) >= zb .and. B%zc(i1,j,k) <=ze) then
		     Vn= Aw*sin( 2.d0*PI*wx*(B%xc(i1,j,k)-xb)/(xe-xb)) *sin ( 2.d0*PI*wz*(B%zc(i1,j,k)-zb)/(ze-zb))
           else
		     Vn=0.d0
		   endif

		    n1=B%ni1(i,j,k) ; n2=B%ni2(i,j,k); n3=B%ni3(i,j,k)   ! 归一化法方向  
            u2= 2.d0*Vn*n1-u1
            v2= 2.d0*Vn*n2-v1
			w2= 2.d0*Vn*n3-w1
			          
           B%U(1,i2,j,k)= d2
           B%U(2,i2,j,k)= d2*u2
           B%U(3,i2,j,k)= d2*v2
           B%U(4,i2,j,k)= d2*w2
           B%U(5,i2,j,k)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)
		 
		 enddo
       enddo

     else if(Bc%face .eq. 2  .or. Bc%face .eq. 5) then   ! j- or j+
       if(Bc%face .eq. 2) then
         j=jb; j1=jb; j2=jb-1 
       else
         j=je; j1=je-1; j2=je 
       endif

       do k=kb,ke-1
         do i=ib,ie-1
 
 		    d1=B%U(1,i,j1,k); u1=B%U(2,i,j1,k)/d1; v1=B%U(3,i,j1,k)/d1; w1=B%U(4,i,j1,k)/d1
            p1=(B%U(5,i,j1,k)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)      
            T1=gamma*Ma*Ma*p1/d1 
    		p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            
			if(Twall >0) then
		       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
		    else
			   T2=T1
		    endif
	        d2=gamma*Ma*Ma*p2/T2 
		
		   if( B%xc(i,j1,k) >= xb  .and. B%xc(i,j1,k) <=xe .and. B%zc(i,j1,k) >= zb .and. B%zc(i,j1,k) <=ze) then
		     Vn= Aw*sin( 2.d0*PI*wx*(B%xc(i,j1,k)-xb)/(xe-xb)) *sin ( 2.d0*PI*wz*(B%zc(i,j1,k)-zb)/(ze-zb))
           else
		     Vn=0.d0
		   endif

            n1=B%nj1(i,j,k) ; n2=B%nj2(i,j,k); n3=B%nj3(i,j,k)     
            u2= 2.d0*Vn*n1-u1
            v2= 2.d0*Vn*n2-v1
			w2= 2.d0*Vn*n3-w1
			          
           B%U(1,i,j2,k)= d2
           B%U(2,i,j2,k)= d2*u2
           B%U(3,i,j2,k)= d2*v2
           B%U(4,i,j2,k)= d2*w2
           B%U(5,i,j2,k)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)

          enddo
       enddo

     else if(Bc%face .eq. 3 .or. Bc%face .eq. 6) then   ! k- or k+
       if(Bc%face .eq. 3) then
         k=kb; k1=kb; k2=kb-1 
       else
         k=ke; k1=ke-1; k2=ke
       endif

       do j=jb,je-1
         do i=ib,ie-1

 		    d1=B%U(1,i,j,k1); u1=B%U(2,i,j,k1)/d1; v1=B%U(3,i,j,k1)/d1; w1=B%U(4,i,j,k1)/d1
            p1=(B%U(5,i,j,k1)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))*(gamma-1.d0)      
            T1=gamma*Ma*Ma*p1/d1 
    		p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            
			if(Twall >0) then
		       T2=2.d0*Twall-T1    ! 等温壁，温度外插    0.5*(T1+T2)=Twall
		    else
			   T2=T1
		    endif
	        d2=gamma*Ma*Ma*p2/T2 
		
		   if( B%xc(i,j,k1) >= xb  .and. B%xc(i,j,k1) <=xe .and. B%zc(i,j,k1) >= zb .and. B%zc(i,j,k1) <=ze) then
		     Vn= Aw*sin( 2.d0*PI*wx*(B%xc(i,j,k1)-xb)/(xe-xb)) *sin ( 2.d0*PI*wz*(B%zc(i,j,k1)-zb)/(ze-zb))
           else
		     Vn=0.d0
		   endif

		    n1=B%nk1(i,j,k) ; n2=B%nk2(i,j,k); n3=B%nk3(i,j,k)     
            u2= 2.d0*Vn*n1-u1
            v2= 2.d0*Vn*n2-v1
			w2= 2.d0*Vn*n3-w1
			          
           B%U(1,i,j,k2)= d2
           B%U(2,i,j,k2)= d2*u2
           B%U(3,i,j,k2)= d2*v2
           B%U(4,i,j,k2)= d2*w2
           B%U(5,i,j,k2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2+w2*w2)
		 enddo
       enddo

     endif
        
    end 





  end module mod_struct_bc
