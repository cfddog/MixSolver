!===============================================================================
! mod_struct_grid.f90 -- structured solver: grid geometry / format conversion
! Encapsulates sub_convert_inp + sub_geometry + sub_debug.
! The legacy Type_def1 helper module is kept verbatim above mod_struct_grid;
! its BC_MSG_TYPE shadows Mod_Type_Def's, same as in the pristine build.
!===============================================================================
!----Boundary message (bc3d.inp, Gridgen general format)----------------------------------------------------------
! 读取边界格式文件(bc3d.inp)，并转化为OpenCFD-EC的内建.inc格式 (见《OpenCFD-EC理论手册》 2.4.3节)
! OpenCFD-EC内建的存储格式比bc3d.inp多了一些冗余信息 （有些类似BXCFD的.in格式），例如多了子面号f_no,
! 面类型face 以及连接的子面号f_no1,连接的面类型face1
! 以及连接次序L1, L2, L3  (例如L1=1表示该维与连接块的第1维正连接， L1=-1表示与连接块的第1为反向连接).
! 这些冗余信息为块-块之间的通信（尤其是MPI并行通信）提供了便利，有利于简化代码
 
 module  Type_def1
    TYPE BC_MSG_TYPE              ! 边界链接信息
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! 边界区域（子面）的定义， .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! 连接区域
	 integer:: L1,L2,L3                     ! 子面号，连接顺序描述符
   END TYPE BC_MSG_TYPE


    TYPE Block_TYPE1                          ! 数据结构：仅包含Bc_msg 
	 integer::  nx,ny,nz                      ! 网格数nx,ny,nz
	 integer::  subface                       ! 子面数
 	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! 边界链接信息 
    END TYPE Block_TYPE1   
 End module  Type_def1

 module mod_struct_grid

 contains

!------------------------------------------------------
  subroutine convert_inp_inc
   use Type_def1
   implicit none
  
   integer:: NB,m,ksub,nx,ny,nz,k,j,k1,ksub1
   integer:: kb(3),ke(3),kb1(3),ke1(3),s(3),p(3),Lp(3)
   TYPE(Block_TYPE1),Pointer,dimension(:):: Block
   Type (Block_TYPE1),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc,Bc1
     
   

   print*, "Convert bc3d.inp to bc3d.inc ..."
   open(88,file="bc3d.inp")
   read(88,*)
   read(88,*) NB
   allocate(Block(NB))
   do m=1,NB
    B => Block(m)
    read(88,*) B%nx,B%ny,B%nz
	read(88,*)
    read(88,*) B%subface   !number of the subface in the Block m
     
	allocate(B%bc_msg(B%subface))   ! 边界描述

    do ksub=1, B%subface
      Bc => B%bc_msg(ksub)
      Bc%f_no=ksub                        ! 子面号
	  read(88,*)  kb(1),ke(1),kb(2),ke(2),kb(3),ke(3),Bc%bc
	 
	  if(Bc%bc .lt. 0) then
 !  --------有连接的情况 (内边界)--------------------------------------------------------	    
	   read(88,*) kb1(1),ke1(1),kb1(2),ke1(2),kb1(3),ke1(3),Bc%nb1
      else
 !---------无连接情况 (物理边界)----------------------------------------
            kb1(:)=0; ke1(:)=0; Bc%nb1=0
	  endif
	 call Convert_bc(Bc,kb,ke,kb1,ke1)
  enddo
  enddo
   
   close(88)

!  搜索连接块的块号 f_no1  (便于MPI并行通信是使用)
   do m=1,NB
     B => Block(m)
     do ksub=1, B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .lt. 0) then
         Bc%f_no1=0
		 B1=>Block(Bc%nb1)         ! 指向连接块
         
		 do ksub1=1,B1%subface
		 Bc1=>B1%bc_msg(ksub1)
		 if(Bc%ib1==Bc1%ib .and. Bc%ie1==Bc1%ie .and. Bc%jb1==Bc1%jb .and. Bc%je1==Bc1%je   &
		    .and. Bc%kb1==Bc1%kb .and. Bc%ke1==Bc1%ke  .and. Bc1%nb1==m) then
		    Bc%f_no1=ksub1
		  exit
		 endif 	 
         enddo

         if(Bc%f_no1 ==0) then
		  print*, " Error in find linked block number !!!"
		  print*, "Block, subface=",m, ksub
		  stop
		 endif
	   endif
      enddo
	enddo

   open(99,file="bc3d.inc")
    write(99,*) " Inp-liked file of OpenCFD-EC"
    write(99,*) NB
	do m=1,NB
     B => Block(m)
     write(99,*) B%nx, B%ny, B%nz
	 write(99,*) "Block ", m
	 write(99,*) B%subface
	  do ksub=1,B%subface
        Bc => B%bc_msg(ksub)
	   write(99,"(9I6)") Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
	   write(99,"(12I6)") Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
      enddo
	enddo
	close(99)

   print*, "Convert bc3d.inp to bc3d.inc OK"


  end  subroutine convert_inp_inc

!-----------------------------------------------------------------------------------




!   将Gridgen格式 转换为 OpenCFD-EC 的边界连接格式
!     计算面号、 连接次序 (L1,L2,L3)等
      subroutine Convert_bc(Bc,kb,ke,kb1,ke1)
       use Type_Def1
       implicit none
       TYPE (BC_MSG_TYPE),pointer:: Bc
	   integer,dimension(3):: kb,ke,kb1,ke1,s,p,LP
       integer:: k,j,k1

	   if(Bc%bc .ge. 0) then
	     Bc%ib1=0; Bc%ie1=0; Bc%jb1=0; Bc%je1=0; Bc%kb1=0; Bc%ke1=0; Bc%nb1=0
         Bc%L1=0; Bc%L2=0; Bc%L3=0; Bc%face1=0; Bc%f_no1=0
       endif
     

!   判断该面的类型 (i-, i+, j-,j+, k-,k+)     
       do k=1,3
  	     if(kb(k) .eq. ke(k) ) then 
	       s(k)=0                           ! 连接维
	     else if (kb(k) .gt. 0) then 
	       s(k)=1                           ! 正
	     else
	       s(k)=-1                          ! 负
	     endif
       enddo

!    边界子面的大小     
	 Bc%ib=min(abs(kb(1)),abs(ke(1))) ;  Bc%ie=max(abs(kb(1)),abs(ke(1)))
     Bc%jb=min(abs(kb(2)),abs(ke(2))) ;  Bc%je=max(abs(kb(2)),abs(ke(2)))
     Bc%kb=min(abs(kb(3)),abs(ke(3))) ;  Bc%ke=max(abs(kb(3)),abs(ke(3))) 


!   判断该面的类型 (i-, i+, j-,j+, k-,k+)     
      if(s(1) .eq. 0) then
	     if (Bc%ib .eq. 1) then
	      Bc%face=1               ! i-
	     else
	      Bc%face=4               ! i+
	    endif
      else if(s(2) .eq. 0) then
	    if(Bc%jb .eq. 1) then
	      Bc%face=2                 ! j-
	    else
	      Bc%face=5                 ! j+
	    endif
      else
 	   if(Bc%kb  .eq. 1) then
	    Bc%face=3                 ! k-
	   else
	    Bc%face=6                 ! k+
	   endif
     endif 

!---------------------------------------------------------------------------
!------内边界的情况，建立连接描述
  if( Bc%bc .lt. 0) then            ! 内边界
!      计算连接顺序描述符L1,L2,L3
!      计算各维之间的连接关系      
     do k=1,3  
	   if(kb1(k) .eq. ke1(k) ) then 
	       p(k)=0                      ! 
       else if (kb1(k) .gt. 0) then
	       p(k)=1                      ! .inp 文件的 正数
       else
	       p(k)=-1
       endif
     enddo
 	   

!    对应连接子面的大小     
	 Bc%ib1=min(abs(kb1(1)),abs(ke1(1))) ;  Bc%ie1=max(abs(kb1(1)),abs(ke1(1)))
     Bc%jb1=min(abs(kb1(2)),abs(ke1(2))) ;  Bc%je1=max(abs(kb1(2)),abs(ke1(2)))
     Bc%kb1=min(abs(kb1(3)),abs(ke1(3))) ;  Bc%ke1=max(abs(kb1(3)),abs(ke1(3))) 
    
 	  
!   判断该面连接面的类型 (i-, i+, j-,j+, k-,k+)     
      if(p(1) .eq. 0) then
	     if (Bc%ib1 .eq. 1) then
	      Bc%face1=1               ! i-
	     else
	      Bc%face1=4               ! i+
	    endif
      else if(p(2) .eq. 0) then
	    if(Bc%jb1 .eq. 1) then
	      Bc%face1=2                 ! j-
	    else
	      Bc%face1=5                 ! j+
	    endif
      else
 	   if(Bc%kb1 .eq. 1) then
	    Bc%face1=3
	   else
	    Bc%face1=6
	   endif
     endif  

!  计算“连接对” 描述符  bc%L1, bc%L2, bc%L3 
 	   do k=1,3
	     do j=1,3
	       if(s(k) .eq. p(j)) Lp(k)=j          ! .inp文件的连接格式： 正对正、 负对负、 0对0； 
	     enddo
	   enddo    
	   
!    计算连接次序 （正为顺序；负为拟序）	  
	  do k=1,3
	   if(s(k) .ne. 0) then
	     k1=Lp(k)
	     if( (ke(k)-kb(k))*(ke1(k1)-kb1(k1)) .lt. 0) Lp(k)=-Lp(k)      ! 逆序连接
	   else
         k1=Lp(k)
!		 if( (mod(Bc%face,2)-mod(Bc%face1,2))==0) Lp(k)=-Lp(k)         ! 逆序连接 （单面） ! Bug 2012-7-13
! 正-正连接 Lp为负 (例， i+ 面连接到 j+面， 则为逆序连接)
		 if( (Bc%face-1)/3 .eq. (Bc%face1-1)/3 ) Lp(k)=-Lp(k)        ! (Bc%face=1,2,3为 +面，4,5,6为-面)  ! 逆序连接 （单面）

	   endif
	 enddo 
     
	 Bc%L1=Lp(1); Bc%L2=Lp(2); Bc%L3=Lp(3)   ! 连接次序描述符 （详见《理论手册》）
   endif
  end subroutine Convert_bc

!-----------------------------------------------------
!   计算几何量：控制体体积和中心点坐标，Jocabian系数
!   2013-4-26:  可处理退化线 （面积为0的面）
!   2013-5-3: 修改粘性项Jocabian系数计算方法，与物理量导数方法一致

  subroutine Comput_Goemetric_var(nMesh)
   use   Global_Var
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nMesh
   real(PRE_EC):: t1x,t1y,t1z,t2x,t2y,t2z,s1x,s1y,s1z,xa,ya,za,ss
   real(PRE_EC):: xi,yi,zi,xj,yj,zj,xk,yk,zk,Jac,Jac1
   real(PRE_EC):: xi1,xi2,yi1,yi2,zi1,zi2,xj1,xj2,yj1,yj2,zj1,zj2,xk1,xk2,yk1,yk2,zk1,zk2
   real(PRE_EC),allocatable,dimension(:,:,:)::Vi,Vj,Vk
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
!  计算控制体的体积
!  计算控制体各表面的面积 （为了避免内存占用过多，表面的法方向、切方向在计算中求出，不进行存储）
   MP=>Mesh(nMesh)
   do m=1,MP%Num_Block   
     B => MP%Block(m)
     allocate(Vi(B%nx,B%ny,B%nz),Vj(B%nx,B%ny,B%nz),Vk(B%nx,B%ny,B%nz))
!  Area of surface i, j and k     
     do k=1,B%nz
       do j=1,B%ny
         do i=1,B%nx
!-------------------------------------------
           t1x=B%x(i,j+1,k)-B%x(i,j,k+1); t1y=B%y(i,j+1,k)-B%y(i,j,k+1); t1z=B%z(i,j+1,k)-B%z(i,j,k+1)   ! 对角线1
           t2x=B%x(i,j+1,k+1)-B%x(i,j,k); t2y=B%y(i,j+1,k+1)-B%y(i,j,k) ; t2z=B%z(i,j+1,k+1)-B%z(i,j,k)  ! 对角线2
           s1x=t1y*t2z-t1z*t2y ; s1y=t1z*t2x-t1x*t2z ; s1z=t1x*t2y-t1y*t2x   ! 法向量 （对角线向量叉乘得到）
           ss=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)  ! 长度
           B%Si(i,j,k)=ss*0.5d0
           
		   if(ss .ge. Lim_Zero) then
		     B%ni1(i,j,k)=s1x/ss; B%ni2(i,j,k)=s1y/ss ; B%ni3(i,j,k)=s1z/ss
           else
		     B%ni1(i,j,k)=1.d0; B%ni2(i,j,k)=0.d0 ; B%ni3(i,j,k)=0.d0          ! 退化线；任意确定法方向
		   endif
		     

           xa=(B%x(i,j,k)+B%x(i,j+1,k)+B%x(i,j,k+1)+B%x(i,j+1,k+1))*0.25d0
           ya=(B%y(i,j,k)+B%y(i,j+1,k)+B%y(i,j,k+1)+B%y(i,j+1,k+1))*0.25d0
           za=(B%z(i,j,k)+B%z(i,j+1,k)+B%z(i,j,k+1)+B%z(i,j+1,k+1))*0.25d0
           Vi(i,j,k)=(s1x*xa+s1y*ya+s1z*za)*0.5d0
!----------------------------------------------------------------------------------
           t1x=B%x(i+1,j,k+1)-B%x(i,j,k); t1y=B%y(i+1,j,k+1)-B%y(i,j,k) ; t1z=B%z(i+1,j,k+1)-B%z(i,j,k)  ! 对角线1
           t2x=B%x(i+1,j,k)-B%x(i,j,k+1); t2y=B%y(i+1,j,k)-B%y(i,j,k+1);  t2z=B%z(i+1,j,k)-B%z(i,j,k+1)   ! 对角线2
           s1x=t1y*t2z-t1z*t2y ; s1y=t1z*t2x-t1x*t2z ; s1z=t1x*t2y-t1y*t2x   ! 法向量 （对角线向量叉乘得到）
           ss=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)
           B%Sj(i,j,k)=ss*0.5d0
    
	   	   if(ss .ge. Lim_Zero) then
    	     B%nj1(i,j,k)=s1x/ss; B%nj2(i,j,k)=s1y/ss ; B%nj3(i,j,k)=s1z/ss
           else
    	     B%nj1(i,j,k)=1.d0; B%nj2(i,j,k)=0.d0 ; B%nj3(i,j,k)=0.d0
		   endif


           xa=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j,k+1)+B%x(i+1,j,k+1))*0.25d0
           ya=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j,k+1)+B%y(i+1,j,k+1))*0.25d0
           za=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j,k+1)+B%z(i+1,j,k+1))*0.25d0
           Vj(i,j,k)=(s1x*xa+s1y*ya+s1z*za)*0.5d0
!----------------------=----------------------------------------------------
           t1x=B%x(i+1,j+1,k)-B%x(i,j,k); t1y=B%y(i+1,j+1,k)-B%y(i,j,k) ; t1z=B%z(i+1,j+1,k)-B%z(i,j,k)  ! 对角线1
           t2x=B%x(i,j+1,k)-B%x(i+1,j,k); t2y=B%y(i,j+1,k)-B%y(i+1,j,k) ; t2z=B%z(i,j+1,k)-B%z(i+1,j,k)   ! 对角线2
           s1x=t1y*t2z-t1z*t2y ; s1y=t1z*t2x-t1x*t2z ; s1z=t1x*t2y-t1y*t2x   ! 法向量 （对角线向量叉乘得到）
           ss=sqrt(s1x*s1x+s1y*s1y+s1z*s1z)
           B%Sk(i,j,k)=ss*0.5d0  
           if(ss .ge. Lim_Zero) then
             B%nk1(i,j,k)=s1x/ss; B%nk2(i,j,k)=s1y/ss; B%nk3(i,j,k)=s1z/ss
           else
             B%nk1(i,j,k)=1.d0; B%nk2(i,j,k)=0.d0; B%nk3(i,j,k)=0.d0
		   endif

           xa=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j+1,k)+B%x(i+1,j+1,k))*0.25d0
           ya=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j+1,k)+B%y(i+1,j+1,k))*0.25d0
           za=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j+1,k)+B%z(i+1,j+1,k))*0.25d0
           Vk(i,j,k)=(s1x*xa+s1y*ya+s1z*za)*0.5d0
!-----------------------------------------------
         enddo
       enddo
     enddo

! 控制体 体积    
     do k=1,B%nz-1
       do j=1,B%ny-1
         do i=1,B%nx-1
           B%vol(i,j,k)=(Vi(i+1,j,k)-Vi(i,j,k)+Vj(i,j+1,k)-Vj(i,j,k)+Vk(i,j,k+1)-Vk(i,j,k))/3.d0
         enddo
       enddo
     enddo

!  网格中心点坐标
     do k=0,B%nz
       do j=0,B%ny
         do i=0,B%nx
           B%xc(i,j,k)=(B%x(i,j,k)+B%x(i,j+1,k)+B%x(i,j,k+1)+B%x(i,j+1,k+1)+ &
                       B%x(i+1,j,k)+B%x(i+1,j+1,k)+B%x(i+1,j,k+1)+B%x(i+1,j+1,k+1))*0.125
           B%yc(i,j,k)=(B%y(i,j,k)+B%y(i,j+1,k)+B%y(i,j,k+1)+B%y(i,j+1,k+1)+ &
                       B%y(i+1,j,k)+B%y(i+1,j+1,k)+B%y(i+1,j,k+1)+B%y(i+1,j+1,k+1))*0.125
           B%zc(i,j,k)=(B%z(i,j,k)+B%z(i,j+1,k)+B%z(i,j,k+1)+B%z(i,j+1,k+1)+ &
                       B%z(i+1,j,k)+B%z(i+1,j+1,k)+B%z(i+1,j,k+1)+B%z(i+1,j+1,k+1))*0.125
         enddo
       enddo
     enddo              
     deallocate(Vi,Vj,Vk)
 ! -----------计算 Jocabian系数 （粘性项计算导数时使用）------------------------
 !  (I+1/2,J,K)点
 ! revised, 2013-5-3:  Jocabian系数与 物理量导数 计算方法一致，避免额外误差；
 ! revised, 2013-5-4:  避免使用角点坐标（计算域立方体的棱），以免出现不稳定性
       do k=1,B%nz-1 
       do j=1,B%ny-1
       do i=1,B%nx
        
		  xi=B%xc(i,j,k)-B%xc(i-1,j,k) 
          yi=B%yc(i,j,k)-B%yc(i-1,j,k)
          zi=B%zc(i,j,k)-B%zc(i-1,j,k)
		 

          if( (i==1 .or. i==B%nx) .and. (j==1 .or. j==B%ny-1) ) then
		   xj1=B%xc(i,j-1,k)
		   yj1=B%yc(i,j-1,k)
		   zj1=B%zc(i,j-1,k)
		   xj2=B%xc(i,j+1,k)
		   yj2=B%yc(i,j+1,k)
		   zj2=B%zc(i,j+1,k)
          else
		   xj1=0.5d0*(B%xc(i,j-1,k)+B%xc(i-1,j-1,k))
		   yj1=0.5d0*(B%yc(i,j-1,k)+B%yc(i-1,j-1,k))
		   zj1=0.5d0*(B%zc(i,j-1,k)+B%zc(i-1,j-1,k))
		   xj2=0.5d0*(B%xc(i,j+1,k)+B%xc(i-1,j+1,k))
		   yj2=0.5d0*(B%yc(i,j+1,k)+B%yc(i-1,j+1,k))
		   zj2=0.5d0*(B%zc(i,j+1,k)+B%zc(i-1,j+1,k))
          endif

          if( (i==1 .or. i==B%nx) .and. (k==1 .or. k==B%nz-1) ) then
            xk1=B%xc(i,j,k-1)
            yk1=B%yc(i,j,k-1)
            zk1=B%zc(i,j,k-1)
            xk2=B%xc(i,j,k+1)
            yk2=B%yc(i,j,k+1)
            zk2=B%zc(i,j,k+1)
		  else
            xk1=0.5d0*(B%xc(i,j,k-1)+B%xc(i-1,j,k-1))
            yk1=0.5d0*(B%yc(i,j,k-1)+B%yc(i-1,j,k-1))
            zk1=0.5d0*(B%zc(i,j,k-1)+B%zc(i-1,j,k-1))
            xk2=0.5d0*(B%xc(i,j,k+1)+B%xc(i-1,j,k+1))
            yk2=0.5d0*(B%yc(i,j,k+1)+B%yc(i-1,j,k+1))
            zk2=0.5d0*(B%zc(i,j,k+1)+B%zc(i-1,j,k+1))
          endif
	      
 		   xj=0.5d0*(xj2-xj1)
		   yj=0.5d0*(yj2-yj1)
		   zj=0.5d0*(zj2-zj1)
 		   xk=0.5d0*(xk2-xk1)
		   yk=0.5d0*(yk2-yk1)
		   zk=0.5d0*(zk2-zk1)

           Jac1=(xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk)
           if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
		   Jac=1.d0/Jac1

!   Jac=B%Jaci(i,j,k)
!   9个Jocabian变换系数    
          B%ix1(i,j,k)=Jac*(yj*zk-zj*yk)
          B%iy1(i,j,k)=Jac*(zj*xk-xj*zk)
          B%iz1(i,j,k)=Jac*(xj*yk-yj*xk)
          B%jx1(i,j,k)=Jac*(yk*zi-zk*yi)
          B%jy1(i,j,k)=Jac*(zk*xi-xk*zi)
          B%jz1(i,j,k)=Jac*(xk*yi-yk*xi)
          B%kx1(i,j,k)=Jac*(yi*zj-zi*yj)
          B%ky1(i,j,k)=Jac*(zi*xj-xi*zj)
          B%kz1(i,j,k)=Jac*(xi*yj-yi*xj)
      enddo
	  enddo
	  enddo
 
 ! (I,J-1/2,K) 点的值， 即 (i+1/2,j,k+1/2)点的值
 
! Revised, 2013-5-3, 坐标的导数与物理量的导数 计算方法相同
      do k=1,B%nz-1 
      do j=1,B%ny
      do i=1,B%nx-1
       xj=B%xc(i,j,k)-B%xc(i,j-1,k)
       yj=B%yc(i,j,k)-B%yc(i,j-1,k)
       zj=B%zc(i,j,k)-B%zc(i,j-1,k)
      
! Revised, 2013-5-4, 避免使用角点（棱）坐标
	   if( (j==1 .or. j==B%ny) .and. (i==1 .or. i==B%nx-1) ) then
		xi1=B%xc(i-1,j,k)
		yi1=B%yc(i-1,j,k)
		zi1=B%zc(i-1,j,k)
		xi2=B%xc(i+1,j,k)
		yi2=B%yc(i+1,j,k)
		zi2=B%zc(i+1,j,k)
	   else
		xi1=0.5d0*(B%xc(i-1,j,k)+B%xc(i-1,j-1,k))
		yi1=0.5d0*(B%yc(i-1,j,k)+B%yc(i-1,j-1,k))
		zi1=0.5d0*(B%zc(i-1,j,k)+B%zc(i-1,j-1,k))
		xi2=0.5d0*(B%xc(i+1,j,k)+B%xc(i+1,j-1,k))
		yi2=0.5d0*(B%yc(i+1,j,k)+B%yc(i+1,j-1,k))
		zi2=0.5d0*(B%zc(i+1,j,k)+B%zc(i+1,j-1,k))
       endif
	
	   if( (j==1 .or. j==B%ny) .and. (k==1 .or. k==B%nz-1) ) then
		 xk1=B%xc(i,j,k-1)
		 yk1=B%yc(i,j,k-1)
		 zk1=B%zc(i,j,k-1)
		 xk2=B%xc(i,j,k+1)
		 yk2=B%yc(i,j,k+1)
		 zk2=B%zc(i,j,k+1)
	   else	 
		 xk1=0.5d0*(B%xc(i,j,k-1)+B%xc(i,j-1,k-1))
		 yk1=0.5d0*(B%yc(i,j,k-1)+B%yc(i,j-1,k-1))
		 zk1=0.5d0*(B%zc(i,j,k-1)+B%zc(i,j-1,k-1))
		 xk2=0.5d0*(B%xc(i,j,k+1)+B%xc(i,j-1,k+1))
		 yk2=0.5d0*(B%yc(i,j,k+1)+B%yc(i,j-1,k+1))
		 zk2=0.5d0*(B%zc(i,j,k+1)+B%zc(i,j-1,k+1))
       endif
	    xi=0.5d0*(xi2-xi1)
	    yi=0.5d0*(yi2-yi1)
	    zi=0.5d0*(zi2-zi1)
	    xk=0.5d0*(xk2-xk1)
	    yk=0.5d0*(yk2-yk1)
	    zk=0.5d0*(zk2-zk1)
	
      
	  Jac1=xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk
      if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
	  Jac=1.d0/Jac1
	 
	  B%ix2(i,j,k)=Jac*(yj*zk-zj*yk)
      B%iy2(i,j,k)=Jac*(zj*xk-xj*zk)
      B%iz2(i,j,k)=Jac*(xj*yk-yj*xk)
      B%jx2(i,j,k)=Jac*(yk*zi-zk*yi)
      B%jy2(i,j,k)=Jac*(zk*xi-xk*zi)
      B%jz2(i,j,k)=Jac*(xk*yi-yk*xi)
      B%kx2(i,j,k)=Jac*(yi*zj-zi*yj)
      B%ky2(i,j,k)=Jac*(zi*xj-xi*zj)
      B%kz2(i,j,k)=Jac*(xi*yj-yi*xj)
	 enddo
	 enddo
	 enddo
	
! (I,J,K-1/2) 点的值， 即 (i+1/2,j+1/2,k)点的值
! Revised, 2013-5-3 
 
     do k=1,B%nz 
     do j=1,B%ny-1
     do i=1,B%nx-1
      xk=B%xc(i,j,k)-B%xc(i,j,k-1)
      yk=B%yc(i,j,k)-B%yc(i,j,k-1)
      zk=B%zc(i,j,k)-B%zc(i,j,k-1)
 
 
 	 if( (k==1 .or. k==B%nz) .and. (i==1 .or. i==B%nx-1) ) then
	  xi1=B%xc(i-1,j,k)
	  yi1=B%yc(i-1,j,k)
	  zi1=B%zc(i-1,j,k)
	  xi2=B%xc(i+1,j,k)
	  yi2=B%yc(i+1,j,k)
	  zi2=B%zc(i+1,j,k)
     else
	  xi1=0.5d0*(B%xc(i-1,j,k)+B%xc(i-1,j,k-1))
	  yi1=0.5d0*(B%yc(i-1,j,k)+B%yc(i-1,j,k-1))
	  zi1=0.5d0*(B%zc(i-1,j,k)+B%zc(i-1,j,k-1))
	  xi2=0.5d0*(B%xc(i+1,j,k)+B%xc(i+1,j,k-1))
	  yi2=0.5d0*(B%yc(i+1,j,k)+B%yc(i+1,j,k-1))
	  zi2=0.5d0*(B%zc(i+1,j,k)+B%zc(i+1,j,k-1))
     endif

  	 if( (k==1 .or. k==B%nz) .and. (j==1 .or. j==B%ny-1) ) then
      xj1=B%xc(i,j-1,k)
      yj1=B%yc(i,j-1,k)
      zj1=B%zc(i,j-1,k)
      xj2=B%xc(i,j+1,k)
      yj2=B%yc(i,j+1,k)
      zj2=B%zc(i,j+1,k)
	 else
      xj1=0.5d0*(B%xc(i,j-1,k)+B%xc(i,j-1,k-1))
      yj1=0.5d0*(B%yc(i,j-1,k)+B%yc(i,j-1,k-1))
      zj1=0.5d0*(B%zc(i,j-1,k)+B%zc(i,j-1,k-1))
      xj2=0.5d0*(B%xc(i,j+1,k)+B%xc(i,j+1,k-1))
      yj2=0.5d0*(B%yc(i,j+1,k)+B%yc(i,j+1,k-1))
      zj2=0.5d0*(B%zc(i,j+1,k)+B%zc(i,j+1,k-1))
	 endif
	  xi=0.5d0*(xi2-xi1)
	  yi=0.5d0*(yi2-yi1)
	  zi=0.5d0*(zi2-zi1)
	  xj=0.5d0*(xj2-xj1)
	  yj=0.5d0*(yj2-yj1)
	  zj=0.5d0*(zj2-zj1)
    
	 

	  Jac1=xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk
      if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
	  Jac=1.d0/Jac1
	 
	 B%ix3(i,j,k)=Jac*(yj*zk-zj*yk)
     B%iy3(i,j,k)=Jac*(zj*xk-xj*zk)
     B%iz3(i,j,k)=Jac*(xj*yk-yj*xk)
     B%jx3(i,j,k)=Jac*(yk*zi-zk*yi)
     B%jy3(i,j,k)=Jac*(zk*xi-xk*zi)
     B%jz3(i,j,k)=Jac*(xk*yi-yk*xi)
     B%kx3(i,j,k)=Jac*(yi*zj-zi*yj)
     B%ky3(i,j,k)=Jac*(zi*xj-xi*zj)
     B%kz3(i,j,k)=Jac*(xi*yj-yi*xj)
     enddo
	 enddo
	 enddo
!  (I,J,K)点的值, 标量方程的源项需要 (仅最密的网格使用)
    if(nMesh .eq. 1) then
    do k=1,B%nz-1
    do j=1,B%ny-1
    do i=1,B%nx-1
	 xi=B%xc(i+1,j,k)-B%xc(i-1,j,k) 
     yi=B%yc(i+1,j,k)-B%yc(i-1,j,k)
     zi=B%zc(i+1,j,k)-B%zc(i-1,j,k)
     xj=B%xc(i,j+1,k)-B%xc(i,j-1,k)   
     yj=B%yc(i,j+1,k)-B%yc(i,j-1,k)
     zj=B%zc(i,j+1,k)-B%zc(i,j-1,k)
     xk=B%xc(i,j,k+1)-B%xc(i,j,k-1) 
     yk=B%yc(i,j,k+1)-B%yc(i,j,k-1)
     zk=B%zc(i,j,k+1)-B%zc(i,j,k-1)
	  Jac1=xi*yj*zk+yi*zj*xk+zi*xj*yk-xi*zj*yk-yi*xj*zk-zi*yj*xk
      if(abs(Jac1) .lt. Lim_Zero) Jac1=Lim_Zero
	  Jac=1.d0/Jac1

     B%ix0(i,j,k)=Jac*(yj*zk-zj*yk)
     B%iy0(i,j,k)=Jac*(zj*xk-xj*zk)
     B%iz0(i,j,k)=Jac*(xj*yk-yj*xk)
     B%jx0(i,j,k)=Jac*(yk*zi-zk*yi)
     B%jy0(i,j,k)=Jac*(zk*xi-xk*zi)
     B%jz0(i,j,k)=Jac*(xk*yi-yk*xi)
     B%kx0(i,j,k)=Jac*(yi*zj-zi*yj)
     B%ky0(i,j,k)=Jac*(zi*xj-xi*zj)
     B%kz0(i,j,k)=Jac*(xi*yj-yi*xj)
    enddo
    enddo
	enddo
	endif
	
	enddo

 
  end subroutine Comput_Goemetric_var

!-------------------------------------------------
  subroutine check_mesh_quality
   use   Global_Var
   implicit none
   integer::nMesh
   do nMesh=1,Num_Mesh
   call check_mesh_quality_onemesh(nMesh)
   enddo
   end


! 检查网格质量
! 检查方法： 网格的连续性 （体积的连续性、法方向的连续性）
  subroutine check_mesh_quality_onemesh(nMesh)
   use   Global_Var
   use, intrinsic :: ieee_arithmetic
   implicit none
   integer :: i,j,k,m,nx,ny,nz,nMesh,i1,j1,k1,i2,j2,k2,im,jm,km,in,jn,kn
   real(PRE_EC):: fl,ft,flmax,ftmax,fl1,fl2,fl3,ft1,ft2,ft3,Af,Aa
   real(PRE_EC):: x1,x2,y1,y2,z1,z2,Vmax,Vmin
   character(len=50):: filename

   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B

   ! Initialise worst-cell indices: on a perfectly uniform Cartesian grid
   ! fl==1 and ft==0 everywhere, so the search loops below never update them
   ! and they would stay undefined (index 0 -> out-of-bounds write).
   i1=1; j1=1; k1=1
   i2=1; j2=1; k2=1

   if(my_id .eq. 0)    print*, "Check mesh quality ......"

   MP=>Mesh(nMesh)
   do m=1,MP%Num_Block   
     B => MP%Block(m)
     nx=B%nx; ny=B%ny ; nz=B%nz
 
 ! ----统计最大、最小网格 ----------------------
      Vmax=B%vol(1,1,1)
	  Vmin=B%vol(1,1,1)
	  im=1; jm=1; km=1
	  in=1; jn=1; kn=1
      
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	   if(B%vol(i,j,k) > Vmax) then
	    Vmax=B%vol(i,j,k)
		im=i ; jm=j ; km=k
	   endif
	   if(B%vol(i,j,k) < Vmin) then
	    Vmin=B%vol(i,j,k)
		in=i; jn=j; kn=k
	   endif
	   enddo
	   enddo
	   enddo



 !-----计算网格长度比、网格线的转角 ------------
 
 	 flmax=1.d0
	 ftmax=0.d0

!  长度比


     do k=1,nz-1
	 do j=1,ny-1
 	 do i=1,nx-1
      fl1=sqrt((B%xc(i,j,k)-B%xc(i-1,j,k))**2+(B%yc(i,j,k)-B%yc(i-1,j,k))**2+(B%zc(i,j,k)-B%zc(i-1,j,k))**2) &
	     /sqrt((B%xc(i,j,k)-B%xc(i+1,j,k))**2+(B%yc(i,j,k)-B%yc(i+1,j,k))**2+(B%zc(i,j,k)-B%zc(i+1,j,k))**2)
      fl2=sqrt((B%xc(i,j,k)-B%xc(i,j-1,k))**2+(B%yc(i,j,k)-B%yc(i,j-1,k))**2+(B%zc(i,j,k)-B%zc(i,j-1,k))**2) &
	     /sqrt((B%xc(i,j,k)-B%xc(i,j+1,k))**2+(B%yc(i,j,k)-B%yc(i,j+1,k))**2+(B%zc(i,j,k)-B%zc(i,j+1,k))**2)
      fl3=sqrt((B%xc(i,j,k)-B%xc(i,j,k-1))**2+(B%yc(i,j,k)-B%yc(i,j,k-1))**2+(B%zc(i,j,k)-B%zc(i,j,k-1))**2) &
	     /sqrt((B%xc(i,j,k)-B%xc(i,j,k+1))**2+(B%yc(i,j,k)-B%yc(i,j,k+1))**2+(B%zc(i,j,k)-B%zc(i,j,k+1))**2)
	        
         if(fl1 .lt. 1.d0) fl1=1.d0/fl1
         if(fl2 .lt. 1.d0) fl2=1.d0/fl2
         if(fl3 .lt. 1.d0) fl3=1.d0/fl3
         fl=max(fl1,fl2,fl3)
     
	    x1=B%xc(i,j,k)-B%xc(i-1,j,k) ; x2= B%xc(i+1,j,k)-B%xc(i,j,k)
		y1=B%yc(i,j,k)-B%yc(i-1,j,k) ; y2= B%yc(i+1,j,k)-B%yc(i,j,k)
		z1=B%zc(i,j,k)-B%zc(i-1,j,k) ; z2= B%zc(i+1,j,k)-B%zc(i,j,k)
        ft1=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
    
	    x1=B%xc(i,j,k)-B%xc(i,j-1,k) ; x2= B%xc(i,j+1,k)-B%xc(i,j,k)
		y1=B%yc(i,j,k)-B%yc(i,j-1,k) ; y2= B%yc(i,j+1,k)-B%yc(i,j,k)
		z1=B%zc(i,j,k)-B%zc(i,j-1,k) ; z2= B%zc(i,j+1,k)-B%zc(i,j,k)
        ft2=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))

	    x1=B%xc(i,j,k)-B%xc(i,j,k-1) ; x2= B%xc(i,j,k+1)-B%xc(i,j,k)
		y1=B%yc(i,j,k)-B%yc(i,j,k-1) ; y2= B%yc(i,j,k+1)-B%yc(i,j,k)
		z1=B%zc(i,j,k)-B%zc(i,j,k-1) ; z2= B%zc(i,j,k+1)-B%zc(i,j,k)
        ft3=(x1*x2+y1*y2+z1*z2)/sqrt((x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2))
        ft=min(1.d0,1.d0*min(ft1,ft2,ft3))
        ft=acos(ft)           ! 网格线折角 (容易出现NaN)
        
        Af=0.9d0*exp(-4.d0*(fl-1.d0)**2)+0.1d0
        Aa=0.9d0*exp(-(4.d0/3.1415926535d0*ft)**2)+0.1d0
        B%dtime_mesh(i,j,k)=min(Af,Aa)
       
	    if(ieee_is_nan(Af) .or. ieee_is_nan(Aa)) then
		 print*, "--------Find bad grid ------------"
		 print*, "Af, Aa=", Af,Aa
		 print*, "Block_no, i,j,k=",B%block_no, i,j,k
		 write(*,"(7E30.20)") fl,ft,ft1,ft2,ft3,min(ft1,ft2,ft3),acos(min(ft1,ft2,ft3))
		endif

!---------找出质量最差的网格，输出-------------------	   
       if(fl .gt. flmax) then
	    flmax=fl
	    i1=i
	    j1=j
	    k1=k
	   endif
   
       if(ft .gt. ftmax) then
	   ftmax=ft
	   i2=i
	   j2=j
	   k2=k
	   endif


	  enddo
      enddo
	  enddo

       if(nMesh .eq. 1) then
        open(103,file="mesh-quality.dat",position="append")
        write(103,*) "------------------Block ", m, "---------------------------"
        write(103,*) "Cell Number=", (nx-1)*(ny-1)*(nz-1)
        write(103,*) "Max volume=", Vmax, im,jm,km
		write(103,*) "Min Volume=", Vmin, in,jn,kn
	    write(103,*) "============="
	    write(103,*) "max grid factor=",flmax, i1,j1,k1
		write(103,*) "dtime_mesh=", B%dtime_mesh(i1,j1,k1)
	    write(103,*) "max grid-line turn angle (degree)=",ftmax*180.d0/3.1415926535d0, i2,j2,k2
	    write(103,*) "dtime_mesh=", B%dtime_mesh(i2,j2,k2)
        close(103)
	   endif
     
    enddo

    if(nMesh .eq. 1) then
       if(IF_Debug .eq. 1) then
	    write(filename, "('mesh-quality-',I5.5,'.dat')") my_id
		open(106, file=filename)
        write(106,*) "variables=x,y,z,dtfact"
	    do m=1,MP%Num_Block   
         B => MP%Block(m)
         nx=B%nx; ny=B%ny ; nz=B%nz
         write(106,*) "zone i=", nx-1, " j= ", ny-1, " k= ",nz-1
		do k=1,nz-1
		do j=1,ny-1
        do i=1,nx-1
		write(106,"(4f20.10)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%dtime_mesh(i,j,k)
	    enddo
		enddo
		enddo
	    enddo
	    close(106)
       endif

!	  close(103)
     if(my_id .eq. 0)    print*, "Check mesh quality OK"

    endif


	end  subroutine check_mesh_quality_onemesh


!------------------------------------------------------------------    
!  设定各重网格上的控制信息
  subroutine set_control_para
   use Global_var
   implicit none
   integer nMesh
   TYPE (Mesh_TYPE),pointer:: MP
   MP=>Mesh(1)            ! 最细的网格
!  最细网格上的控制参数与主控制参数相同
   MP%Iflag_turbulence_model=Iflag_turbulence_model
   MP%Iflag_Scheme=Iflag_Scheme
   MP%IFlag_flux=IFlag_flux
   MP%IFlag_Reconstruction=IFlag_Reconstruction
   MP%Bound_Scheme=Bound_scheme   !  边界格式

!  设定粗网格上的控制参数
   do nMesh=2,Num_Mesh
     MP=>Mesh(nMesh)
     MP%Iflag_turbulence_model=Turbulence_NONE    ! 粗网格不使用湍流模型
     MP%Iflag_Scheme=Scheme_UD1                   ! 粗网格使用1阶迎风格式
     MP%IFlag_flux=IFlag_flux                     ! 粗网格的通量分裂技术、时间推进近似及重构技术与细网格相同
     MP%IFlag_Reconstruction=IFlag_Reconstruction
     MP%Bound_Scheme=Scheme_UD1                   ! 粗网格边界点使用1阶格式
   enddo
  
  end subroutine set_control_para




  subroutine check_mesh_multigrid 
   use Global_var
   implicit none
   integer,allocatable,dimension(:):: NI,NJ,NK
   integer:: NB,NST,m,k,NN,Km,Km_grid,N_Cell,Ntmp,Bsub,Ksub
   integer:: ib,ie,jb,je,kb,ke,bc,ist,iend,jst,jend,kst,kend
   print*, "Check if Multi-Grid can be used ..."

   if( Mesh_File_Format .eq. 1) then   ! 格式文件
     open(99,file="Mesh3d.x")
     read(99,*) NB
     allocate(NI(NB),NJ(NB),NK(NB))
     read(99,*) (NI(m), NJ(m), NK(m), m=1,NB)
     close(99)
   else
     open(99,file="Mesh3d.x",form="unformatted")
     read(99) NB
     allocate(NI(NB),NJ(NB),NK(NB))
     read(99) (NI(m), NJ(m), NK(m), m=1,NB)
     close(99)
  endif


   N_Cell=0
   Km_grid=NI(1)  ! 初始值    
   do m=1,NB 
	 N_Cell=N_Cell+(NI(m)-1)*(NJ(m)-1)*(NK(m)-1)  ! 统计总网格单元数 
!  判断可使用的网格重数      
 	 Km=1
	 NN=2
!  判断准则： 网格数-1 能被2**km 整除， 且最稀的网格单元数不小于2
     do while( mod((NI(m)-1),NN) .eq. 0 .and. (NI(m)-1)/NN .ge. 2     &
		     .and. mod((NJ(m)-1),NN) .eq. 0 .and. (NJ(m)-1)/NN .ge. 2    &
		     .and. mod((NK(m)-1),NN) .eq. 0 .and. (NK(m)-1)/NN .ge. 2) 
       Km=Km+1              ! 所允许的网格重数
	   NN=NN*2
     enddo
     Km_grid=min(Km_grid,Km)
   enddo
!   Print*, " Finished check Mesh3d.x,  Most stage is ", Km_grid
!   print*,  "Check bc3d.inp ..." 
   open(88,file="bc3d.inp")
   read(88,*)
   read(88,*) 
   do m=1,NB
     read(88,*)
     read(88,*)
     read(88,*) Bsub    !number of the subface in the Block m
     do ksub=1, Bsub
       read(88,*)  ib,ie,jb,je,kb,ke,bc
	   if(bc .lt. 0) read(88,*)
	   ist=min(abs(ib),abs(ie)) ; iend=max(abs(ib),abs(ie))
       jst=min(abs(jb),abs(je));  jend=max(abs(jb),abs(je))
       kst=min(abs(kb),abs(ke));  kend=max(abs(kb),abs(ke))

	   NN=1
	   Km=1
       do while( mod((ist-1),NN) .eq. 0 .and. mod((iend-1),NN) .eq.0       &
		       .and. mod((jst-1),NN) .eq. 0 .and. mod((jend-1),NN) .eq. 0      &
		       .and. mod((kst-1),NN) .eq. 0 .and. mod((kend-1),NN) .eq. 0    ) 
         NN=NN*2
		 Km=Km+1
	   enddo
       Km_grid=min(Km_grid,Km)
     enddo
   enddo
   close(88)

!--------------------------------------------------------- 
   print*, "Total Block number is ", NB, "Total Cell number is " , N_Cell
!   print*, "Most stage number of multi-grid is ", Km_grid
!   print*, "-------------------------------------------------"
!-------------------------------------------------------------
   if(Num_Mesh .gt. Km_grid .or. Num_mesh .gt. 3) then
     print*, "Wrong !, Stage number of multi-grid error !!!"
	 stop 
   endif
   print*, "Check multigrid OK"

!-------------------------------------------------------
   deallocate(NI,NJ,NK)

  end subroutine check_mesh_multigrid


   subroutine Output_mesh_debug
   use Global_var
   implicit none
   integer:: m,nx,ny,nz,i,j,k
   Type (Block_TYPE),pointer:: B
   character(len=50):: filename
   write(filename,"('Meshdeb-',I5.5,'.dat')") my_id
   open(50,file=filename)
    
!  输出网格文件 tecplot格式；   
   write(50,*) "variables=x,y,z"
   do m=1,Mesh(1)%Num_Block
     B=>Mesh(1)%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
     write(50,*) "zone i=",nx+2, " j= ", ny+2, " k= ", nz+2
	 do k=0,nz+1
	 do j=0,ny+1
	 do i=0,nx+1
	 write(50,"(3f20.10)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k)
	 enddo
	 enddo
	 enddo
	enddo
	close(50)
   end
 

 
   subroutine debug1
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: nx,ny,nz,m,i,j,k
   real(PRE_EC):: u1,u2   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
    
	MP=>Mesh(1)
    u1=0.d0
	u2=0.d0
	do m=1,MP%Num_Block
    B=>MP%Block(m)
	nx=B%nx; ny=B%ny; nz=B%nz
	do k=1-LAP,nz+LAP-1
	do j=1-LAP,ny+LAP-1
	do i=1-LAP,nx+LAP-1
	 u1=u1+B%U(1,i,j,k)**2
	 u2=u2+B%U(5,i,j,k)**2
	enddo
	enddo
	enddo
	enddo
	print*, "************ u1,u2=",u1,u2
	end



   subroutine debug2(nx,ny,nz)
   use Global_Var
   use Flow_Var 
   implicit none
   integer:: nx,ny,nz,i,j,k
   real(PRE_EC):: f1,f2   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
    
	MP=>Mesh(1)
    f1=0.d0
	f2=0.d0
	do k=1,nz
	do j=1,ny
	do i=1,nx
	 f1=f1+Flux(1,i,j,k)**2
     f2=f2+Flux(5,i,j,k)**2
	enddo
	enddo
	enddo

	print*, "************ f1,f2=",f1,f2
		end

!===============================================================================
! register_bc_interfaces -- scan Mesh(1) for gridgen "generic: 8" coupling
! faces (bc == BC_INTERFACE) and register them into mod_interface's
! Interface_List.
!
! Phase 4 change: each *cell* of an interface subface is registered as its
! own entry (not one aggregate entry per subface), with the 4 quad vertex
! coordinates stored in verts(3,4).  This lets the phase-4 matcher project
! an unstructured face centroid onto the structured quad and compute
! bilinear interpolation weights.
!
! Pure marker/registration: BC_MSG_TYPE is left untouched, no coupling
! algorithm and no SI conversion here.  Called once from init() after
! Comput_Goemetric_var(1) so block coordinates x/y/z are available.
!
! Face orientation (Convert_bc convention): face=1/4 -> i-/i+ (constant i),
! face=2/5 -> j-/j+, face=3/6 -> k-/k+.  The face grid is the two varying
! directions over [jb..je-1] x [kb..ke-1] etc., using node coordinates.
!===============================================================================
   subroutine register_bc_interfaces
    use Global_var
    use mod_interface, only: Interface_FACE_TYPE, Interface_List, Num_Interface, &
                             BC_INTERFACE, PEER_STRUCT
    implicit none
    integer:: mBlock, ksub, cnt, nquad
    integer:: nx, ny, nz
    integer:: i0, j0, k0
    integer:: a0_, a1_, b0_, b1_, fa, a2, b2
    real(PRE_EC):: xa, ya, za, xb, yb, zb, xc_, yc_, zc_
    real(PRE_EC):: xd, yd, zd
    real(PRE_EC):: ux, uy, uz, vx, vy, vz
    real(PRE_EC):: cnx, cny, cnz, tri_area, area_tot, nlen
    real(PRE_EC):: cx, cy, cz
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc

    MP => Mesh(1)

    ! ---- first pass: count total interface quad cells ----
    nquad = 0
    do mBlock = 1, MP%Num_Block
       B => MP%Block(mBlock)
       do ksub = 1, B%subface
          Bc => B%bc_msg(ksub)
          if (Bc%bc /= BC_INTERFACE) cycle
          select case (Bc%face)
          case (1, 4)
             nquad = nquad + (Bc%je - Bc%jb)*(Bc%ke - Bc%kb)
          case (2, 5)
             nquad = nquad + (Bc%ie - Bc%ib)*(Bc%ke - Bc%kb)
          case default
             nquad = nquad + (Bc%ie - Bc%ib)*(Bc%je - Bc%jb)
          end select
       enddo
    enddo
    Num_Interface = nquad
    if (allocated(Interface_List)) deallocate(Interface_List)
    allocate(Interface_List(Num_Interface))

    ! stay silent when no coupling face exists so existing cases keep
    ! byte-identical logs; only announce when interfaces are present.
    if (Num_Interface == 0) return

    ! ---- second pass: one entry per interface quad cell ----
    cnt = 0
    do mBlock = 1, MP%Num_Block
       B => MP%Block(mBlock)
       nx = B%nx; ny = B%ny; nz = B%nz
       do ksub = 1, B%subface
          Bc => B%bc_msg(ksub)
          if (Bc%bc /= BC_INTERFACE) cycle

          fa = Bc%face
          select case (fa)
          case (1, 4)                      ! i- / i+ : constant i, vary (j,k)
             i0 = Bc%ib
             a0_ = Bc%jb; a1_ = Bc%je      ! j range
             b0_ = Bc%kb; b1_ = Bc%ke      ! k range
          case (2, 5)                      ! j- / j+ : constant j, vary (i,k)
             j0 = Bc%jb
             a0_ = Bc%ib; a1_ = Bc%ie      ! i range
             b0_ = Bc%kb; b1_ = Bc%ke      ! k range
          case default                     ! k- / k+ : constant k, vary (i,j)
             k0 = Bc%kb
             a0_ = Bc%ib; a1_ = Bc%ie      ! i range
             b0_ = Bc%jb; b1_ = Bc%je      ! j range
          end select

          do a2 = a0_, a1_ - 1
          do b2 = b0_, b1_ - 1
             cnt = cnt + 1

             ! four corners of this quad in (i,j,k)
             call iface_node(fa, i0, j0, k0, a2,   b2,   B, xa, ya, za)
             call iface_node(fa, i0, j0, k0, a2+1, b2,   B, xb, yb, zb)
             call iface_node(fa, i0, j0, k0, a2+1, b2+1, B, xc_, yc_, zc_)
             call iface_node(fa, i0, j0, k0, a2,   b2+1, B, xd, yd, zd)

             ! identity
             Interface_List(cnt)%solver   = PEER_STRUCT
             Interface_List(cnt)%block_no = B%block_no
             Interface_List(cnt)%face     = Bc%face
             Interface_List(cnt)%f_no     = Bc%f_no
             Interface_List(cnt)%ib = Bc%ib; Interface_List(cnt)%ie = Bc%ie
             Interface_List(cnt)%jb = Bc%jb; Interface_List(cnt)%je = Bc%je
             Interface_List(cnt)%kb = Bc%kb; Interface_List(cnt)%ke = Bc%ke
             Interface_List(cnt)%match_state = 0
             Interface_List(cnt)%peer_id     = 0

             ! ghost-cell indices for BC application (phase 6)
             Interface_List(cnt)%ic = a2; Interface_List(cnt)%jc = b2; Interface_List(cnt)%kc = 0
             Interface_List(cnt)%ig = a2; Interface_List(cnt)%jg = b2; Interface_List(cnt)%kg = 0
             select case (fa)
             case (1)  ! i- face: inner cell at ib, ghost at ib-1
                Interface_List(cnt)%ic = i0; Interface_List(cnt)%jc = a2; Interface_List(cnt)%kc = b2
                Interface_List(cnt)%ig = i0-1; Interface_List(cnt)%jg = a2; Interface_List(cnt)%kg = b2
             case (4)  ! i+ face: inner cell at ie-1, ghost at ie
                Interface_List(cnt)%ic = i0-1; Interface_List(cnt)%jc = a2; Interface_List(cnt)%kc = b2
                Interface_List(cnt)%ig = i0; Interface_List(cnt)%jg = a2; Interface_List(cnt)%kg = b2
             case (2)  ! j- face: inner cell at jb, ghost at jb-1
                Interface_List(cnt)%ic = a2; Interface_List(cnt)%jc = j0; Interface_List(cnt)%kc = b2
                Interface_List(cnt)%ig = a2; Interface_List(cnt)%jg = j0-1; Interface_List(cnt)%kg = b2
             case (5)  ! j+ face: inner cell at je-1, ghost at je
                Interface_List(cnt)%ic = a2; Interface_List(cnt)%jc = j0-1; Interface_List(cnt)%kc = b2
                Interface_List(cnt)%ig = a2; Interface_List(cnt)%jg = j0; Interface_List(cnt)%kg = b2
             case (3)  ! k- face: inner cell at kb, ghost at kb-1
                Interface_List(cnt)%ic = a2; Interface_List(cnt)%jc = b2; Interface_List(cnt)%kc = k0
                Interface_List(cnt)%ig = a2; Interface_List(cnt)%jg = b2; Interface_List(cnt)%kg = k0-1
             case (6)  ! k+ face: inner cell at ke-1, ghost at ke
                Interface_List(cnt)%ic = a2; Interface_List(cnt)%jc = b2; Interface_List(cnt)%kc = k0-1
                Interface_List(cnt)%ig = a2; Interface_List(cnt)%jg = b2; Interface_List(cnt)%kg = k0
             end select

             ! 4 vertices, CCW around the face normal
             Interface_List(cnt)%nv = 4
             allocate( Interface_List(cnt)%verts(3,4) )
             Interface_List(cnt)%verts(:,1) = (/ xa,  ya,  za  /)
             Interface_List(cnt)%verts(:,2) = (/ xb,  yb,  zb  /)
             Interface_List(cnt)%verts(:,3) = (/ xc_, yc_, zc_ /)
             Interface_List(cnt)%verts(:,4) = (/ xd,  yd,  zd  /)

             ! bounding box over the four corners
             Interface_List(cnt)%bbox_min = min((/ xa, ya, za /), &
                                                 (/ xb, yb, zb /), &
                                                 (/ xc_, yc_, zc_ /), &
                                                 (/ xd, yd, zd /))
             Interface_List(cnt)%bbox_max = max((/ xa, ya, za /), &
                                                 (/ xb, yb, zb /), &
                                                 (/ xc_, yc_, zc_ /), &
                                                 (/ xd, yd, zd /))

             ! area, centroid, normal from the two triangles (A,B,C),(A,C,D)
             area_tot = 0.d0
             cx = 0.d0; cy = 0.d0; cz = 0.d0
             cnx = 0.d0; cny = 0.d0; cnz = 0.d0
             ux = xb-xa; uy = yb-ya; uz = zb-za
             vx = xc_-xa; vy = yc_-ya; vz = zc_-za
             call iface_tri(ux,uy,uz,vx,vy,vz, xa,ya,za, xb,yb,zb, xc_,yc_,zc_, &
                            tri_area, cnx,cny,cnz, cx,cy,cz)
             area_tot = area_tot + tri_area
             ux = xc_-xa; uy = yc_-ya; uz = zc_-za
             vx = xd-xa;  vy = yd-ya;  vz = zd-za
             call iface_tri(ux,uy,uz,vx,vy,vz, xa,ya,za, xc_,yc_,zc_, xd,yd,zd, &
                            tri_area, cnx,cny,cnz, cx,cy,cz)
             area_tot = area_tot + tri_area

             Interface_List(cnt)%area = area_tot
             if (area_tot > 0.d0) then
                Interface_List(cnt)%centroid = (/ cx/area_tot, cy/area_tot, cz/area_tot /)
             else
                Interface_List(cnt)%centroid = (/ 0.d0, 0.d0, 0.d0 /)
             endif
             nlen = sqrt(cnx*cnx + cny*cny + cnz*cnz)
             if (nlen > 0.d0) then
                Interface_List(cnt)%normal = (/ cnx/nlen, cny/nlen, cnz/nlen /)
             else
                Interface_List(cnt)%normal = (/ 0.d0, 0.d0, 0.d0 /)
             endif
          enddo
          enddo
       enddo
    enddo

    if (my_id == 0) then
       print*, "Registered", Num_Interface, &
                               " gridgen generic:8 coupling interface cell-face(s)."
       ! report aggregate stats per (block, face) pair
       call report_struct_ifaces
    end if

   contains

      ! map the (constant-dir, two varying dirs) triplet back to (i,j,k) and
      ! return the node coordinates from block B.
      subroutine iface_node(fa, ci, cj, ck, av, bv, B, xo, yo, zo)
       implicit none
       integer, intent(in):: fa, ci, cj, ck, av, bv
       Type (Block_TYPE), intent(in):: B
       real(PRE_EC), intent(out):: xo, yo, zo
       integer:: ii, jj, kk
       select case (fa)
       case (1, 4);  ii = ci; jj = av; kk = bv    ! const i, vary (j,k)
       case (2, 5);  ii = av; jj = cj; kk = bv    ! const j, vary (i,k)
       case default; ii = av; jj = bv; kk = ck    ! const k, vary (i,j)
       end select
       xo = B%x(ii, jj, kk); yo = B%y(ii, jj, kk); zo = B%z(ii, jj, kk)
      end subroutine iface_node

      ! accumulate one triangle: area-weighted centroid and un-normalised
      ! normal (cross product magnitude == 2*area; consistent weighting).
      subroutine iface_tri(ux,uy,uz, vx,vy,vz, x1,y1,z1, x2,y2,z2, x3,y3,z3, &
                           tri_area, cnx,cny,cnz, cx,cy,cz)
       implicit none
       real(PRE_EC), intent(in):: ux,uy,uz, vx,vy,vz
       real(PRE_EC), intent(in):: x1,y1,z1, x2,y2,z2, x3,y3,z3
       real(PRE_EC), intent(out):: tri_area
       real(PRE_EC), intent(inout):: cnx,cny,cnz, cx,cy,cz
       real(PRE_EC):: px, py, pz
       px = uy*vz - uz*vy
       py = uz*vx - ux*vz
       pz = ux*vy - uy*vx
       tri_area = 0.5d0*sqrt(px*px + py*py + pz*pz)
       cnx = cnx + px; cny = cny + py; cnz = cnz + pz
       cx = cx + tri_area*(x1+x2+x3)/3.d0
       cy = cy + tri_area*(y1+y2+y3)/3.d0
       cz = cz + tri_area*(z1+z2+z3)/3.d0
      end subroutine iface_tri

      ! aggregate geometry report per (block, face) pair (rank 0 only)
      subroutine report_struct_ifaces
       implicit none
       integer:: i, j, blk, fc, nf
       real(PRE_EC):: garea, gc(3), gn(3), gbmin(3), gbmax(3), nrm
       logical, allocatable:: seen(:)
       allocate(seen(Num_Interface), source = .false.)
       do i = 1, Num_Interface
          if (seen(i)) cycle
          blk = Interface_List(i)%block_no
          fc  = Interface_List(i)%face
          garea = 0.d0; gc = 0.d0; gn = 0.d0
          gbmin =  huge(1.d0); gbmax = -huge(1.d0)
          nf = 0
          do j = 1, Num_Interface
             if (Interface_List(j)%solver /= PEER_STRUCT) cycle
             if (Interface_List(j)%block_no /= blk .or. &
                 Interface_List(j)%face /= fc) cycle
             seen(j) = .true.
             nf = nf + 1
             garea = garea + Interface_List(j)%area
             gc = gc + Interface_List(j)%area * Interface_List(j)%centroid
             gn = gn + Interface_List(j)%normal * Interface_List(j)%area
             gbmin = min(gbmin, Interface_List(j)%bbox_min)
             gbmax = max(gbmax, Interface_List(j)%bbox_max)
          enddo
          write(*,'(a,i0,a,i0,a,i0)') '  [struct] block=', blk, &
             '  face=', fc, '  nfaces=', nf
          nrm = sqrt(gn(1)**2 + gn(2)**2 + gn(3)**2)
          if (garea > 0.d0) then
             write(*,'(a,3(1x,es12.4))') '           centroid : ', gc/garea
             if (nrm > 0.d0) write(*,'(a,3(1x,es12.4))') '           normal   : ', gn/nrm
          endif
          write(*,'(a,es12.4)')       '           area     : ', garea
          write(*,'(a,3(1x,es12.4))') '           bbox min : ', gbmin
          write(*,'(a,3(1x,es12.4))') '           bbox max : ', gbmax
       enddo
       deallocate(seen)
      end subroutine report_struct_ifaces

   end subroutine register_bc_interfaces

  end module mod_struct_grid
