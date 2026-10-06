!===============================================================================
! mod_struct_io.f90 -- structured solver: mesh/flow I/O, wall distance,
! forces, smoothing, time-averaged output.
! Encapsulates sub_IO + sub_comput_dw + sub_Post + sub_Post_timeAverage
! (phase 2b, batch B4).  The legacy Wall_dist helper module is kept verbatim
! above mod_struct_io.
!===============================================================================

  module Wall_dist
   use precision_EC
   implicit none
   integer,save:: Npw_total                ! 总壁面网格数
   real(PRE_EC),allocatable,dimension(:,:):: Xw   ! 壁面点的坐标
  end module Wall_dist

module mod_struct_io
   ! Running sample counter for time-averaging; init_average resets it and
   ! Time_average increments it.  Relocated from Global_Var (phase 2b s3):
   ! every access is inside this module.
   integer, save :: Istep_average = 0
contains


!  读取网格信息
!  根进程(0号进程)读取数据，其他进程接收
    subroutine read_main_mesh
    use Global_var
    implicit none
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    real(PRE_EC),allocatable,dimension(:,:,:,:):: Ux
    integer:: NB,m,nx,ny,nz,i,j,k,mt,Num_data
	integer:: Send_to_ID,tag,ierr, status(MPI_status_size)
    integer,allocatable,dimension(:):: NI,NJ,NK

     MP=>Mesh(1)
 
  if(my_id .eq. 0) then 
    print*, " read main mesh ..."
	if( Mesh_File_Format .eq. 1) then   ! 格式文件
     open(99,file="Mesh3d.x")
     read(99,*) NB   ! Block number
     allocate( NI(NB),NJ(NB),NK(NB) )
     read(99,*) (NI(k), NJ(k), NK(k), k=1,NB)
	else                                ! 无格式文件
     open(99,file="Mesh3d.x",form="unformatted")
     read(99) NB                       ! 总块数
     allocate( NI(NB),NJ(NB),NK(NB) )
     read(99) (NI(k), NJ(k), NK(k), k=1,NB)
    endif

!----------------------------------------
   do m=1,NB
!     print*, "block=",m
	 nx=NI(m); ny=NJ(m); nz=NK(m)
	 allocate(Ux(nx,ny,nz,3))
	   
     if( Mesh_File_Format .eq. 1) then
       read(99,*) (((Ux(i,j,k,1),i=1,nx),j=1,ny),k=1,nz) , &
                  (((Ux(i,j,k,2),i=1,nx),j=1,ny),k=1,nz) , &
                  (((Ux(i,j,k,3),i=1,nx),j=1,ny),k=1,nz)
	 else
       read(99)   (((Ux(i,j,k,1),i=1,nx),j=1,ny),k=1,nz) , &
                  (((Ux(i,j,k,2),i=1,nx),j=1,ny),k=1,nz) , &
                  (((Ux(i,j,k,3),i=1,nx),j=1,ny),k=1,nz)
	 endif
   if(B_proc(m) .eq. 0) then            ! 这些块属于根进程
      mt=B_n(m)                          ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	  do k=1,nz
	  do j=1,ny
	  do i=1,nx
	   B%x(i,j,k)=Ux(i,j,k,1)
	   B%y(i,j,k)=Ux(i,j,k,2)
	   B%z(i,j,k)=Ux(i,j,k,3)
	  enddo
	  enddo
 	  enddo
     else                        ! 将该块数据发送出
	   Num_data=nx*ny*nz*3
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
!	   call MPI_Bsend(Ux,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
	   call MPI_send(Ux,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
     endif
     deallocate(Ux)
   enddo

   deallocate(NI,NJ,NK)
   close(99)
  
  else     ! 非根节点
    do m=1,MP%Num_Block     ! 本进程包含的块
     B=>MP%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
   	 allocate(Ux(nx,ny,nz,3))
   	 Num_data=nx*ny*nz*3
 	 tag=m
	 call MPI_Recv(Ux,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,Status,ierr )
      do k=1,nz
	  do j=1,ny
	  do i=1,nx
	   B%x(i,j,k)=Ux(i,j,k,1)
	   B%y(i,j,k)=Ux(i,j,k,2)
	   B%z(i,j,k)=Ux(i,j,k,3)
	  enddo
	  enddo
 	  enddo
      deallocate(Ux)
     enddo
   endif
  
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "read Mesh3d.x OK"
 end subroutine read_main_mesh

!-------------------------------------------------------------------------------------
!  读取几何量： 到壁面的距离
!  根进程(0号进程)读取数据，其他进程接收
    subroutine read_dw
    use Global_var
    implicit none
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    real(PRE_EC),allocatable,dimension(:,:,:):: dw
    integer:: NB,m,nx,ny,nz,i,j,k,mt,Num_data
	integer:: Send_to_ID,tag,ierr, status(MPI_status_size)
	logical:: Ext

     MP=>Mesh(1)

!------根进程读取数据----------------------------- 
  if(my_id .eq. 0) then 
    print*, " read distance to the wall:  wall_dist.dat"
  
   open(99,file="wall_dist.dat",form="unformatted")

   do m=1,Total_block
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(dw(nx-1,ny-1,nz-1))
     read(99) (((dw(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1) 
     
	 if(B_proc(m) .eq. 0) then            ! 这些块属于根进程
      mt=B_n(m)                          ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	   B%dw(i,j,k)=dw(i,j,k)
	  enddo
	  enddo
 	  enddo
     else                        ! 将该块数据发送出
	   Num_data=(nx-1)*(ny-1)*(nz-1)
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
	   call MPI_send(dw,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
     endif
     deallocate(dw)
   enddo
   close(99)
  
  else     ! 非根节点
    do m=1,MP%Num_Block     ! 本进程包含的块
     B=>MP%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
   	 allocate(dw(nx-1,ny-1,nz-1))
	 Num_data=(nx-1)*(ny-1)*(nz-1)
 	 tag=m
	 call MPI_Recv(dw,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,Status,ierr )
      do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	   B%dw(i,j,k)=dw(i,j,k)
	  enddo
	  enddo
 	  enddo
      deallocate(dw)
     enddo
   endif
  
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "read wall_dist.dat OK"
 end subroutine read_dw




!----------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!  读取流场: d,u,v,w,T;  SA, SST 中的标量
!  根进程(0号进程)读取数据，其他进程接收
    subroutine read_flow_data
    use Global_var
    implicit none
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    real(PRE_EC),allocatable,dimension(:,:,:,:):: U
    integer:: NB,NVAR1,m,m1,nx,ny,nz,i,j,k,mt,Num_data
	integer:: Send_to_ID,tag,ierr, status(MPI_status_size)
	logical:: Ex
     real(PRE_EC):: d1,u1,v1,w1,T1

     MP=>Mesh(1)
     NVAR1=MP%NVAR

!------根进程读取数据，发送到其他进程----------------------------- 
  if(my_id .eq. 0) then 
    print*, " read flow data:  flow3d.dat"
 
     open(99,file="flow3d.dat",form="unformatted")
    
 	 if(NVAR1 .eq. 6) then        ! 6个自变量
      Inquire(file="SA3d.dat",exist=Ex)
      if(Ex) then
        open(100,file="SA3d.dat",form="unformatted")
      endif
	 endif
     
	 if(NVAR1 .eq. 7) then          ! 7个自变量
     Inquire(file="SST3d.dat",exist=Ex)
     if(Ex) then
       open(101,file="SST3d.dat",form="unformatted")
     endif
     endif

 
   do m=1,Total_block
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
 	  allocate(U(0:nx,0:ny,0:nz,NVAR1))
      read(99)   ((((U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=1,5) 
     if(NVAR1 .eq. 6) then
       if(Ex) then
           read(100)   (((U(i,j,k,6),i=0,nx),j=0,ny),k=0,nz)
	   else
           do k=0,nz
		   do j=0,ny
		   do i=0,nx
		     U(i,j,k,6)=1.d0/Re
		   enddo
		   enddo
		   enddo
	   endif
     endif

     if(NVAR1 .eq. 7) then
      if(Ex) then
         read(101)   ((((U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=6,7)
	  else
	       do k=0,nz
		   do j=0,ny
		   do i=0,nx
		     U(i,j,k,6)=Kt_Inf
             U(i,j,k,7)=Wt_Inf
		   enddo
		   enddo
		   enddo
  	  endif
     endif
!------------发送--------------------------------
     
	 if(B_proc(m) .eq. 0) then            ! 这些块属于根进程
      mt=B_n(m)                           ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	  do k=0,nz
	  do j=0,ny
	  do i=0,nx
	  do m1=1,NVAR1
	    B%U(m1,i,j,k)=U(i,j,k,m1)
	  enddo
	  enddo
	  enddo
 	  enddo
     else                        ! 将该块数据发送出
	   Num_data=(nx+1)*(ny+1)*(nz+1)*NVAR1
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
	   call MPI_send(U,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
     endif
     deallocate(U)
   enddo
   close(99)

   if(EX) then
    if(NVAR1 .eq. 6) then
      close (100)
    else 
	  close(101)
    endif
   endif

  else     ! 非根节点
   
    do m=1,MP%Num_Block     ! 本进程包含的块
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz,NVAR1))
	  Num_data=(nx+1)*(ny+1)*(nz+1)*NVAR1
 	  tag=m
	 
	 call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,Status,ierr )
      
	   do k=0,nz
	   do j=0,ny
	   do i=0,nx
	   do m1=1,NVAR1
	     B%U(m1,i,j,k)=U(i,j,k,m1)
	   enddo
	   enddo
 	   enddo
	   enddo
       deallocate(U)
     enddo
   
   endif
   
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "read flow3d.dat OK"

!----------------------------------Transform data----------------
! 读入的数据位d,u,v,w,T, 转化为守恒变量
 do m=1,MP%Num_Block     ! 本进程包含的块
    B=>MP%Block(m)
    nx=B%nx; ny=B%ny; nz=B%nz
 
   do k=0,nz
   do j=0,ny
   do i=0,nx
        d1=B%U(1,i,j,k)
        u1=B%U(2,i,j,k)
        v1=B%U(3,i,j,k)
        w1=B%U(4,i,j,k)
        T1=B%U(5,i,j,k)

	    B%U(1,i,j,k)=d1
        B%U(2,i,j,k)=d1*u1
        B%U(3,i,j,k)=d1*v1
        B%U(4,i,j,k)=d1*w1
        B%U(5,i,j,k)=Cv*d1*T1+0.5d0*d1*(u1**2+v1**2+w1**2)
   enddo
   enddo
   enddo
 
 enddo

!---------------读入时间步-----------------------
  if(my_id .eq. 0) then
   Inquire(file="Step_mess.dat",exist=Ex)
    if(Ex) then
     open(88,file="Step_mess.dat")
     read(88,*) Mesh(1)%Kstep, Mesh(1)%tt
    endif
  print*, "Init data OK, Kstep,tt=", Mesh(1)%Kstep, Mesh(1)%tt
 endif

 call MPI_bcast(Mesh(1)%Kstep,1,MPI_INTEGER,0,  Struct_Comm,ierr)
 call MPI_bcast(Mesh(1)%tt,1,OCFD_DATA_TYPE,0,  Struct_Comm,ierr)

!---------------------------------------------------
 end subroutine read_flow_data








!----------------------------------------------------------------------
!  输出几何及物理量 （Plot3d格式）, 最细网格flow3d.dat  

  subroutine output_flow
   use Global_Var
   implicit none
   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   TYPE(BC_MSG_TYPE),pointer::Bc
   real(PRE_EC),allocatable,dimension(:,:,:,:):: U
   integer:: NB,NVAR1,m,m1,nx,ny,nz,i,j,k,mt,Num_data
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)

   character(len=50):: filename   
   
   MP=>Mesh(1)
   NVAR1=MP%NVAR


!   输出计算步数及时间信息
   if(my_id .eq. 0) then
     open(88,file="Step_mess.dat")
     write(88,*) MP%Kstep, MP%tt
     write(88,*) "---time Step, time ----"
     close(88)
   endif

!---------------------------------------------------------------
! print*, "write 3D data file ......"


 if(my_id .eq. 0) then
  
   print*, "write flow3d.dat ......"
   
   if(Iflag_savefile==0 ) then
     open(99,file="flow3d.dat",form="unformatted")    ! d,u,v,w,T
   else 
     write(filename, "('flow3d-',I8.8,'.dat')") Mesh(1)%Kstep
     open(99,file=filename,form="unformatted")    ! d,u,v,w,T
   endif




   if(NVAR1 .eq. 6) open(100,file="SA3d.dat",form="unformatted")    ! U6
   if(NVAR1 .eq. 7) open(101,file="SST3d.dat",form="unformatted")   ! U6,U7
  
   do m=1, Total_block   ! 全部块
     
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz,NVAR1))

	if(B_proc(m) .eq. 0) then            ! 这些块属于根进程
      mt=B_n(m)                           ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	 
	  do k=0,nz
	  do j=0,ny
	  do i=0,nx
		 U(i,j,k,1)=B%U(1,i,j,k)                     ! d
         U(i,j,k,2)=B%U(2,i,j,k)/B%U(1,i,j,k)        ! u
         U(i,j,k,3)=B%U(3,i,j,k)/B%U(1,i,j,k)        ! v
         U(i,j,k,4)=B%U(4,i,j,k)/B%U(1,i,j,k)        ! w
         U(i,j,k,5)=(B%U(5,i,j,k)-0.5d0*U(i,j,k,1)*(U(i,j,k,2)**2+U(i,j,k,3)**2+U(i,j,k,4)**2) )/(Cv*U(i,j,k,1))    ! T
        if(NVAR1 .eq. 6) then
		 U(i,j,k,6)=B%U(6,i,j,k)
		endif
		if(NVAR1 .eq. 7) then
		 U(i,j,k,6)=B%U(6,i,j,k)
		 U(i,j,k,7)=B%U(7,i,j,k)
		endif 
	  enddo
	  enddo
 	  enddo

    else                        ! 接收该块信息
	   Num_data=NVAR1*(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! 在该块中的编号
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, Struct_Comm,Status,ierr )
    endif
! write Data ....
    
	write(99) (((( U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=1,5)    
    
	if(NVAR1 .eq. 6) then
	  write(100) ((( U(i,j,k,6),i=0,nx),j=0,ny),k=0,nz)    
	endif
	if(NVAR1 .eq. 7) then
	  write(101) (((( U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=6,7)    
	endif 

	deallocate(U)

   enddo
   close(99)
   close(100)
   close(101)
 
 else     ! 非0节点

    do m=1,MP%Num_Block     ! 本进程包含的块
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz,NVAR1))
	  Num_data=(nx+1)*(ny+1)*(nz+1)*NVAR1
 	  tag=m
	 
	   do k=0,nz
	   do j=0,ny
	   do i=0,nx
		 U(i,j,k,1)=B%U(1,i,j,k)
         U(i,j,k,2)=B%U(2,i,j,k)/B%U(1,i,j,k)
         U(i,j,k,3)=B%U(3,i,j,k)/B%U(1,i,j,k)
         U(i,j,k,4)=B%U(4,i,j,k)/B%U(1,i,j,k)
         U(i,j,k,5)=(B%U(5,i,j,k)-0.5d0*U(i,j,k,1)*(U(i,j,k,2)**2+U(i,j,k,3)**2+U(i,j,k,4)**2) )/(Cv*U(i,j,k,1))
        if(NVAR1 .eq. 6) then
		 U(i,j,k,6)=B%U(6,i,j,k)
		endif
		if(NVAR1 .eq. 7) then
		 U(i,j,k,6)=B%U(6,i,j,k)
		 U(i,j,k,7)=B%U(7,i,j,k)
		endif 
	   enddo
	   enddo
 	   enddo
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "write flow3d.dat OK"

  end subroutine output_flow

!----------------------------------------------------------------------
!  Node-interpolated Plot3D function output (phase 10).
!  Cell-centred primitive variables (d,u,v,w,T) are averaged onto the
!  Mesh3d.x node distribution: interior nodes average the 8 adjacent
!  cells, boundary nodes fewer (clamped cell ranges).  File format
!  follows Mesh3d.x (Mesh_File_Format: 1 = ascii, otherwise unformatted).
!  Plot3D function layout: NB, then (NI,NJ,NK) per block, then
!  ((((Q(n,i,j,k),i=1,nx),j=1,ny),k=1,nz),n=1,5) per block.

  subroutine output_flow_nodes
   use Global_Var
   implicit none

   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:,:,:):: Uc   ! cell primitive
   real(PRE_EC),allocatable,dimension(:,:,:,:):: Qn   ! node primitive
   integer:: NB,m,nx,ny,nz,i,j,k,mt,Num_data,nv
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)
   integer:: i0,i1,j0,j1,k0,k1,ii,jj,kk,cnt

   MP=>Mesh(1)

 if(my_id .eq. 0) then
   print*, "write flow3d_node.dat (node-interpolated Plot3D) ......"
   if( Mesh_File_Format .eq. 1 ) then
     open(99,file="flow3d_node.dat",form="formatted")
     write(99,*) Total_block
     write(99,*) (bNi(m), bNj(m), bNk(m), m=1,Total_block)
   else
     open(99,file="flow3d_node.dat",form="unformatted")
     write(99) Total_block
     write(99) (bNi(m), bNj(m), bNk(m), m=1,Total_block)
   endif

   do m=1, Total_block   ! all blocks
     nx=bNi(m); ny=bNj(m); nz=bNk(m)
     allocate(Uc(0:nx,0:ny,0:nz,5))

    if(B_proc(m) .eq. 0) then            ! blocks on the root rank
      mt=B_n(m)                          ! local block index
      B=>MP%Block(mt)
      do k=0,nz
      do j=0,ny
      do i=0,nx
        Uc(i,j,k,1)=B%U(1,i,j,k)                     ! d
        Uc(i,j,k,2)=B%U(2,i,j,k)/B%U(1,i,j,k)        ! u
        Uc(i,j,k,3)=B%U(3,i,j,k)/B%U(1,i,j,k)        ! v
        Uc(i,j,k,4)=B%U(4,i,j,k)/B%U(1,i,j,k)        ! w
        Uc(i,j,k,5)=(B%U(5,i,j,k)-0.5d0*Uc(i,j,k,1)*(Uc(i,j,k,2)**2+Uc(i,j,k,3)**2+Uc(i,j,k,4)**2))/(Cv*Uc(i,j,k,1))  ! T
      enddo
      enddo
      enddo
    else                                 ! receive from the owner rank
       Num_data=5*(nx+1)*(ny+1)*(nz+1)
       Recv_from_ID=B_proc(m)
       tag=B_n(m)
      call MPI_Recv(Uc,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, Struct_Comm,Status,ierr )
    endif

    ! ---- cell -> node interpolation (average of adjacent cells) ----------
    ! Interior cells are 1..nx-1 (0 and nx are ghost).  Node i is bracketed
    ! by cells i-1 and i, clamped to the interior range.
    allocate(Qn(nx,ny,nz,5))
    do k=1,nz
      k0=max(k-1,1); k1=min(k,nz-1)
      do j=1,ny
        j0=max(j-1,1); j1=min(j,ny-1)
        do i=1,nx
          i0=max(i-1,1); i1=min(i,nx-1)
          do nv=1,5
            Qn(i,j,k,nv)=0.0d0
            cnt=0
            do kk=k0,k1
            do jj=j0,j1
            do ii=i0,i1
              Qn(i,j,k,nv)=Qn(i,j,k,nv)+Uc(ii,jj,kk,nv)
              cnt=cnt+1
            enddo
            enddo
            enddo
            Qn(i,j,k,nv)=Qn(i,j,k,nv)/real(cnt,PRE_EC)
          enddo
        enddo
      enddo
    enddo

    if( Mesh_File_Format .eq. 1 ) then
      write(99,*) (((( Qn(i,j,k,nv),i=1,nx),j=1,ny),k=1,nz),nv=1,5)
    else
      write(99) (((( Qn(i,j,k,nv),i=1,nx),j=1,ny),k=1,nz),nv=1,5)
    endif

    deallocate(Uc)
    deallocate(Qn)
   enddo
   close(99)

 else     ! non-root ranks send their cell primitive state

    do m=1,MP%Num_Block     ! blocks owned by this rank
      B=>MP%Block(m)
      nx=B%nx; ny=B%ny; nz=B%nz
      allocate(Uc(0:nx,0:ny,0:nz,5))
       do k=0,nz
       do j=0,ny
       do i=0,nx
        Uc(i,j,k,1)=B%U(1,i,j,k)
        Uc(i,j,k,2)=B%U(2,i,j,k)/B%U(1,i,j,k)
        Uc(i,j,k,3)=B%U(3,i,j,k)/B%U(1,i,j,k)
        Uc(i,j,k,4)=B%U(4,i,j,k)/B%U(1,i,j,k)
        Uc(i,j,k,5)=(B%U(5,i,j,k)-0.5d0*Uc(i,j,k,1)*(Uc(i,j,k,2)**2+Uc(i,j,k,3)**2+Uc(i,j,k,4)**2))/(Cv*Uc(i,j,k,1))
       enddo
       enddo
       enddo
      Num_data=5*(nx+1)*(ny+1)*(nz+1)
      tag=m
      call MPI_Send(Uc,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,ierr )
      deallocate(Uc)
    enddo

  endif

   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "write flow3d_node.dat OK"

  end subroutine output_flow_nodes




!----------------------------------------------------------------------
!  输出湍流粘性系数vt （Plot3d格式）, 最细网格vt.dat  

  subroutine output_vt
   use Global_Var
   implicit none
   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:,:):: U
   integer:: NB,m,m1,nx,ny,nz,i,j,k,mt,Num_data
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)

   character(len=50):: filename   
   
   MP=>Mesh(1)
!---------------------------------------------------------------
 if(my_id .eq. 0) then
   print*, "write vt.dat ......"
   open(99,file="vt.dat",form="unformatted")                    ! d,u,v,w,T
 
   do m=1, Total_block   ! 全部块
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz))

	if(B_proc(m) .eq. 0) then            ! 这些块属于根进程
      mt=B_n(m)                           ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	 
	  do k=0,nz
	  do j=0,ny
	  do i=0,nx
		 U(i,j,k)=B%mu_t(i,j,k)*Re                     ! mu_t
	  enddo
	  enddo
 	  enddo
    else                        ! 接收该块信息
	   Num_data=(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! 在该块中的编号
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, Struct_Comm,Status,ierr )
    endif
! write Data ....
    
	write(99) ((( U(i,j,k),i=0,nx),j=0,ny),k=0,nz)    

	deallocate(U)
   enddo
   close(99)
   close(100)
   close(101)
 
 else     ! 非0节点

    do m=1,MP%Num_Block     ! 本进程包含的块
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz))
	  Num_data=(nx+1)*(ny+1)*(nz+1)
 	  tag=m
	 
	   do k=0,nz
	   do j=0,ny
	   do i=0,nx
		 U(i,j,k)=B%mu_t(i,j,k)*Re
	   enddo
	   enddo
 	   enddo
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "write vt.dat OK"

  end subroutine output_vt



!----------------------------------------------------------------------
!  输出到壁面的距离  

  subroutine write_dw
   use Global_Var
   implicit none
   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:,:):: U
   integer:: NB,m,m1,nx,ny,nz,i,j,k,mt,Num_data
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)

   character(len=50):: filename   
   
     MP=>Mesh(1)
  

 if(my_id .eq. 0) then
   print*, "write wall_dist.dat ......"
   open(99,file="wall_dist.dat",form="unformatted")                    ! dw
  
   do m=1, Total_block   ! 全部块
 	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(nx-1,ny-1,nz-1))

	if(B_proc(m) .eq. 0) then            ! 这些块属于根进程
      mt=B_n(m)                          ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	 
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
        U(i,j,k)=B%dw(i,j,k)
	  enddo
	  enddo
 	  enddo
    else                        ! 接收该块信息
	   Num_data=(nx-1)*(ny-1)*(nz-1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! 在该块中的编号
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, Struct_Comm,Status,ierr )
    endif
! write Data ....
 	write(99) ((( U(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)    
	deallocate(U)
   enddo
   close(99)
 
 else     ! 非0节点

    do m=1,MP%Num_Block     ! 本进程包含的块
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(nx-1,ny-1,nz-1))
	  Num_data=(nx-1)*(ny-1)*(nz-1)
 	  tag=m
	 
	   do k=1,nz-1
	   do j=1,ny-1
	   do i=1,nx-1
		 U(i,j,k)=B%dw(i,j,k)
 	   enddo
	   enddo
 	   enddo
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "write wall_dist.dat OK"

  end subroutine write_dw




!----Boundary message (bc3d.inc, OpenCFD-EC Build-in format)----------------------------------------------------------
! OpenCFD-EC内建的.inc边界连接格式是在 Gridgen 的.inp格式基础上发展来的。
!  比.inp格式多了一些冗余信息，例如多了子面号f_no,
!  面类型face 以及连接的子面号f_no1,连接的面类型face1
!  以及连接次序L1, L2, L3  (例如L1=1表示该维与连接块的第1维正连接， L1=-1表示与连接块的第1为反向连接).
!  这些冗余信息为块-块之间的通信（尤其是MPI并行通信）提供了便利，有利于简化通信代码
!--------------------------------------------------------------------------------------------------------------------
  subroutine read_inc
   use Global_Var
   use mod_struct_grid, only: convert_inp_inc
   implicit none
   integer,parameter:: NC=21         ! .inc文件每单元有21个元素
   integer:: nx,ny,nz,NB,Nsub,m,mt,k,j
   integer:: Send_to_ID,tag,ierr,Status(MPI_Status_SIZE)
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer,pointer,dimension(:,:):: Bs 

 !  将Gridgen .inp 格式转化为 .inc格式   
	if(my_id .eq. 0) then
	  call  convert_inp_inc 
	endif


!  read bc3d.inc, 根进程读入，并发送至其他进程
   if(my_id .eq. 0) then
     print*, "read bc3d.inc (Link/boundary file)......"
     open(88,file="bc3d.inc")
     read(88,*)
     read(88,*) NB
 !   读取.inc文件中的元素   
    do m=1,NB
     read(88,*) nx,ny,nz
     read(88,*)
     read(88,*) Nsub   ! m块的子面数
 	 allocate(Bs(NC,Nsub))
	 do k=1,Nsub
	  read(88,*) (Bs(j,k),j=1,9)
	  read(88,*) (Bs(j,k),j=10,21)
!	 read(88,*) Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
!	 read(88,*) Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
     enddo
!     将该块信息发送出去
      
     if(B_Proc(m) .eq. 0) then
!          该块在0进程
       	 mt=B_n(m)   ! 在0进程中的内部编号
		 B=>Mesh(1)%Block(mt)
		 B%subface=Nsub   
	     allocate(B%bc_msg(B%subface))   ! 边界描述
         do k=1,Nsub
		 Bc=>B%bc_msg(k)
          Bc%ib=Bs(1,k); Bc%ie=Bs(2,k); Bc%jb=Bs(3,k); Bc%je=Bs(4,k)
		  Bc%kb=Bs(5,k); Bc%ke=Bs(6,k); Bc%bc=Bs(7,k); Bc%face=Bs(8,k); Bc%f_no=Bs(9,k)

          Bc%ib1=Bs(10,k); Bc%ie1=Bs(11,k); Bc%jb1=Bs(12,k); Bc%je1=Bs(13,k)
		  Bc%kb1=Bs(14,k); Bc%ke1=Bs(15,k); Bc%nb1=Bs(16,k); Bc%face1=Bs(17,k)
		  Bc%f_no1=Bs(18,k); Bc%L1=Bs(19,k); Bc%L2=Bs(20,k); Bc%L3=Bs(21,k)
         enddo
      else
!        将nsub 及Bs 发送出去
	     Send_to_ID=B_proc(m)              ! 发送目标块所在的进程号
	     tag=B_n(m)                        ! 标记
  	     call MPI_send(Nsub,1,MPI_INTEGER, Send_to_ID, tag, Struct_Comm,ierr )        !子面数
  	     call MPI_send(Bs,Nsub*Nc,MPI_INTEGER, Send_to_ID, tag, Struct_Comm,ierr )    !子面连接信息
     endif
	  deallocate(Bs)
    enddo
	  close(88)
   endif

! 非根进程
   if(my_id .ne. 0) then
      do m=1,Mesh(1)%Num_Block
	    B=>Mesh(1)%Block(m)
	    call MPI_Recv(Nsub,1,MPI_INTEGER,0,m,Struct_Comm,status,ierr)
 	    allocate(Bs(Nc,Nsub))
 	    call MPI_Recv(Bs,Nsub*Nc,MPI_INTEGER,0,m,Struct_Comm,status,ierr)
	    B%subface=Nsub   
      	allocate(B%bc_msg(B%subface))   ! 边界描述
 	     do k=1,Nsub
		  Bc=>B%bc_msg(k)
          Bc%ib=Bs(1,k); Bc%ie=Bs(2,k); Bc%jb=Bs(3,k); Bc%je=Bs(4,k)
		  Bc%kb=Bs(5,k); Bc%ke=Bs(6,k); Bc%bc=Bs(7,k); Bc%face=Bs(8,k); Bc%f_no=Bs(9,k)

          Bc%ib1=Bs(10,k); Bc%ie1=Bs(11,k); Bc%jb1=Bs(12,k); Bc%je1=Bs(13,k)
		  Bc%kb1=Bs(14,k); Bc%ke1=Bs(15,k); Bc%nb1=Bs(16,k); Bc%face1=Bs(17,k)
		  Bc%f_no1=Bs(18,k); Bc%L1=Bs(19,k); Bc%L2=Bs(20,k); Bc%L3=Bs(21,k)
         enddo
        deallocate(Bs)
	  enddo
	endif
	 call MPI_Barrier(Struct_Comm,ierr)
	 if(my_id .eq. 0) print*, "read bc3d.inc OK"

  end  subroutine read_inc 

!-----------------------------------------------------
!  计算到壁面的距离    
!  Code by Li Xinliang


 subroutine Comput_dist_wall
    use  Global_Var
    use  Wall_dist
    implicit none
    logical EXT
	integer:: Iext,ierr

!  判断 'wall_dist.dat'文件是否存在， 如存在则读取； 如不存在，则计算到壁面的距离

   if(my_id .eq. 0) then
     Inquire(file="wall_dist.dat",Exist=EXT)
 	 if(EXT) then
	   print*, "Find 'wall_dist.dat', read it ..."
	   Iext=1
	 else
	   print*, "Can not find the file 'wall_dist.dat', Comput wall distance ... "
	   Iext=0
     endif
   endif
   call MPI_Bcast(Iext,1,MPI_Integer,0,Struct_Comm,ierr)

   if(Iext .eq. 1) then
     call read_dw
   else
     call wall_point         ! 收集壁面网格点
     call comput_wall_dist          ! 计算各点到壁面的距离
     call write_dw           ! 写入文件
   endif

   end subroutine Comput_dist_wall
     
   subroutine comput_wall_dist
    use  Global_Var
    use  Wall_dist
    implicit none
    integer::  mb,i,j,k,k1,nx,ny,nz,ierr
	real(PRE_EC),allocatable,dimension(:,:,:):: dis
	real(PRE_EC):: d1, Dis_max,  Dis_min, Dis_max0,  Dis_min0
	real(PRE_EC),parameter:: d_init=1.d8
	Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    TYPE(BC_MSG_TYPE),pointer:: BC

     Dis_max=0.d0; Dis_min=d_init
     MP=>Mesh(1)
     do mb=1,MP%Num_Block
       B=>MP%Block(mb)
	    nx=B%nx; ny=B%ny; nz=B%nz
        allocate(dis(nx,ny,nz))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE (i,j,k,k1,d1)
		 do k=1,nz
          do j=1,ny
           do i=1,nx
		    dis(i,j,k)=d_init
		    do k1=1,Npw_total 
             d1=sqrt((B%x(i,j,k)-Xw(1,k1))**2+(B%y(i,j,k)-Xw(2,k1))**2+(B%z(i,j,k)-Xw(3,k1))**2)
             if(d1 .lt. dis(i,j,k)) dis(i,j,k)=d1
             enddo
            enddo
           enddo
		  enddo
 !$OMP END PARALLEL DO      
	   
     do k=1,nz-1
      do j=1,ny-1
       do i=1,nx-1
        B%dw(i,j,k)=(dis(i,j,k)+dis(i,j+1,k)+dis(i,j,k+1)+dis(i,j+1,k+1)+ &
            dis(i+1,j,k)+dis(i+1,j+1,k)+dis(i+1,j,k+1)+dis(i+1,j+1,k+1))*0.125    ! 格心点上的值
        Dis_max=max(Dis_max,B%dw(i,j,k))
        Dis_min=min(Dis_min,B%dw(i,j,k))
       enddo
      enddo
     enddo

    deallocate(dis)
   enddo
    call MPI_ALLREDUCE(Dis_max,Dis_max0,1,OCFD_DATA_TYPE,MPI_MAX,Struct_Comm,ierr)
    call MPI_ALLREDUCE(Dis_min,Dis_min0,1,OCFD_DATA_TYPE,MPI_MIN,Struct_Comm,ierr)
    if(my_id .eq. 0) then
	  print*, "----------comput wall dist OK ----------------"
	  print*, "Maximum distance to the wall is :" , Dis_max0
	  print*, "Minimum distance to the wall is :" , Dis_min0
	  print*, "----------------------------------------------"
	endif

    call MPI_Barrier(Struct_Comm,ierr)
   

  end subroutine comput_wall_dist
!-----------------------------------------
  
  
   
   
! 收集壁面上的网格点  
  subroutine wall_point
    use  Global_Var
    use  Wall_dist
    implicit none
    integer:: Npw1, mb,ksub,i,j,k,m,k0,ierr
	Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B
    TYPE(BC_MSG_TYPE),pointer:: BC
    integer,allocatable,dimension(:):: Npw,Nc,displs  ! 各进程壁面点的数目
    real(PRE_EC),allocatable,dimension(:,:):: Xw1  ! 壁面点（本进程）
    allocate(Npw(0:Total_proc-1),Nc(0:Total_proc-1),displs(0:Total_proc-1))

    MP=> Mesh(1)
!     统计壁面网格点的数目
      Npw1=0
      do mb=1,MP%Num_Block
          B=>Mp%Block(mb)
           do ksub=1,B%subface
             Bc=>B%bc_msg(ksub)
              if(Bc%bc .eq. BC_Wall) then
                Npw1=Npw1+(Bc%ie-Bc%ib+1)*(Bc%je-Bc%jb+1)*(Bc%ke-Bc%kb+1)
              endif
            enddo
	   enddo
       allocate(Xw1(3,Npw1))

!    读入本进程的壁面网格坐标
      k0=1
      do mb=1,MP%Num_Block
          B=>Mp%Block(mb)
           do ksub=1,B%subface
             Bc=>B%bc_msg(ksub)
              if(Bc%bc .eq. BC_Wall) then
                do k=Bc%kb,Bc%ke
                  do j=Bc%jb,Bc%je
                    do i=Bc%ib,Bc%ie
                      Xw1(1,k0)=B%x(i,j,k)
                      Xw1(2,k0)=B%y(i,j,k)
                      Xw1(3,k0)=B%z(i,j,k)
                      k0=k0+1
                    enddo
                   enddo
                 enddo
              endif
            enddo
	   enddo

! 全部进程的总壁面网格数
       call MPI_ALLREDUCE(Npw1,Npw_Total,1,MPI_INTEGER,MPI_SUM,Struct_Comm,ierr)
       allocate(Xw(3,Npw_total))
! 进行全搜集操作，得到全部的壁面网格坐标
       call MPI_Allgather(Npw1,1,MPI_Integer,Npw,1,MPI_Integer,Struct_Comm,ierr)
	   Nc(:)=3*Npw(:)   !数据量 
	   displs(0)=0
	   do m=1,Total_proc-1
	     displs(m)=displs(m-1)+Nc(m-1)
	   enddo
	   call MPI_Allgatherv(Xw1,3*Npw1,OCFD_DATA_TYPE,Xw,Nc,displs,OCFD_DATA_TYPE,Struct_Comm,ierr)
       deallocate(Xw1,Npw,Nc,displs)

     end

!  后处理模块,  计算力和力矩
!  A bug removed, 2017-3-13 
!--------------------------------------------------------
  subroutine comput_force      
   use Global_Var
   implicit none
   integer:: i,j,k,m,mB,nf,nx,ny,nz,NM,ierr
   real(PRE_EC):: Fx,Fy,Fz,Mx,My,Mz   ! 6分量
   real(PRE_EC):: Px,Py,Pz,Cfx,Cfy,Cfz   ! 压力和摩擦力
   real(PRE_EC):: fx0,fy0,fz0,xc,yc,zc,p1,p2
   real(PRE_EC):: P_inf,Pw  
   real(PRE_EC):: CL,CD,CS      ! 升力系数、阻力系数、侧向力系数
   real(PRE_EC),dimension(:),allocatable:: Mx1,My1,Mz1,Px1,Py1,Pz1,Cfx1,Cfy1,Cfz1    ! 各块的气动力、力矩
   real(PRE_EC),dimension(9):: Ft,Ft0  ! mpi reduce 汇总
   
   integer:: nMesh
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   character(len=50):: filename

!   print*, "comput force ..."    ! 如没有该语句在SW上运行出错 ???

   p_inf=1.d0/(gamma*Ma*Ma)
! 搜索细网格所有块的所有子面，如果发现固壁边界条件，则统计气动力及力矩
   NM=Mesh(1)%Num_Block
   allocate(Mx1(NM),My1(NM),Mz1(NM),Px1(NM),Py1(NM),Pz1(NM),Cfx1(NM),Cfy1(NM),Cfz1(NM))
   
   Mx1=0.d0; My1=0.d0; Mz1=0.d0
   Px1=0.d0; Py1=0.d0; Pz1=0.d0; Cfx1=0.d0; Cfy1=0.d0; Cfz1=0.d0
   
   do mB=1,NM
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz
     do nf=1,B%subface  
       Bc=> B%bc_msg(nf)

       if(Bc%bc .eq. BC_WALL) then
         if(Bc%face .eq. 1 ) then              ! i- 面
           do k=Bc%kb,Bc%ke-1
             do j=Bc%jb,Bc%je-1
!-------------------------------表面压力----------------------------------------------------  
               p1=(gamma-1.d0)*(B%U(5,1,j,k)-(B%U(2,1,j,k)**2+B%U(3,1,j,k)**2+B%U(4,1,j,k)**2)/B%U(1,1,j,k))
               p2=(gamma-1.d0)*(B%U(5,0,j,k)-(B%U(2,0,j,k)**2+B%U(3,0,j,k)**2+B%U(4,0,j,k)**2)/B%U(1,0,j,k))
               Pw=-(0.5d0*(p1+p2)-p_inf )*B%si(1,j,k)    ! 外法向
! -------------------------积分表面压力及表面摩擦阻力 -----------------------------------          
               Px1(mB)=Px1(mB)+Pw*B%ni1(1,j,k) ;   Py1(mB)=Py1(mB)+Pw*B%ni2(1,j,k) ;  Pz1(mB)=Pz1(mB)+Pw*B%ni3(1,j,k)   
!              粘性力 （i-, j-, k- 为正； i+ , j+, k+ 为负）
               Cfx1(mB)=Cfx1(mB)+B%Surf1(j,k,1) ;  Cfy1(mB)=Cfy1(mB)+B%Surf1(j,k,2) ;  Cfz1(mB)=Cfz1(mB)+B%Surf1(j,k,3)     

! --------------------------计算力矩 (矩心坐标 centroid(1:3) )------------------------            
               fx0=Pw*B%ni1(1,j,k)+B%Surf1(j,k,1)             
               fy0=Pw*B%ni2(1,j,k)+B%Surf1(j,k,2)
               fz0=Pw*B%ni3(1,j,k)+B%Surf1(j,k,3)
               xc=(B%xc(1,j,k)+B%xc(0,j,k))*0.5d0 -centroid(1)            
               yc=(B%yc(1,j,k)+B%yc(0,j,k))*0.5d0 -centroid(2)             
               zc=(B%zc(1,j,k)+B%zc(0,j,k))*0.5d0 -centroid(3)            
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo


         else if(Bc%face .eq. 2 ) then       ! j- 面
           do k=Bc%kb,Bc%ke-1
             do i=Bc%ib,Bc%ie-1
               p1=(gamma-1.d0)*(B%U(5,i,1,k)-(B%U(2,i,1,k)**2+B%U(3,i,1,k)**2+B%U(4,i,1,k)**2)/B%U(1,i,1,k))
               p2=(gamma-1.d0)*(B%U(5,i,0,k)-(B%U(2,i,0,k)**2+B%U(3,i,0,k)**2+B%U(4,i,0,k)**2)/B%U(1,i,0,k))
               Pw=-(0.5d0*(p1+p2)-p_inf )*B%sj(i,1,k)    
               Px1(mB)=Px1(mB)+Pw*B%nj1(i,1,k);   Py1(mB)=Py1(mB)+Pw*B%nj2(i,1,k) ;   Pz1(mB)=Pz1(mB)+Pw*B%nj3(i,1,k)   
               Cfx1(mB)=Cfx1(mB)+B%Surf2(i,k,1); Cfy1(mB)=Cfy1(mB)+B%Surf2(i,k,2);  Cfz1(mB)=Cfz1(mB)+B%Surf2(i,k,3)

 ! --------------------------计算力矩 (以坐标原点为中心) ---------------------------           
               fx0=Pw*B%nj1(i,1,k)+B%Surf2(i,k,1)                
               fy0=Pw*B%nj2(i,1,k)+B%Surf2(i,k,2)
               fz0=Pw*B%nj3(i,1,k)+B%Surf2(i,k,3)
               xc=(B%xc(i,1,k)+B%xc(i,0,k))*0.5d0 -centroid(1)            
               yc=(B%yc(i,1,k)+B%yc(i,0,k))*0.5d0 -centroid(2)            
               zc=(B%zc(i,1,k)+B%zc(i,0,k))*0.5d0 -centroid(3)            
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo        


         else if(Bc%face .eq. 3 ) then       ! k- 面
           do j=Bc%jb,Bc%je-1
             do i=Bc%ib,Bc%ie-1   
               p1=(gamma-1.d0)*(B%U(5,i,j,1)-(B%U(2,i,j,1)**2+B%U(3,i,j,1)**2+B%U(4,i,j,1)**2)/B%U(1,i,j,1))
               p2=(gamma-1.d0)*(B%U(5,i,j,0)-(B%U(2,i,j,0)**2+B%U(3,i,j,0)**2+B%U(4,i,j,0)**2)/B%U(1,i,j,0))
               Pw=-(0.5d0*(p1+p2)-p_inf )*B%sk(i,j,1)   
               Px1(mB)=Px1(mB)+Pw*B%nk1(i,j,1);  Py1(mB)=Py1(mB)+Pw*B%nk2(i,j,1) ;  Pz1(mB)=Pz1(mB)+Pw*B%nk3(i,j,1)   
               Cfx1(mB)=Cfx1(mB)+B%Surf3(i,j,1);  Cfy1(mB)=Cfy1(mB)+B%Surf3(i,j,2);  Cfz1(mB)=Cfz1(mB)+B%Surf3(i,j,3)         
			     
               fx0=Pw*B%nk1(i,j,1)+B%Surf3(i,j,1)          
               fy0=Pw*B%nk2(i,j,1)+B%Surf3(i,j,2)
               fz0=Pw*B%nk3(i,j,1)+B%Surf3(i,j,3)

               xc=(B%xc(i,j,1)+B%xc(i,j,0))*0.5d0  - centroid(1)
               yc=(B%yc(i,j,1)+B%yc(i,j,0))*0.5d0  - centroid(2)
               zc=(B%zc(i,j,1)+B%zc(i,j,0))*0.5d0  - centroid(3)
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo 

         else if(Bc%face .eq. 4 ) then              ! i+ 面
           do k=Bc%kb,Bc%ke-1
             do j=Bc%jb,Bc%je-1
 !-----------------------------------------表面压力 ---------------------------------------------------               
               p1=(gamma-1.d0)*(B%U(5,nx-1,j,k)-(B%U(2,nx-1,j,k)**2+B%U(3,nx-1,j,k)**2+B%U(4,nx-1,j,k)**2)/B%U(1,nx-1,j,k))
               p2=(gamma-1.d0)*(B%U(5,nx,j,k)-(B%U(2,nx,j,k)**2+B%U(3,nx,j,k)**2+B%U(4,nx,j,k)**2)/B%U(1,nx,j,k))
               Pw= (0.5d0*(p1+p2)-p_inf )*B%si(nx,j,k)  ! 外法向
! ------------------------------积分表面压力及表面摩擦阻力------------------------------------------           
               Px1(mB)=Px1(mB)+Pw*B%ni1(nx,j,k) ;   Py1(mB)=Py1(mB)+Pw*B%ni2(nx,j,k) ;  Pz1(mB)=Pz1(mB)+Pw*B%ni3(nx,j,k)   
!            粘性力 （i+, j+, k+ 面为负， 壁面所受力）
!               Cfx1(mB)=Cfx1(mB)+B%Surf4(j,k,1) ;  Cfy1(mB)=Cfy1(mB)+B%Surf4(j,k,2) ;  Cfz1(mB)=Cfz1(mB)+B%Surf4(j,k,3)         ! Bug removed
               Cfx1(mB)=Cfx1(mB)-B%Surf4(j,k,1) ;  Cfy1(mB)=Cfy1(mB)-B%Surf4(j,k,2) ;  Cfz1(mB)=Cfz1(mB)-B%Surf4(j,k,3)

! -------------------------------计算力矩 (以坐标原点为中心) --------------------------------------           
               fx0=Pw*B%ni1(nx,j,k)-B%Surf4(j,k,1)              ! Bug removed 
               fy0=Pw*B%ni2(nx,j,k)-B%Surf4(j,k,2)
               fz0=Pw*B%ni3(nx,j,k)-B%Surf4(j,k,3)
               xc=(B%xc(nx-1,j,k)+B%xc(nx,j,k))*0.5d0   - centroid(1)           
               yc=(B%yc(nx-1,j,k)+B%yc(nx,j,k))*0.5d0   - centroid(2)          
               zc=(B%zc(nx-1,j,k)+B%zc(nx,j,k))*0.5d0   - centroid(3)         
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo


         else if(Bc%face .eq. 5 ) then       ! j+ 面
           do k=Bc%kb,Bc%ke-1
             do i=Bc%ib,Bc%ie-1
               p1=(gamma-1.d0)*(B%U(5,i,ny-1,k)-(B%U(2,i,ny-1,k)**2+B%U(3,i,ny-1,k)**2+B%U(4,i,ny-1,k)**2)/B%U(1,i,ny-1,k))
               p2=(gamma-1.d0)*(B%U(5,i,ny,k)-(B%U(2,i,ny,k)**2+B%U(3,i,ny,k)**2+B%U(4,i,ny,k)**2)/B%U(1,i,ny,k))
               Pw= (0.5d0*(p1+p2)-p_inf )*B%sj(i,ny,k)
               Px1(mB)=Px1(mB)+Pw*B%nj1(i,ny,k);   Py1(mB)=Py1(mB)+Pw*B%nj2(i,ny,k) ;   Pz1(mB)=Pz1(mB)+Pw*B%nj3(i,ny,k)   
!              Cfx1(mB)=Cfx1(mB)+B%Surf5(i,k,1); Cfy1(mB)=Cfy1(mB)+B%Surf5(i,k,2);  Cfz1(mB)=Cfz1(mB)+B%Surf5(i,k,3)
               Cfx1(mB)=Cfx1(mB)-B%Surf5(i,k,1); Cfy1(mB)=Cfy1(mB)-B%Surf5(i,k,2);  Cfz1(mB)=Cfz1(mB)-B%Surf5(i,k,3)

! --------------------------------计算力矩 (以坐标原点为中心)  ---------------------------------------          
               fx0=Pw*B%nj1(i,ny,k)-B%Surf5(i,k,1)                 ! Bug removed
               fy0=Pw*B%nj2(i,ny,k)-B%Surf5(i,k,2)
               fz0=Pw*B%nj3(i,ny,k)-B%Surf5(i,k,3)
               xc=(B%xc(i,ny-1,k)+B%xc(i,ny,k))*0.5d0  - centroid(1)           
               yc=(B%yc(i,ny-1,k)+B%yc(i,ny,k))*0.5d0  - centroid(2)          
               zc=(B%zc(i,ny-1,k)+B%zc(i,ny,k))*0.5d0  - centroid(3)          
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo        


         else if(Bc%face .eq. 6 ) then       ! k+ 面
           do j=Bc%jb,Bc%je-1
             do i=Bc%ib,Bc%ie-1
               p1=(gamma-1.d0)*(B%U(5,i,j,nz-1)-(B%U(2,i,j,nz-1)**2+B%U(3,i,j,nz-1)**2+B%U(4,i,j,nz-1)**2)/B%U(1,i,j,nz-1))
               p2=(gamma-1.d0)*(B%U(5,i,j,nz)-(B%U(2,i,j,nz)**2+B%U(3,i,j,nz)**2+B%U(4,i,j,nz)**2)/B%U(1,i,j,nz))
               Pw= (0.5d0*(p1+p2)-p_inf )*B%sk(i,j,nz)  
               Px1(mB)=Px1(mB)+Pw*B%nk1(i,j,nz);  Py1(mB)=Py1(mB)+Pw*B%nk2(i,j,nz) ;  Pz1(mB)=Pz1(mB)+Pw*B%nk3(i,j,nz)   
!               Cfx1(mB)=Cfx1(mB)+B%Surf6(i,j,1);  Cfy1(mB)=Cfy1(mB)+B%Surf6(i,j,2);  Cfz1(mB)=Cfz1(mB)+B%Surf6(i,j,3)            ! Bug removed
               Cfx1(mB)=Cfx1(mB)-B%Surf6(i,j,1);  Cfy1(mB)=Cfy1(mB)-B%Surf6(i,j,2);  Cfz1(mB)=Cfz1(mB)-B%Surf6(i,j,3)

               fx0=Pw*B%nk1(i,j,nz)-B%Surf6(i,j,1)        ! Bug removed
               fy0=Pw*B%nk2(i,j,nz)-B%Surf6(i,j,2)
               fz0=Pw*B%nk3(i,j,nz)-B%Surf6(i,j,3)

               xc=(B%xc(i,j,nz-1)+B%xc(i,j,nz))*0.5d0 - centroid(1)
               yc=(B%yc(i,j,nz-1)+B%yc(i,j,nz))*0.5d0 - centroid(2)
               zc=(B%zc(i,j,nz-1)+B%zc(i,j,nz))*0.5d0 - centroid(3)
               Mx1(mB)=Mx1(mB)+ (yc*fz0-zc*fy0)
               My1(mB)=My1(mB)+ (zc*fx0-xc*fz0)
               Mz1(mB)=Mz1(mB)+ (xc*fy0-yc*fx0)
             enddo
           enddo
         endif
       endif
     enddo
   enddo
 
   Fx=0.d0; Fy=0.d0; Fz=0.d0; Mx=0.d0; My=0.d0; Mz=0.d0
   Px=0.d0; Py=0.d0; Pz=0.d0; Cfx=0.d0; Cfy=0.d0; Cfz=0.d0
 
 !  把各块的气动力、力矩加起来
   do mB=1,NM
   Px=Px+Px1(mB); Py=Py+Py1(mB); Pz=Pz+Pz1(mB)
   Cfx=Cfx+Cfx1(mB); Cfy=Cfy+Cfy1(mB); Cfz=Cfz+Cfz1(mB)
   Mx=Mx+Mx1(mB); My=My+My1(mB); Mz=Mz+Mz1(mB)
   enddo

   Ft(1)=Px ; Ft(2)=Py ; Ft(3)=Pz
   Ft(4)=Cfx; Ft(5)=Cfy; Ft(6)=Cfz
   Ft(7)=Mx;  Ft(8)=My; Ft(9)=Mz

!  各进程归约求和
   call MPI_ALLREDUCE(Ft,Ft0,9,OCFD_DATA_TYPE,MPI_SUM,Struct_Comm,ierr)

!   Fx=Px+Cfx; Fy=Py+Cfy; Fz=Pz+Cfz
    Ft0(1:6)=2.d0*Ft0(1:6)/Ref_S
	Ft0(7:9)=2.d0*Ft0(7:9)/(Ref_S*Ref_L)

	Fx=Ft0(1)+Ft0(4)              ! 气动力系数
	Fy=Ft0(2)+Ft0(5)
	Fz=Ft0(3)+Ft0(6)


!------   
 if(Cood_Y_UP ==1) then   ! Y 轴垂直向上 
	CL=Fy*cos(AoA)-Fx*sin(AoA)
	CD=Fx*cos(AoA)+Fy*sin(AoA)
    Cs=Fz
 else                      ! Z轴垂直向上 
	CL=Fz*cos(AoA)-Fx*sin(AoA)
	CD=Fx*cos(AoA)+Fz*sin(AoA)
    Cs=Fy
 endif


   if(my_id .eq. 0) then
   print*, "--------------Force and Moment -----------------------------"
   print*, "Total Force Coefficient (CL, CD, CS) ="
   print*, CL, CD, CS
   print*, "Total Moment Coefficient(CMx, CMy, CMz)="
   print*, Ft0(7),Ft0(8),Ft0(9)
   if(If_Debug==1) then
    print*, "Inviscous force: Fx, Fy, Fz="
    print*, Ft0(1),Ft0(2),Ft0(3)
    print*, "viscous force: Fx, Fy, Fz="
    print*, Ft0(4),Ft0(5),Ft0(6)
   endif



   open(99,file="force.log",position="append")
      write(99,"(I7,6E20.8)") Mesh(1)%Kstep, CL,CD,CS, Ft0(7),Ft0(8),Ft0(9)
   close(99)

!   if(If_Debug==1) then
      open(99,file="force-invis-vis.log",position="append")
      write(99,"(I7,6E20.8)") Mesh(1)%Kstep, Ft0(1), Ft0(2), Ft0(3), Ft0(4), Ft0(5), Ft0(6)
   close(99)
!   endif

   endif

   if(If_Debug == 1) then
      write(filename,"('force-debug-',3I3)") my_id
	  open(88,file=filename)
      do mB=1,NM
	     B => Mesh(1)%Block(mB)                                        
	    write(88,"(I4,6F20.10)")  B%Block_no,   Px1(mB), Py1(mB), Pz1(mB), Cfx1(mB), Cfy1(mB), Cfz1(mB)
	  enddo
      close(88)
   endif

	deallocate(Mx1,My1,Mz1,Px1,Py1,Pz1,Cfx1,Cfy1,Cfz1)


  end  subroutine comput_force   
        
!-----------------------------------------------------------------------------------
!  光顺（滤波） 操作， 耗散很大的滤波操作。 用于对初值的光顺，或者计算异常（如负温度）时的光顺
!  滤波运算可以消除高频振荡，提高计算的稳定性；但也会增加耗散，降低精度
!  2阶精度滤波耗散非常大，只能在处理初值或异常是使用，不可在常规的计算中使用。
!  4阶精度滤波也有一定耗散，计算过程中需谨慎使用 
  subroutine smoothing_oneMesh(nMesh,Smooth_method)
   use Global_Var
   use mod_struct_bc, only: Boundary_condition_onemesh
   use mod_struct_mpi, only: update_buffer_onemesh
   implicit none
   integer:: nMesh,mBlock,Smooth_method
!   print*, "Filtering ......", nMesh
   do mBlock=1,Mesh(nMesh)%Num_Block
   if(Smooth_method .eq. Smooth_2nd)then
    call   smoothing_oneBlock_2nd(nMesh,mBlock)     
   else
    call   smoothing_oneBlock_4th(nMesh,mBlock)     
   endif
   enddo
   call Boundary_condition_onemesh(nMesh)             ! 边界条件 （设定Ghost Cell的值）
   call update_buffer_onemesh(nMesh)                  ! 同步各块的交界区
   end subroutine smoothing_oneMesh

!----------------------------------------------------------
! 低精度滤波（2阶精度）
  subroutine smoothing_oneBlock_2nd(nMesh,mBlock)     
   use Global_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1,istat
   Type (Block_TYPE),pointer:: B
   real(PRE_EC)::tmpa(0:2000),tmpb(0:2000),tmpc(0:2000)
   NVAR1=Mesh(nMesh)%NVAR
   B=>Mesh(nMesh)%block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
!----------------------------------------------------------------------------
! Warning, allocatable不能作为私有变量 !!!   
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)

!$OMP DO   
   do k=1,nz-1
   do j=1,ny-1
   do m=1, NVAR1
    do i=0,nx
    tmpa(i)=B%U(m,i,j,k)
    enddo
   do i=1,nx-1
     B%U(m,i,j,k)=0.25d0*(tmpa(i-1)+tmpa(i+1))+0.5d0*tmpa(i)
   enddo
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do k=1,nz-1
   do i=1,nx-1
   do m=1, NVAR1
   do j=0,ny
    tmpb(j)=B%U(m,i,j,k)
   enddo
   do j=1,ny-1
     B%U(m,i,j,k)=0.25d0*(tmpb(j-1)+tmpb(j+1))+0.5d0*tmpb(j)
   enddo
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do j=1,ny-1
   do i=1,nx-1
   do m=1, NVAR1
   do k=0,nz
    tmpc(k)=B%U(m,i,j,k)
   enddo
   do k=1,nz-1
     B%U(m,i,j,k)=0.25d0*(tmpc(k-1)+tmpc(k+1))+0.5d0*tmpc(k)
   enddo
   enddo
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!-----------------------------------------------------------------------------------------------
   end  subroutine smoothing_oneBlock_2nd
   
   
!------------------------------------------------------------------------------------------------
! 高精度滤波（4阶精度）
  subroutine smoothing_oneBlock_4th(nMesh,mBlock)     
   use Global_Var
   implicit none
   integer:: nMesh,mBlock
   integer:: i,j,k,m,nx,ny,nz,NVAR1
   Type (Block_TYPE),pointer:: B
   real(PRE_EC)::tmpa(0:2000),tmpb(0:2000),tmpc(0:2000)

   NVAR1=Mesh(nMesh)%NVAR
   B=>Mesh(nMesh)%block(mBlock)
   nx=B%nx; ny=B%ny; nz=B%nz
   
!---------------------------------------------------------------------------- 
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,nz,NVAR1,B)

!$OMP DO   
   do k=1,nz-1
   do j=1,ny-1
   do m=1, NVAR1
    do i=0,nx
    tmpa(i)=B%U(m,i,j,k)
    enddo
    do i=2,nx-2
     B%U(m,i,j,k)=(-tmpa(i-2)+4.d0*tmpa(i-1)+10.d0*tmpa(i)+4.d0*tmpa(i+1)-tmpa(i+2))/16.d0
    enddo
     B%U(m,1,j,k)=0.25d0*(tmpa(0)+tmpa(2))+0.5d0*tmpa(1)
     B%U(m,nx-1,j,k)=0.25d0*(tmpa(nx-2)+tmpa(nx))+0.5d0*tmpa(nx-1)
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do k=1,nz-1
   do i=1,nx-1
   do m=1, NVAR1
   do j=0,ny
    tmpb(j)=B%U(m,i,j,k)
   enddo
    do j=2,ny-2
     B%U(m,i,j,k)=(-tmpb(j-2)+4.d0*tmpb(j-1)+10.d0*tmpb(j)+4.d0*tmpb(j+1)-tmpb(j+2))/16.d0
    enddo
     B%U(m,i,1,k)=0.25d0*(tmpb(0)+tmpb(2))+0.5d0*tmpb(1)
     B%U(m,i,ny-1,k)=0.25d0*(tmpb(ny-2)+tmpb(ny))+0.5d0*tmpb(ny-1)
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP DO   
   do j=1,ny-1
   do i=1,nx-1
   do m=1, NVAR1
   do k=0,nz
    tmpc(k)=B%U(m,i,j,k)
   enddo
    do k=2,nz-2
     B%U(m,i,j,k)=(-tmpc(k-2)+4.d0*tmpc(k-1)+10.d0*tmpc(k)+4.d0*tmpc(k+1)-tmpc(k+2))/16.d0
    enddo
     B%U(m,i,j,1)=0.25d0*(tmpc(0)+tmpc(2))+0.5d0*tmpc(1)
     B%U(m,i,j,nz-1)=0.25d0*(tmpc(nz-2)+tmpc(nz))+0.5d0*tmpc(nz-1)
   enddo
   enddo
   enddo
!$OMP END DO

!$OMP END PARALLEL
!-----------------------------------------------------------------------
   end  subroutine smoothing_oneBlock_4th

!  后处理模块： 进行时间平均
!--------------------------------------------------------
  subroutine Time_average      
   use Global_Var
   implicit none
   integer:: i,j,k,m,mB,nf,nx,ny,nz
   integer,save:: Iflag=0
   real(PRE_EC):: d1,u1,v1,w1,T1,tmp

   Type (Block_TYPE),pointer:: B

!---仅被运行1次 -----------------
   if(Iflag == 0) then
     Iflag=1             
    do mB=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz
 	 allocate(B%U_average(0:nx,0:ny,0:nz,5))       ! 时均量 d,u,v,w,T
     enddo

    call init_average        ! 初始化平均场

   endif
!-------------------------------------
! 时间平均  
   if(my_id .eq. 0) print*, "Time Average ......", Istep_average+1

   tmp=1.d0/(Istep_average+1.d0)
   do mB=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz

!$OMP PARALLEL DO PRIVATE(i,j,k,d1,u1,v1,w1,T1) SHARED(nx,ny,nz,B,Cv,tmp,Istep_average)

     do k=0,nz
	 do j=0,ny
	 do i=0,nx
         d1= B%U(1,i,j,k)
         u1= B%U(2,i,j,k)/d1
         v1= B%U(3,i,j,k)/d1
         w1= B%U(4,i,j,k)/d1
         T1=(B%U(5,i,j,k)-0.5d0*d1*(u1*u1+v1*v1+w1*w1))/(Cv*d1)
      
	   B%U_average(i,j,k,1)=(Istep_average*B%U_average(i,j,k,1)+d1)*tmp    
	   B%U_average(i,j,k,2)=(Istep_average*B%U_average(i,j,k,2)+u1)*tmp    
	   B%U_average(i,j,k,3)=(Istep_average*B%U_average(i,j,k,3)+v1)*tmp    
	   B%U_average(i,j,k,4)=(Istep_average*B%U_average(i,j,k,4)+w1)*tmp    
	   B%U_average(i,j,k,5)=(Istep_average*B%U_average(i,j,k,5)+T1)*tmp    

	 enddo
	 enddo
	 enddo
   enddo
    Istep_average=Istep_average+1
  end

!-----------------------------------------------------
! 初始化, 目前版本只支持重新开始平均，暂不支持读取flow3d_average.dat
   subroutine init_average      
   use Global_Var
   implicit none
   integer:: i,j,k,mB,nf,nx,ny,nz

   Type (Block_TYPE),pointer:: B
    Istep_average=0
    do mB=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(mB)                                        
     nx=B%nx; ny=B%ny; nz=B%nz
      do k=0,nz
	  do j=0,ny
	  do i=0,nx
	   B%U_average(i,j,k,1)=0.d0 
	   B%U_average(i,j,k,2)=0.d0    
	   B%U_average(i,j,k,3)=0.d0    
	   B%U_average(i,j,k,4)=0.d0    
	   B%U_average(i,j,k,5)=0.d0    
	  enddo
	  enddo
	  enddo
     enddo
  end
!----------------------------------------------------      

 !  输出平均量 （Plot3d格式）, 最细网格flow3d_average.dat  
  subroutine output_flow_average
   use Global_Var
   implicit none
   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:,:,:):: U
   integer:: NB,m,m1,nx,ny,nz,i,j,k,mt,Num_data
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)

   character(len=50):: filename   
   
   MP=>Mesh(1)

!---------------------------------------------------------------
 if(my_id .eq. 0) then
  
   print*, "write flow3d_average.dat ......"
   
   open(99,file="flow3d_average.dat",form="unformatted")                    ! d,u,v,w,T

   do m=1, Total_block   ! 全部块
     
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz,5))

	if(B_proc(m) .eq. 0) then             ! 这些块属于根进程
      mt=B_n(m)                           ! 该块在进程内部的编号
	  B=>MP%Block(mt)
	 
	   do m1=1,5
	    do k=0,nz
	    do j=0,ny
	    do i=0,nx
		  U(i,j,k,m1)=B%U_average(i,j,k,m1)                     
        enddo
	    enddo
	    enddo
 	   enddo

    else                        ! 接收该块信息
	   Num_data=5*(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! 在该块中的编号
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, Struct_Comm,Status,ierr )
    endif
! write Data ....
    
	write(99) (((( U(i,j,k,m1),i=0,nx),j=0,ny),k=0,nz),m1=1,5)    
    
	deallocate(U)

   enddo
   
   write(99) Istep_average
   close(99)
 
 else     ! 非0节点

    do m=1,MP%Num_Block     ! 本进程包含的块
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz,5))
	  Num_data=(nx+1)*(ny+1)*(nz+1)*5
 	  tag=m
	 
      do m1=1,5
	    do k=0,nz
	    do j=0,ny
	    do i=0,nx
		 U(i,j,k,m1)=B%U_average(i,j,k,m1)
        enddo
	    enddo
	    enddo
 	   enddo
	   
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, Struct_Comm,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(Struct_Comm,ierr)
   if(my_id .eq. 0)  print*, "write flow3d_average.dat OK"

  end subroutine output_flow_average

end module mod_struct_io
