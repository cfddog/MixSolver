!===============================================================================
! mod_struct_mpi.f90 -- structured solver: MPI communication / partitioning
! Encapsulates sub_partation_mpi + sub_update_buffer_mpi + sub_udate_Meshlink_mpi.
! Note: 'use interface_defines' lines removed -- get_U_conner/get_xyz_conner are
! now module procedures in this file; explicit interfaces come free.
!===============================================================================
  module mod_struct_mpi
  contains
!---------------------------------------------------------
!  分区： 确定 “进程”与“块”之间的联系
!  B_proc(m), 给出m块所属的进程号 
!  B_n(m), 给出m块在该进程中的内部编号
!  my_blocks(:), 给出本进程所包含的块号
!    
   subroutine partation
   use Global_var
   implicit none
   integer:: ierr,m,mt,k,TP
   integer,dimension(:),pointer:: Bproc
   logical EXT

   if(my_id .eq. 0) then
    Inquire(file="partation.dat",EXIST=EXT)
     if(EXT) then
	   print*, "Find partation.dat, read it"
	   open(100,file="partation.dat")
	   read(100,*)
       read(100,*)  Total_block, TP          ! 总块数
	   read(100,*) 
       if(TP .ne. Total_proc) then
        print*, "Error ! Total_proc is not the same as that in  'partation.dat' "
        print*, "Tp, Total_proc=", Tp, Total_Proc
        stop
       endif
	 else 
	    print*, "can not find partation.dat,  create partation ..."
        call part(Mesh_File_Format,Total_Proc,Total_Block,Bproc)
	 endif
    endif

	 call MPI_bcast(Total_block,1,MPI_Integer,0,  Struct_Comm,ierr)
     allocate (B_proc(Total_block),B_n(Total_block))   
!   读入B_proc(:)  , 块-->进程 对应关系     
	 
	 if(my_id .eq. 0) then
     if(EXT) then
	   do m=1,  Total_block 
	   read(100,*) mt, B_proc(m)
	   enddo
	   close(100)
	  else
        B_proc(1:Total_Block)=Bproc(1:Total_Block)
!		print*, "------------------"
!		print*, B_proc(:)
	  endif
	  
	  endif

     call MPI_bcast(B_proc(1),Total_block,MPI_Integer,0,  Struct_Comm,ierr)

!  计算B_n(m), 第m块在该进程的内部编号
    do m=1, Total_block
	 B_n(m)=0
	 do k=1,m                 ! m块前面有多数个块在B_proc(m)进程
	  if(B_proc(k) .eq. B_proc(m))  B_n(m)=B_n(m)+1
	 enddo
	enddo 


!  统计每个进程包含的块数
     Num_block=0
	 do m=1, Total_block
	 if(B_proc(m) .eq. my_id )   Num_block=Num_block+1    ! my_id进程包含的块数
     enddo
     
	 allocate (my_blocks(Num_block))       ! 本进程包含的块号（数组）
	 k=1
	 do m=1, Total_block
  	  if(B_proc(m) .eq. my_id )  then
	    my_blocks(k)=m
		k=k+1
	  endif
     enddo

   end subroutine 


     subroutine part(formt,Num_Proc,Num_Block,Bproc)
     implicit none
     integer:: formt,Num_Proc,Num_Block
	 integer,dimension(:),pointer:: Bproc
     integer,allocatable,dimension(:,:):: B_grid
     integer,allocatable,dimension(:):: NI,NJ,NK,Pgrid
	 integer:: m,m0,m1,m2,G0,mb,mg,mp,n,t1,t2,k,Total_cell,Total_grid
	 real:: grid_av
     print*, "Partation ......"
	 if(formt .eq. 0) then
       open(99,file="Mesh3d.x",form="unformatted")
       read(99) Num_Block         
       allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! Size of each block
       read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
	 else
       open(99,file="Mesh3d.x")
       read(99,*) Num_Block         ! 
       allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! Size of each block
       read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
     endif
     close(99)

 !    读入每块网格数，按从多到少次序排列
     allocate(B_grid(Num_Block,2))
	 allocate(Pgrid(Num_proc),Bproc(Num_Block))

     Total_Cell=0
	 Total_grid=0
	  do m=1,Num_Block
       B_grid(m,1)=NI(m)*NJ(m)*NK(m) 
	   B_grid(m,2)=m
	   Total_Cell=Total_Cell+(NI(m)-1)*(NJ(m)-1)*(NK(m)-1)
	   Total_grid=Total_grid+NI(m)*NJ(m)*NK(m)
      enddo
	  print*, "Total Cell number=", Total_Cell
	  print*, "Total Grid number=", Total_grid
 
!    按网格点从多到少的次序排序
     do m=1,Num_block
	   G0=B_grid(m,1)     ! 点数
	    mg=m
	    do n=m+1,Num_block   ! 找出数目最大的
         if(B_grid(n,1) .gt. G0 ) then
		    G0=B_grid(n,1)
			mg=n
         endif
       enddo
!        mg块与m块交换
         t1=B_grid(mg,1)
		 t2=B_grid(mg,2)
		 B_grid(mg,1)=B_grid(m,1)
		 B_grid(mg,2)=B_grid(m,2)
         B_grid(m,1)=t1
		 B_grid(m,2)=t2
	  enddo
!----------------------------
      Pgrid(:)=0
	 do m=1,Num_Block
        mb=B_grid(m,2) 
  !     寻找网格数目最小的进程
        mg=Pgrid(1)
		m0=1
		do mp=1,Num_Proc
         if(Pgrid(mp) .lt. mg) then
		   mg=Pgrid(mp)
		   m0=mp
		 endif
		enddo
       Pgrid(m0)=Pgrid(m0)+B_grid(m,1)
	   Bproc(mb)=m0-1
     enddo
!-------------输出--------------------------------------------------
     open(99,file="partation-auto.dat")
	 write(99,*) " Total_block_number    Total_Proc_number "
	 write(99,*)  Num_block, Num_proc
	 write(99,*) " Block_no, Proc_no "
	 do m=1,Num_block
	 write(99,*) m, Bproc(m)
	 enddo
	 close(99)

     open(100,file="part_grid.dat")
	 do m=1,Num_proc
	  write(100,*) m-1, Pgrid(m)
	 enddo
     close(100)

     m1=Pgrid(1)
	 m2=m1
     do mp=1,Num_Proc
      m1=max(m1,Pgrid(mp))
	  m2=min(m2,Pgrid(mp))
	 enddo
	 
	 grid_av=1.d0*Total_grid/Num_Proc
	 print*, "Max grid number is ", m1,  " ... rato to mean grid is", m1/grid_av
	 print*, "Min grid number is ", m2,  " ... rato to mean grid is", m2/grid_av
	 print*, "rato max to min is ", 1.0*m1/m2
   deallocate(NI,NJ,NK,Pgrid,B_grid)
   end subroutine part




!----多块网格之间的连接--------------------------------------------------------------------------
! Copyright by Li Xinliang
! Ver 0.41 修改了ori的定义方法，以兼容BXCFD
! Ver 0.43b  修改了角点的计算方法 
! Ver 0.46a 改正了计算连接方向子程序 get_ijk_orient()的一个重要Bug   (face1与face2弄反了)
! Ver 0.50  (2010-12-9):  物理量采用双层Ghost Cell, 坐标仍采用单层Ghost Cell
! Ver 0.84  (2012-7-10):  基于Gridgen 的.inp格式的连接
! Ver 0.9-mpi  MPI并行版 （以往为串行版）

!---------Continue boundary (inner boundary) ------------------------------------------------------
! 传送相邻面的物理量（守恒变量）给本块的虚网格； 
! 由于坐标存储在节点，物理量存储在中心点， 因而传输坐标与传输物理量的方法（下标对应方式）略有区别
! MPI并行版； 
!  策略： 如两块(block)位于同一进程(proc),则采用直接传递(不通过MPI),以提高效率；
!  使用MPI传递时，首先使用MPI_Bsend()发送全部信息； 然后使用MPI_recv()接收。 

! 用于多重网格时，也许会有问题 （应当NVAR1） , ...
     subroutine update_buffer_onemesh(nMesh)
     use Global_Var
     implicit none
     integer:: nMesh,ierr
! --------------------------------------------------------------------------------------- 
! 模块边界通信：  MPI 版本
    call Umessage_send_mpi(nMesh)   ! 使用MPI发送全部信息 ；如目标块也在本进程内，则不通过MPI,直接交换信息.
    call Umessage_recv_mpi(nMesh)   ! 接收全部信息
   
    call MPI_Barrier(Struct_Comm,ierr)
    
    
	call Umessage_corner(nMesh)     ! 角区数据信息（采用插值的方法）

  end subroutine update_buffer_onemesh

!---------------------------------------------------------------

! 同一进程内两块之间的通信 （可直接通信）

     subroutine Umessage_send_mpi(nMesh)  
     use Global_Var
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,i1,j1,k1,mb,mBlock,ksub,nMesh,Send_to_ID,Num_Data,tag,ierr
     integer:: kb(3),ke(3),kb1(3),ke1(3),ka(3),L1,L2,L3,P(3),Ks(3)
     real(PRE_EC),allocatable,dimension(:,:,:,:):: Usend

 ! --------------------------------------------------------------------------------------- 
 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)

  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)

!-----------------------------------------------------------------------------------------------
      if(Bc%bc .ge. 0 ) cycle          ! 非内边界 inner boundary
      
!--------------------------------------------------------
!  B块为源数据； B1块为目标数据；   数据从源数据写入目标数据

!      源数据:临近边界的LAP层 内点;   i=kb(1) to ke(1), j=kb(2) to ke(2),  k= kb(3) to ke(3)
!      目标数据 （临近边界的LAP层 Ghost点,  现版本LAP=2）

!      源数据 
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                 ! k=1,2,3 为i,j,k方向
	   if(Bc%face .gt. 3) kb(k)=kb(k)-LAP   ! i+, j+ or k+ 面
       ke(k)=kb(k)+LAP-1

!      目标数据 （边界外扩充的LAP层Ghost点）
       kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
       k=mod(Bc%face1-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face1 .le. 3)  kb1(k)=kb1(k)-LAP   ! i-, j- or k- 面
       ke1(k)=kb1(k)+LAP-1


!  L1,L2,L3 : 描述维的连接 ; P(1),P(2),P(3): 描述连接次序(顺序or 逆序)
       L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)    ! L1=1 意味着（本维）与目标数据的第1维连接， =2 为与第2维连接, ...
       L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
       L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)

!   ks(k) : 目标数据第k维的起始下标
	   do k=1,3  !   目标数据 起始下标     
	    if(P(k) .gt. 0) then
		 ks(k)=kb(k)     ! 顺序，从kb开始
		else
		 ks(k)=ke(k)     ! 拟序， 从ke开始
	    endif
	  enddo
!
!-------------------------------------------------------------
      allocate(Usend(NVAR,kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
       do k=kb1(3),ke1(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
		     ka(1)=i-kb1(1)            
			 ka(2)=j-kb1(2)
			 ka(3)=k-kb1(3)
			 i1=ks(1)+ka(L1)*P(1)        ! 源数据的下标 (i1,j1,k1), 从kb1到ke1 (需考虑 a.顺序或逆序, b. 维的连接), L1,L2,L3控制维的连接，P(k)控制顺序/逆序
			 j1=ks(2)+ka(L2)*P(2)
			 k1=ks(3)+ka(L3)*P(3)
			 do m=1,NVAR
			  Usend(m,i,j,k)=B%U(m,i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
             enddo
           enddo
		 enddo
		enddo
  
  !    
	 Send_to_ID=B_proc(Bc%nb1)              ! 发送目标块所在的进程号

	 if( Send_to_ID .ne. my_id) then        ! 目标块不在本进程内
       Num_data=NVAR*(ke1(1)-kb1(1)+1)*(ke1(2)-kb1(2)+1)*(ke1(3)-kb1(3)+1)      ! 数据量
       tag=  Bc%nb1*1000+Bc%f_no1                                         ! 标记; 块号Bc%nb1, 子面号Bc%f_no1 （向一个进程发送多个数据包时，用于识别）
	   call MPI_Bsend(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
!	   call MPI_send(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
	 else 
!        目标块也在本进程内，直接写入 (不通过MPI)
         mb=B_n(Bc%nb1)          ! 块内部编号
		 B1 =>  Mesh(nMesh)%Block(mb)  ! 相邻块 （目标块）
	     do k=kb1(3),ke1(3)             ! 
	      do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
			 do m=1,NVAR
			  B1%U(m,i,j,k)=Usend(m,i,j,k)    
             enddo
            enddo
		  enddo
		 enddo
	  endif	
	  deallocate(Usend)

  enddo
  enddo
 
  end
!--------------------------------------------------------------------------

! 将全部信息接收

subroutine Umessage_recv_mpi(nMesh) ! 使用MPI发送全部信息
     use Global_Var
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,Recv_from_ID,kb(3),ke(3)
	 integer:: tag,Num_DATA,ierr,Status(MPI_Status_SIZE)
     real(PRE_EC),allocatable,dimension(:,:,:,:):: Urecv

!---------------------------------------------------------------------------------------
! 接收阶段，将全部信息接收
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ge. 0 ) cycle               ! 内边界 inner boundary
   Recv_from_ID=B_proc(Bc%nb1)           ! 相邻块（接收源块）所在的进程号
   if(Recv_from_ID .eq. my_id) cycle     ! 源块在本进程内，不使用MPI通信 (Umessage_send_mpi()已完成写入操作)  
   
!------------------------------------------------------------------------------------------------
!      目标数据 （边界外扩充的LAP层Ghost点）
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3)  kb(k)=kb(k)-LAP   ! i-, j- or k- 面
       ke(k)=kb(k)+LAP-1

       allocate(Urecv(NVAR,kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))    ! 接收数组，按照目标数据的格式
       Num_data=NVAR*(ke(1)-kb(1)+1)*(ke(2)-kb(2)+1)*(ke(3)-kb(3)+1)      ! 数据量
	   tag=B%Block_no*1000+Bc%f_no                                  ! tag 标记，标记块号+子面号
	   
	   call MPI_Recv(Urecv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,Struct_Comm,status,ierr)

        do k=kb(3),ke(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
			 do m=1,NVAR
             B%U(m,i,j,k)=Urecv(m,i,j,k)
			 enddo
           enddo
		  enddo
		enddo
       deallocate(Urecv)

   enddo
   enddo
  end


!------------------------------------------------------------------------
 subroutine Umessage_corner(nMesh)   ! 角区数据信息（采用插值的方法）
   use Global_Var
   implicit none
   integer:: nMesh,mBlock
   Type (Block_TYPE),pointer:: B
 
  do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  !   角点坐标采用外插形式获得
  ! Visual Fortran 与Intel Fortran不兼容， 二者用不同的语法
   call get_U_conner(B%nx,B%ny,B%nz,NVAR,B%U)          ! 使用Intel Fortran
  
  enddo
  end
 

!---------------------------------------------------------------------------------
  subroutine Update_coordinate_buffer
    use Global_Var
    implicit none
    integer nMesh
	do nMesh=1,Num_Mesh
  	  call Update_coordinate_buffer_onemesh(nMesh)
	enddo
   end subroutine Update_coordinate_buffer
!-----------------------------------------------------------------


!---------Continue boundary (inner boundary) -------------------------------
! 传送相邻界面的坐标给本块的虚网格； 
! 由于坐标存储在节点，物理量存储在中心点， 因而传输坐标与传输物理量的方法（下标对应方式）略有区别
! 坐标采用1层虚网格点
    subroutine Update_coordinate_buffer_onemesh(nMesh)
     use Global_Var
     implicit none
     integer:: nMesh,ierr
	 call Coordinate_send_mpi(nMesh)
     call Coordinate_recv_mpi(nMesh)
     call MPI_Barrier(Struct_Comm,ierr)

     call Coordinate_Periodic(nMesh)              ! 周期边界的网格
     call coordinate_boundary_and_corner(nMesh)   ! 非内边界的网格
     call MPI_Barrier(Struct_Comm,ierr)

 !---test----------------------------------------------------------
   end 





    subroutine Coordinate_send_mpi(nMesh)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc,Bc1
     real(PRE_EC),allocatable:: Ux_send(:,:,:,:)   ! 存放界面上的数据（坐标）
     integer:: i,j,k,m,i1,j1,k1,i2,j2,k2,km,mBlock,ksub,m_neighbour,msub,nMesh,mb
	 integer:: Send_to_ID,Num_Data,tag,ierr
     integer:: kb(3),ke(3),kb1(3),ke1(3),Ka(3),L1,L2,L3,P(3),ks(3)
 ! ---------------------------------------------------------------------------------------- 
  do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)
   do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)

   if(Bc%bc .ge. 0 ) cycle    ! 仅处理内边界 inner boundary
!--------------------------------------------------------------------------------------------  
!    将源数据 拷贝 到目标数据； 

!      源数据 （临近边界的1层内点）
       kb(1)=Bc%ib; ke(1)=Bc%ie; kb(2)=Bc%jb; ke(2)=Bc%je; kb(3)=Bc%kb; ke(3)=Bc%ke
       k=mod(Bc%face-1,3)+1         ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3) then      ! i-, j- or k- 方向
	    kb(k)=kb(k)+1       !  i=2 (内点)
	   else
	    kb(k)=kb(k)-1       !  i=nx-1 (内点) 
	   endif
        ke(k)=kb(k)

!      目标数据 （边界外的1层 Ghost）
       kb1(1)=Bc%ib1; ke1(1)=Bc%ie1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1
       k=mod(Bc%face1-1,3)+1
	   if(Bc%face1 .le. 3) then          ! i-, j- or k-
	    kb1(k)=kb1(k)-1                  !  i=0 (Ghost点)
	   else
	    kb1(k)=kb1(k)+1                  ! i=nx+1 (Ghost点)
	   endif
	    ke1(k)=kb1(k)
	   

       L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)
       L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
       L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)
 
 !   目标数据 起始下标     
	   do k=1,3
	    if(P(k) .gt. 0) then
		 ks(k)=kb(k)     ! 顺序，从kb开始
		else
		 ks(k)=ke(k)     ! 拟序， 从ke开始
	    endif
	  enddo
!---------------------------------------------------------
!    开辟数组，存放坐标数据 （以目标数据的格式）
      allocate(Ux_send(3,kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))

!----------------------------------------------------------
! 从源数据 拷贝到 目标数据      
       
       do k=kb1(3),ke1(3)
	     do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
		     ka(1)=i-kb1(1)
			 ka(2)=j-kb1(2)
			 ka(3)=k-kb1(3)
			 i1=ks(1)+ka(L1)*P(1)
			 j1=ks(2)+ka(L2)*P(2)
			 k1=ks(3)+ka(L3)*P(3)
			 Ux_send(1,i,j,k)=B%x(i1,j1,k1)    ! 从源数据（内点） 拷贝入 临时数组
             Ux_send(2,i,j,k)=B%y(i1,j1,k1)
             Ux_send(3,i,j,k)=B%z(i1,j1,k1)
           enddo
		 enddo
		enddo

    	 Send_to_ID=B_proc(Bc%nb1)              ! 发送目标块所在的进程号
    
   
	 if( Send_to_ID .ne. my_id) then        ! 目标块不在本进程内
       Num_data=3*(ke1(1)-kb1(1)+1)*(ke1(2)-kb1(2)+1)*(ke1(3)-kb1(3)+1)      ! 数据量
       tag=  Bc%nb1*1000+Bc%f_no1                                         ! 标记; 块号Bc%nb1, 子面号Bc%f_no1 （向一个进程发送多个数据包时，用于识别）
	   call MPI_Bsend(Ux_send,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
!	   call MPI_send(Ux_send,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
	 else 
       mb=B_n(Bc%nb1)    ! Bc%nb1 块在该进程中的内部编号
       B1 =>  Mesh(nMesh)%Block(mb)    ! 相邻界面
	   do k=kb1(3),ke1(3)
	     do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
  		    B1%x(i,j,k)=Ux_send(1,i,j,k)    ! 临时数组拷贝入目标数据（Ghost 点）
            B1%y(i,j,k)=Ux_send(2,i,j,k)
            B1%z(i,j,k)=Ux_send(3,i,j,k)
           enddo
		 enddo
	   enddo
      endif
	 deallocate(Ux_send)
  
  enddo
  enddo

 end subroutine Coordinate_send_mpi


!---------------------------------------------------------------------------------------
! 接收（坐标信息）
    subroutine Coordinate_recv_mpi(nMesh)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     real(PRE_EC),allocatable:: Ux_recv(:,:,:,:)   ! 存放界面上的数据（坐标）
     integer:: i,j,k,mBlock,ksub,nMesh,Recv_from_ID,Num_data,tag,ierr,Status(MPI_Status_SIZE)

     integer:: kb(3),ke(3)
 ! ---------------------------------------------------------------------------------------- 
  do mBlock=1,Mesh(nMesh)%Num_Block

   B => Mesh(nMesh)%Block(mBlock)
   do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ge. 0 ) cycle    ! 仅处理内边界 inner boundary
!--------------------------------------------------------------------------------------------  
    Recv_from_ID=B_proc(Bc%nb1)           ! 相邻块（接收源块）所在的进程号
    if(Recv_from_ID .eq. my_id) cycle     ! 源块在本进程内，不使用MPI通信 ( Coordinate_send_mpi()已完成写入操作)  

!      目标数据 （边界外的1层 Ghost）
       kb(1)=Bc%ib; ke(1)=Bc%ie; kb(2)=Bc%jb; ke(2)=Bc%je; kb(3)=Bc%kb; ke(3)=Bc%ke
       k=mod(Bc%face-1,3)+1
	   if(Bc%face .le. 3) then          ! i-, j- or k-
	    kb(k)=kb(k)-1                  !  i=0 (Ghost点)
	   else
	    kb(k)=kb(k)+1                  ! i=nx+1 (Ghost点)
	   endif
	    ke(k)=kb(k)
       allocate(Ux_recv(3,kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))        ! 接收数组，按照目标数据的格式
       Num_data=3*(ke(1)-kb(1)+1)*(ke(2)-kb(2)+1)*(ke(3)-kb(3)+1)      ! 数据量
	   tag=B%Block_no*1000+Bc%f_no                                     ! tag 标记，标记块号+子面号
	   call MPI_Recv(Ux_recv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,Struct_Comm,status,ierr)

	   do k=kb(3),ke(3)
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
  		    B%x(i,j,k)=Ux_recv(1,i,j,k)    ! 临时数组拷贝入目标数据（Ghost 点）
            B%y(i,j,k)=Ux_recv(2,i,j,k)
            B%z(i,j,k)=Ux_recv(3,i,j,k)
           enddo
		 enddo
	   enddo
       deallocate(Ux_recv)
	
	enddo
	enddo
   end subroutine Coordinate_recv_mpi





!---------------------------------------------------------------------------------------
!  物理边界，虚网格的坐标采用外插法获得 ;  角部区域，采用内插获得   
	subroutine coordinate_boundary_and_corner(nMesh)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i1,j1,k1,i2,j2,k2,mBlock,ksub,nMesh
 ! ---------------------------------------------------------------------------------------- 
 do mBlock=1,Mesh(nMesh)%Num_Block

   B => Mesh(nMesh)%Block(mBlock)
   do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
 
   if(Bc%bc .ge. 0 ) then    ! 非内边界 Not inner boundary

!   对应非内边界，虚网格点上的坐标采用内点值外插获得 (x0=2*x1-x2)
    i1=Bc%ib; i2=Bc%ie; j1=Bc%jb; j2=Bc%je; k1=Bc%kb; k2=Bc%ke
    if(Bc%face .eq. 1 ) then             !  boundary  i-
      B%x(i1-1,j1:j2,k1:k2)  =2.d0*B%x(i1,j1:j2,k1:k2)-B%x(i1+1,j1:j2,k1:k2)
      B%y(i1-1,j1:j2,k1:k2)  =2.d0*B%y(i1,j1:j2,k1:k2)-B%y(i1+1,j1:j2,k1:k2)
      B%z(i1-1,j1:j2,k1:k2)  =2.d0*B%z(i1,j1:j2,k1:k2)-B%z(i1+1,j1:j2,k1:k2)
    else if(Bc%face .eq. 2 ) then        !  boundary  j-
      B%x(i1:i2,j1-1,k1:k2)  =2.d0*B%x(i1:i2,j1,k1:k2)-B%x(i1:i2,j1+1,k1:k2)
      B%y(i1:i2,j1-1,k1:k2)  =2.d0*B%y(i1:i2,j1,k1:k2)-B%y(i1:i2,j1+1,k1:k2)
      B%z(i1:i2,j1-1,k1:k2)  =2.d0*B%z(i1:i2,j1,k1:k2)-B%z(i1:i2,j1+1,k1:k2)
    else if (Bc%face .eq. 3 ) then       ! k-
      B%x(i1:i2,j1:j2,k1-1)  =2.d0*B%x(i1:i2,j1:j2,k1)-B%x(i1:i2,j1:j2,k1+1)
      B%y(i1:i2,j1:j2,k1-1)  =2.d0*B%y(i1:i2,j1:j2,k1)-B%y(i1:i2,j1:j2,k1+1)
      B%z(i1:i2,j1:j2,k1-1)  =2.d0*B%z(i1:i2,j1:j2,k1)-B%z(i1:i2,j1:j2,k1+1)
    else if (Bc%face .eq. 4 ) then       ! i+
      B%x(i2+1,j1:j2,k1:k2)  =2.d0*B%x(i2,j1:j2,k1:k2)-B%x(i2-1,j1:j2,k1:k2)
      B%y(i2+1,j1:j2,k1:k2)  =2.d0*B%y(i2,j1:j2,k1:k2)-B%y(i2-1,j1:j2,k1:k2)
      B%z(i2+1,j1:j2,k1:k2)  =2.d0*B%z(i2,j1:j2,k1:k2)-B%z(i2-1,j1:j2,k1:k2)
    else if (Bc%face .eq. 5 ) then       ! j+
      B%x(i1:i2,j2+1,k1:k2)  =2.d0*B%x(i1:i2,j2,k1:k2)-B%x(i1:i2,j2-1,k1:k2)
      B%y(i1:i2,j2+1,k1:k2)  =2.d0*B%y(i1:i2,j2,k1:k2)-B%y(i1:i2,j2-1,k1:k2)
      B%z(i1:i2,j2+1,k1:k2)  =2.d0*B%z(i1:i2,j2,k1:k2)-B%z(i1:i2,j2-1,k1:k2)
    else if (Bc%face .eq. 6 ) then       ! k+
      B%x(i1:i2,j1:j2,k2+1)  =2.d0*B%x(i1:i2,j1:j2,k2)-B%x(i1:i2,j1:j2,k2-1)
      B%y(i1:i2,j1:j2,k2+1)  =2.d0*B%y(i1:i2,j1:j2,k2)-B%y(i1:i2,j1:j2,k2-1)
      B%z(i1:i2,j1:j2,k2+1)  =2.d0*B%z(i1:i2,j1:j2,k2)-B%z(i1:i2,j1:j2,k2-1)
    endif
   endif
  enddo
  enddo

  !   角点坐标采用外插形式获得     

  do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     call get_xyz_conner(B%nx,B%ny,B%nz,B%x)   
     call get_xyz_conner(B%nx,B%ny,B%nz,B%y)
     call get_xyz_conner(B%nx,B%ny,B%nz,B%z)
  enddo
  end subroutine coordinate_boundary_and_corner





 !----------------------------------------------------------
 ! 计算角点(及角部区域)坐标  : 立方体的12个棱和8个顶点
    subroutine get_xyz_conner(nx,ny,nz,x)
    use precision_EC
    implicit none 
    integer:: nx,ny,nz
    real(PRE_EC),dimension(:,:,:),pointer::x
    
!    real(PRE_EC):: x(0:nx+1,0:ny+1,0:nz+1)
!   12条棱 （由相邻三个点外插得到）  

!-------------Modified By Li Xinliang ----------------------------------------------------------
     x(0,0,1:nz)=(x(1,0,1:nz)+x(0,1,1:nz))*0.5d0
     x(nx+1,0,1:nz)=(x(nx,0,1:nz)+x(nx+1,1,1:nz))*0.5d0
     x(0,ny+1,1:nz)=(x(1,ny+1,1:nz)+x(0,ny,1:nz))*0.5d0
     x(nx+1,ny+1,1:nz)=(x(nx,ny+1,1:nz)+x(nx+1,ny,1:nz))*0.5d0
   
     x(0,1:ny,0)=(x(1,1:ny,0)+x(0,1:ny,1))*0.5d0
     x(nx+1,1:ny,0)=(x(nx,1:ny,0)+x(nx+1,1:ny,1))*0.5d0
     x(0,1:ny,nz+1)=(x(1,1:ny,nz+1)+x(0,1:ny,nz))*0.5d0
     x(nx+1,1:ny,nz+1)=(x(nx,1:ny,nz+1)+x(nx+1,1:ny,nz))*0.5d0

     x(1:nx,0,0)=(x(1:nx,1,0)+x(1:nx,0,1))*0.5d0
     x(1:nx,ny+1,0)=(x(1:nx,ny,0)+x(1:nx,ny+1,1))*0.5d0
     x(1:nx,0,nz+1)=(x(1:nx,1,nz+1)+x(1:nx,0,nz))*0.5d0
     x(1:nx,ny+1,nz+1)=(x(1:nx,ny,nz+1)+x(1:nx,ny+1,nz))*0.5d0
!-------------------------------------------------------------------------
 ! 8个顶点 （由相邻6个点外插得到）
    x(0,0,0)=(2.d0*(x(1,0,0)+x(0,1,0)+x(0,0,1))-(x(1,1,0)+x(1,0,1)+x(0,1,1)))/3.d0
    x(nx+1,0,0)=(2.d0*(x(nx,0,0)+x(nx+1,1,0)+x(nx+1,0,1))-(x(nx,1,0)+x(nx,0,1)+x(nx+1,1,1)))/3.d0
    x(0,ny+1,0)=(2.d0*(x(1,ny+1,0)+x(0,ny,0)+x(0,ny+1,1))-(x(1,ny,0)+x(1,ny+1,1)+x(0,ny,1)))/3.d0
    x(nx+1,ny+1,0)=(2.d0*(x(nx,ny+1,0)+ x(nx+1,ny,0)+ x(nx+1,ny+1,1))-( x(nx,ny,0)+ x(nx,ny+1,1)+ x(nx+1,ny,1)))/3.d0
    x(0,0,nz+1)=(2.d0*(x(1,0,nz+1)+x(0,1,nz+1)+x(0,0,nz))-(x(1,1,nz+1)+x(1,0,nz)+x(0,1,nz)))/3.d0
    x(nx+1,0,nz+1)=(2.d0*(x(nx,0,nz+1)+x(nx+1,1,nz+1)+x(nx+1,0,nz))-(x(nx,1,nz+1)+x(nx,0,nz)+x(nx+1,1,nz)))/3.d0
    x(0,ny+1,nz+1)=(2.d0*(x(1,ny+1,nz+1)+x(0,ny,nz+1)+x(0,ny+1,nz))-(x(1,ny,nz+1)+x(1,ny+1,nz)+x(0,ny,nz)))/3.d0
    x(nx+1,ny+1,nz+1)=(2.d0*(x(nx,ny+1,nz+1)+ x(nx+1,ny,nz+1)+ x(nx+1,ny+1,nz))   &
                        -( x(nx,ny,nz+1)+ x(nx,ny+1,nz)+ x(nx+1,ny,nz)))/3.d0

     end
!-------------------------------------------------------------------
 ! 计算角点(及角部区域)的物理量  : 立方体的12个棱和8个顶点
    subroutine get_U_conner(nx,ny,nz,NVAR,U)
    use precision_EC
	implicit none 
    integer:: nx,ny,nz,NVAR
    real(PRE_EC),dimension(:,:,:,:),Pointer::U

! 12条棱  内插获得 
    U(:,0,0,1:nz-1)=(U(:,1,0,1:nz-1)+U(:,0,1,1:nz-1))*0.5d0
    U(:,nx,0,1:nz-1)=(U(:,nx-1,0,1:nz-1)+U(:,nx,1,1:nz-1))*0.5d0
    U(:,0,ny,1:nz-1)=(U(:,1,ny,1:nz-1)+U(:,0,ny-1,1:nz-1))*0.5d0
    U(:,nx,ny,1:nz-1)=(U(:,nx-1,ny,1:nz-1)+U(:,nx,ny-1,1:nz-1))*0.5d0
   
    U(:,0,1:ny-1,0)=(U(:,1,1:ny-1,0)+U(:,0,1:ny-1,1))*0.5d0
    U(:,nx,1:ny-1,0)=(U(:,nx-1,1:ny-1,0)+U(:,nx,1:ny-1,1))*0.5d0
    U(:,0,1:ny-1,nz)=(U(:,1,1:ny-1,nz)+U(:,0,1:ny-1,nz-1))*0.5d0
    U(:,nx,1:ny-1,nz)=(U(:,nx-1,1:ny-1,nz)+U(:,nx,1:ny-1,nz-1))*0.5d0

    U(:,1:nx-1,0,0)=(U(:,1:nx-1,1,0)+U(:,1:nx-1,0,1))*0.5d0
    U(:,1:nx-1,ny,0)=(U(:,1:nx-1,ny-1,0)+U(:,1:nx-1,ny,1))*0.5d0
    U(:,1:nx-1,0,nz)=(U(:,1:nx-1,1,nz)+U(:,1:nx-1,0,nz-1))*0.5d0
    U(:,1:nx-1,ny,nz)=(U(:,1:nx-1,ny-1,nz)+U(:,1:nx-1,ny,nz-1))*0.5d0


! 8个顶点
    U(:,0,0,0)=(2.d0*(U(:,1,0,0)+U(:,0,1,0)+U(:,0,0,1))-(U(:,1,1,0)+U(:,1,0,1)+U(:,0,1,1)))/3.d0
    U(:,nx,0,0)=(2.d0*(U(:,nx-1,0,0)+U(:,nx,1,0)+U(:,nx,0,1))-(U(:,nx-1,1,0)+U(:,nx-1,0,1)+U(:,nx,1,1)))/3.d0
    U(:,0,ny,0)=(2.d0*(U(:,1,ny,0)+U(:,0,ny-1,0)+U(:,0,ny,1))-(U(:,1,ny-1,0)+U(:,1,ny,1)+U(:,0,ny-1,1)))/3.d0
    U(:,nx,ny,0)=(2.d0*(U(:,nx-1,ny,0)+ U(:,nx,ny-1,0)+ U(:,nx,ny,1))-( U(:,nx-1,ny-1,0)+ U(:,nx-1,ny,1)+ U(:,nx,ny-1,1)))/3.d0
    U(:,0,0,nz)=(2.d0*(U(:,1,0,nz)+U(:,0,1,nz)+U(:,0,0,nz-1))-(U(:,1,1,nz)+U(:,1,0,nz-1)+U(:,0,1,nz-1)))/3.d0
    U(:,nx,0,nz)=(2.d0*(U(:,nx-1,0,nz)+U(:,nx,1,nz)+U(:,nx,0,nz-1))-(U(:,nx-1,1,nz)+U(:,nx-1,0,nz-1)+U(:,nx,1,nz-1)))/3.d0
    U(:,0,ny,nz)=(2.d0*(U(:,1,ny,nz)+U(:,0,ny-1,nz)+U(:,0,ny,nz-1))-(U(:,1,ny-1,nz)+U(:,1,ny,nz-1)+U(:,0,ny-1,nz-1)))/3.d0
    U(:,nx,ny,nz)=(2.d0*(U(:,nx-1,ny,nz)+ U(:,nx,ny-1,nz)+ U(:,nx,ny,nz-1))   & 
                   -( U(:,nx-1,ny-1,nz)+ U(:,nx-1,ny,nz-1)+ U(:,nx,ny-1,nz-1)))/3.d0

     end


! 接收（坐标信息）
    subroutine Coordinate_Periodic(nMesh)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,mBlock,ksub,nMesh
     integer:: kb(3),ke(3)
    real(PRE_EC):: Xp,Yp,Zp

 ! ---------------------------------------------------------------------------------------- 
  do mBlock=1,Mesh(nMesh)%Num_Block

   B => Mesh(nMesh)%Block(mBlock)
   do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)

!--------------------------------------------------------------------------------------------  
    if(Bc%bc .ne. BC_PeriodicL .and. Bc%bc .ne. BC_PeriodicR ) cycle               ! 非周期性边界
    if(Bc%bc == Bc_PeriodicL ) then    ! 左周期边界
	    Xp = - Periodic_dX;  Yp=- Periodic_dY; Zp= - Periodic_dZ
	endif

	if(Bc%bc == BC_PeriodicR)  then
      Xp =  Periodic_dX;  Yp= Periodic_dY; Zp=  Periodic_dZ
	endif

!      目标数据 （边界外的1层 Ghost）
       kb(1)=Bc%ib; ke(1)=Bc%ie; kb(2)=Bc%jb; ke(2)=Bc%je; kb(3)=Bc%kb; ke(3)=Bc%ke
       k=mod(Bc%face-1,3)+1
	   if(Bc%face .le. 3) then          ! i-, j- or k-
	    kb(k)=kb(k)-1                  !  i=0 (Ghost点)
	   else
	    kb(k)=kb(k)+1                  ! i=nx+1 (Ghost点)
	   endif
	    ke(k)=kb(k)

        do k=kb(3),ke(3)             
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             B%x(i,j,k)=B%x(i,j,k)+Xp    ! 根据周期条件， 添加一个增量
             B%y(i,j,k)=B%y(i,j,k)+Yp
             B%z(i,j,k)=B%z(i,j,k)+Zp
           enddo
		  enddo
		enddo


	
	enddo
	enddo
   end 


!----多块网格之间的连接--------------------------------------------------------------------------
! 交换LAP层 网格 (格心点) 信息


!---------Continue boundary (inner boundary) ------------------------------------------------------
!  传送相邻面的中心坐标给本块的虚网格  （便于使用周期性边界条件）； 
!  由于坐标存储在节点，物理量存储在中心点， 因而传输坐标与传输物理量的方法（下标对应方式）略有区别
!  MPI并行版； 
!  策略： 如两块(block)位于同一进程(proc),则采用直接传递(不通过MPI),以提高效率；
!  使用MPI传递时，首先使用MPI_Bsend()发送全部信息； 然后使用MPI_recv()接收。 

     subroutine update_Mesh_Center(nMesh)
     use Global_Var
     implicit none
     integer:: nMesh,ierr
! --------------------------------------------------------------------------------------- 
! 模块边界通信：  MPI 版本
    call Mesh_send_mpi(nMesh)   ! 使用MPI发送全部信息 ；如目标块也在本进程内，则不通过MPI,直接交换信息.
    call Mesh_recv_mpi(nMesh)   ! 接收全部信息
    call MPI_Barrier(Struct_Comm,ierr)
    call Mesh_Center_Periodic(nMesh)     ! 利用周期条件，修改坐标 （旋转一定角度, 或增加Dx,Dy,Dz）

  end subroutine update_Mesh_Center

!---------------------------------------------------------------

! 同一进程内两块之间的通信 （可直接通信）

     subroutine Mesh_send_mpi(nMesh)  
     use Global_Var
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,i1,j1,k1,mb,mBlock,ksub,nMesh,Send_to_ID,Num_Data,tag,ierr
     integer:: kb(3),ke(3),kb1(3),ke1(3),ka(3),L1,L2,L3,P(3),Ks(3)
     real(PRE_EC),allocatable,dimension(:,:,:,:):: Usend

 ! --------------------------------------------------------------------------------------- 
 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)

  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)

!-----------------------------------------------------------------------------------------------
      if(Bc%bc .ge. 0 ) cycle          ! 非内边界 inner boundary
      
!--------------------------------------------------------
!  B块为源数据； B1块为目标数据；   数据从源数据写入目标数据

!      源数据:临近边界的LAP层 内点;   i=kb(1) to ke(1), j=kb(2) to ke(2),  k= kb(3) to ke(3)
!      目标数据 （临近边界的LAP层 Ghost点,  现版本LAP=2）

!      源数据 
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                 ! k=1,2,3 为i,j,k方向
	   if(Bc%face .gt. 3) kb(k)=kb(k)-LAP   ! i+, j+ or k+ 面
       ke(k)=kb(k)+LAP-1

!      目标数据 （边界外扩充的LAP层Ghost点）
       kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
       k=mod(Bc%face1-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face1 .le. 3)  kb1(k)=kb1(k)-LAP   ! i-, j- or k- 面
       ke1(k)=kb1(k)+LAP-1


!  L1,L2,L3 : 描述维的连接 ; P(1),P(2),P(3): 描述连接次序(顺序or 逆序)
       L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)    ! L1=1 意味着（本维）与目标数据的第1维连接， =2 为与第2维连接, ...
       L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
       L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)

!   ks(k) : 目标数据第k维的起始下标
	   do k=1,3  !   目标数据 起始下标     
	    if(P(k) .gt. 0) then
		 ks(k)=kb(k)     ! 顺序，从kb开始
		else
		 ks(k)=ke(k)     ! 拟序， 从ke开始
	    endif
	  enddo
!
!-------------------------------------------------------------
      allocate(Usend(3,kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
       do k=kb1(3),ke1(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
		     ka(1)=i-kb1(1)            
			 ka(2)=j-kb1(2)
			 ka(3)=k-kb1(3)
			 i1=ks(1)+ka(L1)*P(1)        ! 源数据的下标 (i1,j1,k1), 从kb1到ke1 (需考虑 a.顺序或逆序, b. 维的连接), L1,L2,L3控制维的连接，P(k)控制顺序/逆序
			 j1=ks(2)+ka(L2)*P(2)
			 k1=ks(3)+ka(L3)*P(3)
			 Usend(1,i,j,k)=B%xc(i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
			 Usend(2,i,j,k)=B%yc(i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
			 Usend(3,i,j,k)=B%zc(i1,j1,k1)   ! 从源数据（内点） 拷贝入 目标数据 (临时数组)   
           
		   enddo
		 enddo
		enddo
  
  !    
	 Send_to_ID=B_proc(Bc%nb1)              ! 发送目标块所在的进程号

	 if( Send_to_ID .ne. my_id) then        ! 目标块不在本进程内
       Num_data=3*(ke1(1)-kb1(1)+1)*(ke1(2)-kb1(2)+1)*(ke1(3)-kb1(3)+1)      ! 数据量
       tag=  Bc%nb1*1000+Bc%f_no1                                         ! 标记; 块号Bc%nb1, 子面号Bc%f_no1 （向一个进程发送多个数据包时，用于识别）
	   call MPI_Bsend(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
!	   call MPI_send(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, Struct_Comm,ierr )
	 else 
!        目标块也在本进程内，直接写入 (不通过MPI)
         mb=B_n(Bc%nb1)          ! 块内部编号
		 B1 =>  Mesh(nMesh)%Block(mb)  ! 相邻块 （目标块）
	     do k=kb1(3),ke1(3)             ! 
	      do j=kb1(2),ke1(2)
		   do i=kb1(1),ke1(1)
			  B1%xc(i,j,k)=Usend(1,i,j,k)    
			  B1%yc(i,j,k)=Usend(2,i,j,k)    
			  B1%zc(i,j,k)=Usend(3,i,j,k)    
            enddo
		  enddo
		 enddo
	  endif	
	  deallocate(Usend)

  enddo
  enddo
 
  end
!--------------------------------------------------------------------------

! 将全部信息接收

subroutine Mesh_recv_mpi(nMesh) ! 使用MPI发送全部信息
     use Global_Var
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,Recv_from_ID,kb(3),ke(3)
	 integer:: tag,Num_DATA,ierr,Status(MPI_Status_SIZE)
     real(PRE_EC),allocatable,dimension(:,:,:,:):: Urecv

!---------------------------------------------------------------------------------------
! 接收阶段，将全部信息接收
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ge. 0 ) cycle               ! 内边界 inner boundary 或周期性边界
   Recv_from_ID=B_proc(Bc%nb1)           ! 相邻块（接收源块）所在的进程号
   if(Recv_from_ID .eq. my_id) cycle     ! 源块在本进程内，不使用MPI通信 (Umessage_send_mpi()已完成写入操作)  
   
!------------------------------------------------------------------------------------------------
!      目标数据 （边界外扩充的LAP层Ghost点）
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3)  kb(k)=kb(k)-LAP   ! i-, j- or k- 面
       ke(k)=kb(k)+LAP-1

       allocate(Urecv(3,kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))    ! 接收数组，按照目标数据的格式
       Num_data=3*(ke(1)-kb(1)+1)*(ke(2)-kb(2)+1)*(ke(3)-kb(3)+1)      ! 数据量
	   tag=B%Block_no*1000+Bc%f_no                                  ! tag 标记，标记块号+子面号
	   
	   call MPI_Recv(Urecv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,Struct_Comm,status,ierr)

        do k=kb(3),ke(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             B%xc(i,j,k)=Urecv(1,i,j,k)
             B%yc(i,j,k)=Urecv(2,i,j,k)
             B%zc(i,j,k)=Urecv(3,i,j,k)
           enddo
		  enddo
		enddo
       deallocate(Urecv)
   enddo
   enddo
  end


!------------------------------------------------------------------------
! 根据周期性条件， 对Ghost 点的坐标进行调整
! Left (BC_PeriodicL) ghost points shifted by -dX/dY/dZ; right by +dX/dY/dZ
 
  subroutine Mesh_Center_Periodic(nMesh)
     use Global_Var
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,kb(3),ke(3)
     real(PRE_EC):: Xp,Yp,Zp
!---------------------------------------------------------------------------------------
! 
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ne. BC_PeriodicL .and. Bc%bc .ne. BC_PeriodicR ) cycle               ! 非周期性边界
    if(Bc%bc == Bc_PeriodicL ) then    ! 左周期边界
	    Xp = - Periodic_dX;  Yp=- Periodic_dY; Zp= - Periodic_dZ
	endif

	if(Bc%bc == BC_PeriodicR)  then
      Xp =  Periodic_dX;  Yp= Periodic_dY; Zp=  Periodic_dZ
	endif

!------------------------------------------------------------------------------------------------
!    边界外扩充的LAP层Ghost点
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3)  kb(k)=kb(k)-LAP   ! i-, j- or k- 面
       ke(k)=kb(k)+LAP-1
	   
        do k=kb(3),ke(3)             
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             B%xc(i,j,k)=B%xc(i,j,k)+Xp
             B%yc(i,j,k)=B%yc(i,j,k)+Yp
             B%zc(i,j,k)=B%zc(i,j,k)+Zp
           enddo
		  enddo
		enddo


   enddo
   enddo
  end

  end module mod_struct_mpi
