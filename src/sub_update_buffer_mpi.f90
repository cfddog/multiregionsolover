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
     use interface_defines
     implicit none
     integer:: nMesh,ierr
! --------------------------------------------------------------------------------------- 
! 模块边界通信：  MPI 版本
    call Umessage_send_mpi(nMesh)   ! 使用MPI发送全部信息 ；如目标块也在本进程内，则不通过MPI,直接交换信息.
    call Umessage_recv_mpi(nMesh)   ! 接收全部信息
   
!   AC 逐块伪时间步循环内不做全局屏障：各进程的 AC 块数可能不等（np=2 的 1+2），
!   同一外层步内本子程序的调用次数就不同，MPI_Barrier 计数不匹配会直接死锁（实测
!   np=2 挂死）。U/p 鬼点交换全是成对的 MPI_Bsend/MPI_Recv，不依赖屏障定序。
    if(AC_In_Pseudo_Loop .eq. 0) call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
	if(If_TurboMachinary ==1) then
        call Umessage_Turbo_Periodic(nMesh)   !  处理周期边界条件 
    endif
    
!	call Umessage_corner(nMesh)     ! 角区数据信息（采用插值的方法） [DISABLED: corner treatment causes instability at T-junction]

  end subroutine update_buffer_onemesh

!---------------------------------------------------------------

! 同一进程内两块之间的通信 （可直接通信）

     subroutine Umessage_send_mpi(nMesh)  
     use Global_Var
     use interface_defines
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
      !  BC_Inner (-1) / Periodic (-2,-3) mark SAME-CLASS block interfaces: the
!  ghost buffer is exchanged for EVERY block type (compressible / low-speed
!  / porous / solid).  Cross-class interfaces (11..19) are physical
!  coupling faces -- they are handled by the couple_* routines and never
!  enter this buffer exchange (their code is >= 0 and is skipped above).
      
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
	   call MPI_Bsend(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
!	   call MPI_send(Usend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
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
!  B%p (low-speed / porous) is NOT exchanged here any more: p lives outside
!  U(1:5) and therefore cannot be packed into Usend.  Its halo is updated by
!  p_send_mpi/p_recv_mpi (update_p_buffer_onemesh), which handles the
!  same-process (direct write) and the cross-process (MPI_Bsend/Recv) cases in
!  ONE code path with ONE block-type filter.  The old fused copy that lived
!  here only covered same-process LOWSPEED<->LOWSPEED pairs, so the p halo
!  coverage depended on who called update_buffer_onemesh and the porous pairs
!  were silently skipped (see docs/程序说明 changelog).
  deallocate(Usend)

  enddo
  enddo
 
  end
!--------------------------------------------------------------------------

! 将全部信息接收

subroutine Umessage_recv_mpi(nMesh) ! 使用MPI发送全部信息
     use Global_Var
     use interface_defines
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
      !  BC_Inner (-1) / Periodic (-2,-3) mark SAME-CLASS block interfaces: the
!  ghost buffer is exchanged for EVERY block type (compressible / low-speed
!  / porous / solid).  Cross-class interfaces (11..19) are physical
!  coupling faces -- they are handled by the couple_* routines and never
!  enter this buffer exchange (their code is >= 0 and is skipped above).
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
	   
	   call MPI_Recv(Urecv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,MPI_COMM_WORLD,status,ierr)

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
   use interface_defines
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
     use interface_defines
     implicit none
     integer:: nMesh,ierr
	 call Coordinate_send_mpi(nMesh)
     call Coordinate_recv_mpi(nMesh)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     call Coordinate_Periodic(nMesh)              ! 周期边界的网格
     call coordinate_boundary_and_corner(nMesh)   ! 非内边界的网格
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

 !---test----------------------------------------------------------
   end 





    subroutine Coordinate_send_mpi(nMesh)
     use Global_Var
     use interface_defines
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
	   call MPI_Bsend(Ux_send,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
!	   call MPI_send(Ux_send,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
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
     use interface_defines
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
	   call MPI_Recv(Ux_recv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,MPI_COMM_WORLD,status,ierr)

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
     use interface_defines
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


!------------------------------------------------------------------------
! 根据周期性条件， 对Ghost 点的物理量调整  （径向、周向和轴向速度具有周期性，而直角坐标系下的速度分量并不具有周期性）

! 左侧(BC_PeriodicL)的点 -Turbo_Seta;  右侧的点 + Turbo_Seta
 
 subroutine Umessage_Turbo_Periodic(nMesh)
     use Global_Var
     use interface_defines
     implicit none
!---------------------------------------------    
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,Recv_from_ID,kb(3),ke(3)
     real(PRE_EC):: Seta,SetaP,seta0,ur,us,rr
!---------------------------------------------------------------------------------------
! 
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
   if(Bc%bc .ne. BC_PeriodicL .and. Bc%bc .ne. BC_PeriodicR ) cycle               ! 非周期性边界
     if(Bc%bc == Bc_PeriodicL ) then
	    setaP= -Turbo_Periodic_Seta    !  叶轮机周向计算域的周角
 	 else 
	    SetaP=  Turbo_Periodic_Seta
	 endif

!------------------------------------------------------------------------------------------------
!    边界外扩充的LAP层Ghost点
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1                    ! k=1,2,3 为i,j,k方向
	   if(Bc%face .le. 3)  kb(k)=kb(k)-LAP   ! i-, j- or k- 面
       ke(k)=kb(k)+LAP-1
   
!  径向速度ur,切向速度us按照新的坐标方向重建； 其余物理量不变
	   
        do k=kb(3),ke(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             rr=sqrt(B%yc(i,j,k)**2+B%zc(i,j,k)**2)
			 seta=acos(B%yc(i,j,k)/rr)
			 if(B%zc(i,j,k) < 0) seta=-seta
 			 seta0=seta-SetaP     ! 原角度
 			 ur=B%U(3,i,j,k)*cos(seta0)+B%U(4,i,j,k)*sin(seta0)          ! 径向速度 (乘以密度)
		     us=-B%U(3,i,j,k)*sin(seta0)+B%U(4,i,j,k)*cos(seta0)         ! 周向速度 (乘以密度)
		     
			 B%U(3,i,j,k)=ur*cos(seta)-us*sin(seta)                       ! 新值 
			 B%U(4,i,j,k)=ur*sin(seta)+us*cos(seta)
 		   enddo
		  enddo
		enddo

   enddo
   enddo
  end

! 接收（坐标信息）
    subroutine Coordinate_Periodic(nMesh)
     use Global_Var
     use interface_defines
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,mBlock,ksub,nMesh
     integer:: kb(3),ke(3)
     real(PRE_EC):: Seta,SetaP,rr,Xp,Yp,Zp

 ! ---------------------------------------------------------------------------------------- 
  do mBlock=1,Mesh(nMesh)%Num_Block

   B => Mesh(nMesh)%Block(mBlock)
   do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)

!--------------------------------------------------------------------------------------------  
    if(Bc%bc .ne. BC_PeriodicL .and. Bc%bc .ne. BC_PeriodicR ) cycle               ! 非周期性边界
    if(Bc%bc == Bc_PeriodicL ) then    ! 左周期边界
	    setaP= -Turbo_Periodic_Seta    !  叶轮机周向计算域的周角
	    Xp = - Periodic_dX;  Yp=- Periodic_dY; Zp= - Periodic_dZ
	endif

	if(Bc%bc == BC_PeriodicR)  then
	  SetaP=  Turbo_Periodic_Seta
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

     if( IF_TurboMachinary == 1) then  ! 叶轮机模式
        do k=kb(3),ke(3)             ! 目标数据的下标 (i,j,k) , 从kb到ke
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             rr=sqrt(B%y(i,j,k)**2+B%z(i,j,k)**2)
			 seta=acos(B%y(i,j,k)/rr)
			 if(B%z(i,j,k) < 0) seta=-seta
 			 seta=seta+SetaP     ! 根据周期性条件， 旋转 SetaP 角度
             B%y(i,j,k)=rr*cos(seta)
             B%z(i,j,k)=rr*sin(seta)
           enddo
		  enddo
		enddo
	 else             ! 非叶轮机模式
        do k=kb(3),ke(3)             
	     do j=kb(2),ke(2)
		   do i=kb(1),ke(1)
             B%x(i,j,k)=B%x(i,j,k)+Xp    ! 根据周期条件， 添加一个增量
             B%y(i,j,k)=B%y(i,j,k)+Yp
             B%z(i,j,k)=B%z(i,j,k)+Zp
           enddo
		  enddo
		enddo
     endif


	
	enddo
	enddo
   end 


!----------------------------------------------------------------------
! Ts buffer exchange: orchestrates Ts send/recv using bc_msg2
!----------------------------------------------------------------------
     subroutine update_Ts_buffer_onemesh(nMesh)
     use Global_Var
     use interface_defines
     implicit none
     integer:: nMesh,ierr
     call Ts_send_mpi(nMesh)
     call Ts_recv_mpi(nMesh)
!    no MPI_Barrier here: called inside the per-block solid GS (some ranks only)
   end subroutine update_Ts_buffer_onemesh

!----------------------------------------------------------------------
     subroutine Ts_send_mpi(nMesh)
     use Global_Var
     use const_var
     use interface_defines
     implicit none
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,i1,j1,k1,mb,mBlock,ksub,nMesh,Send_to_ID,Num_Data,tag,ierr
     integer:: kb(3),ke(3),kb1(3),ke1(3),ka(3),L1,L2,L3,P(3),Ks(3)
     real(PRE_EC),allocatable,dimension(:,:,:):: Tsend

 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)
   if(.not. associated(B%bc_msg2)) cycle
   do ksub=1,B%subface
     Bc=> B%bc_msg2(ksub)
     if(.not. is_interface_bc(Bc%bc)) cycle
     if(.not. (B%Block_type==BLOCK_SOLID .or. B%Block_type==BLOCK_POROUS)) cycle
!     Source data: LAP layers near boundary
     kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
     k=mod(Bc%face-1,3)+1
     if(Bc%face .gt. 3) kb(k)=kb(k)-LAP
     ke(k)=kb(k)+LAP-1
!     Target data: ghost layers
     kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
     k=mod(Bc%face1-1,3)+1
     if(Bc%face1 .le. 3) kb1(k)=kb1(k)-LAP
     ke1(k)=kb1(k)+LAP-1
     L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)
     L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
     L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)
     do k=1,3
       if(P(k) .gt. 0) then
         ks(k)=kb(k)
       else
         ks(k)=ke(k)
       endif
     enddo
     allocate(Tsend(kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
     do k=kb1(3),ke1(3)
       do j=kb1(2),ke1(2)
         do i=kb1(1),ke1(1)
           ka(1)=i-kb1(1)
           ka(2)=j-kb1(2)
           ka(3)=k-kb1(3)
           i1=ks(1)+ka(L1)*P(1)
           j1=ks(2)+ka(L2)*P(2)
           k1=ks(3)+ka(L3)*P(3)
           Tsend(i,j,k)=B%Ts(i1,j1,k1)
         enddo
       enddo
     enddo
     if(.not. (Block_Type_List(Bc%nb1)==BLOCK_SOLID .or. Block_Type_List(Bc%nb1)==BLOCK_POROUS)) cycle
     Send_to_ID=B_proc(Bc%nb1)
     if(Send_to_ID .ne. my_id) then
       Num_data=(ke1(1)-kb1(1)+1)*(ke1(2)-kb1(2)+1)*(ke1(3)-kb1(3)+1)
       tag=Bc%nb1*1000+Bc%f_no1
       call MPI_Bsend(Tsend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr)
     else
       mb=B_n(Bc%nb1)
       B1 => Mesh(nMesh)%Block(mb)
       do k=kb1(3),ke1(3)
         do j=kb1(2),ke1(2)
           do i=kb1(1),ke1(1)
             B1%Ts(i,j,k)=Tsend(i,j,k)
           enddo
         enddo
       enddo
     endif
     deallocate(Tsend)
   enddo
 enddo
 end subroutine Ts_send_mpi

!----------------------------------------------------------------------
 subroutine Ts_recv_mpi(nMesh)
     use Global_Var
     use const_var
     use interface_defines
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,mBlock,ksub,nMesh,Recv_from_ID,kb(3),ke(3)
     integer:: tag,Num_DATA,ierr,Status(MPI_Status_SIZE)
     real(PRE_EC),allocatable,dimension(:,:,:):: Trecv

 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)
   if(.not. associated(B%bc_msg2)) cycle
   do ksub=1,B%subface
     Bc=> B%bc_msg2(ksub)
     if(.not. is_interface_bc(Bc%bc)) cycle
     if(.not. (B%Block_type==BLOCK_SOLID .or. B%Block_type==BLOCK_POROUS)) cycle
     if(.not. (Block_Type_List(Bc%nb1)==BLOCK_SOLID .or. Block_Type_List(Bc%nb1)==BLOCK_POROUS)) cycle
     Recv_from_ID=B_proc(Bc%nb1)
     if(Recv_from_ID .eq. my_id) cycle
     kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
     k=mod(Bc%face-1,3)+1
     if(Bc%face .le. 3) kb(k)=kb(k)-LAP
     ke(k)=kb(k)+LAP-1
     allocate(Trecv(kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))
     Num_data=(ke(1)-kb(1)+1)*(ke(2)-kb(2)+1)*(ke(3)-kb(3)+1)
     tag=B%Block_no*1000+Bc%f_no
     call MPI_Recv(Trecv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,MPI_COMM_WORLD,status,ierr)
     do k=kb(3),ke(3)
       do j=kb(2),ke(2)
         do i=kb(1),ke(1)
           B%Ts(i,j,k)=Trecv(i,j,k)
         enddo
       enddo
     enddo
     deallocate(Trecv)
   enddo
 enddo
 end subroutine Ts_recv_mpi

!----------------------------------------------------------------------
! Same-class interface pressure halo for the AC low-speed / porous blocks
!
!   BLOCK_LOWSPEED / BLOCK_POROUS keep the primitive state in
!   U(1:4) = (rho,u,v,w) and the pressure OUTSIDE U in B%p.  Umessage_send_mpi
!   refreshes the U ghost layers of every same-class interface and already
!   copies B%p on the SAME-PROCESS path (see its `Send_to_ID .eq. my_id`
!   branch), which is why np=1 has always been consistent.  Across processes
!   the p ghost layer was never sent: the AC flux reconstructs XW(1)=p/rho
!   from the ghost cells, so a remote interface froze the interface pressure
!   at its initial level.
!
!   These two routines cover exactly that missing path and use the same
!   subface list, index orientation (L1..L3 / P) and tag scheme
!   (nb1*1000+f_no1) as Umessage_send/recv_mpi, so a pair of processes
!   matches link by link.  Only low-speed/porous pairs are exchanged here;
!   cross-class coupling faces (11..19) belong to the couple_* routines.
!
!   Called from the AC pseudo loop right after update_buffer_onemesh, which
!   means both sides must iterate that loop the same number of times -- that
!   is guaranteed by AC_MPI_Sync=1 (fixed AC_Max_Iter steps on every rank).
!----------------------------------------------------------------------
     subroutine update_p_buffer_onemesh(nMesh)
     use Global_Var
     implicit none
     integer:: nMesh
     call p_send_mpi(nMesh)
     call p_recv_mpi(nMesh)
!    no MPI_Barrier: called inside the per-block AC pseudo loop
     end subroutine update_p_buffer_onemesh

!----------------------------------------------------------------------
     subroutine p_send_mpi(nMesh)
     use Global_Var
     use const_var
     implicit none
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,m,i1,j1,k1,mb,mBlock,ksub,nMesh,Send_to_ID,Num_data,tag,ierr
     integer:: kb(3),ke(3),kb1(3),ke1(3),ka(3),L1,L2,L3,P(3),Ks(3)
     real(PRE_EC),allocatable,dimension(:,:,:):: psend

 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)
   if(.not. (B%Block_type==BLOCK_LOWSPEED .or. B%Block_type==BLOCK_POROUS)) cycle
   if(.not. associated(B%p)) cycle      ! 非 AC 场（未分配 p）：本子程序现在每个外层步都被调用
   do ksub=1,B%subface
     Bc=> B%bc_msg(ksub)
     if(Bc%bc .ge. 0) cycle                ! same-class (BC_Inner/periodic) only
     if(.not. (Block_Type_List(Bc%nb1)==BLOCK_LOWSPEED .or. &
               Block_Type_List(Bc%nb1)==BLOCK_POROUS)) cycle
!    Source data: LAP interior layers next to the boundary
     kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
     k=mod(Bc%face-1,3)+1
     if(Bc%face .gt. 3) kb(k)=kb(k)-LAP
     ke(k)=kb(k)+LAP-1
!    Target data: ghost layers
     kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
     k=mod(Bc%face1-1,3)+1
     if(Bc%face1 .le. 3) kb1(k)=kb1(k)-LAP
     ke1(k)=kb1(k)+LAP-1
     L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)
     L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
     L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)
     do k=1,3
       if(P(k) .gt. 0) then
         ks(k)=kb(k)
       else
         ks(k)=ke(k)
       endif
     enddo
     Send_to_ID=B_proc(Bc%nb1)
     if(Send_to_ID .eq. my_id) then
!      Target block is local -> write the ghost layers directly (no temporary
!      buffer, no MPI).  This is the ONLY same-process path for p: it used to
!      be duplicated/fused inside Umessage_send_mpi, where the filter was
!      LOWSPEED<->LOWSPEED only, so porous pairs never got a p halo and the
!      refresh cadence followed the U-exchange call sites.
       mb=B_n(Bc%nb1)
       B1 => Mesh(nMesh)%Block(mb)
       if(associated(B1%p)) then
         do k=kb1(3),ke1(3)
         do j=kb1(2),ke1(2)
         do i=kb1(1),ke1(1)
           ka(1)=i-kb1(1); ka(2)=j-kb1(2); ka(3)=k-kb1(3)
           i1=ks(1)+ka(L1)*P(1)
           j1=ks(2)+ka(L2)*P(2)
           k1=ks(3)+ka(L3)*P(3)
           B1%p(i,j,k)=B%p(i1,j1,k1)
         enddo; enddo; enddo
       endif
     else
       allocate(psend(kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
       do k=kb1(3),ke1(3)
         do j=kb1(2),ke1(2)
           do i=kb1(1),ke1(1)
             ka(1)=i-kb1(1); ka(2)=j-kb1(2); ka(3)=k-kb1(3)
             i1=ks(1)+ka(L1)*P(1)
             j1=ks(2)+ka(L2)*P(2)
             k1=ks(3)+ka(L3)*P(3)
             psend(i,j,k)=B%p(i1,j1,k1)
           enddo
         enddo
       enddo
       Num_data=(ke1(1)-kb1(1)+1)*(ke1(2)-kb1(2)+1)*(ke1(3)-kb1(3)+1)
       tag=Bc%nb1*1000+Bc%f_no1
       call MPI_Bsend(psend,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr)
       deallocate(psend)
     endif
   enddo
 enddo
 end subroutine p_send_mpi

!----------------------------------------------------------------------
 subroutine p_recv_mpi(nMesh)
     use Global_Var
     use const_var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: i,j,k,mBlock,ksub,nMesh,Recv_from_ID,kb(3),ke(3)
     integer:: tag,Num_data,ierr,Status(MPI_Status_SIZE)
     real(PRE_EC),allocatable,dimension(:,:,:):: precv

 do mBlock=1,Mesh(nMesh)%Num_Block
   B => Mesh(nMesh)%Block(mBlock)
   if(.not. (B%Block_type==BLOCK_LOWSPEED .or. B%Block_type==BLOCK_POROUS)) cycle
   if(.not. associated(B%p)) cycle      ! 非 AC 场（未分配 p）
   do ksub=1,B%subface
     Bc=> B%bc_msg(ksub)
     if(Bc%bc .ge. 0) cycle
     if(.not. (Block_Type_List(Bc%nb1)==BLOCK_LOWSPEED .or. &
               Block_Type_List(Bc%nb1)==BLOCK_POROUS)) cycle
     Recv_from_ID=B_proc(Bc%nb1)
     if(Recv_from_ID .eq. my_id) cycle     ! same process: p_send_mpi already copied
!    own ghost layers
     kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
     k=mod(Bc%face-1,3)+1
     if(Bc%face .le. 3) kb(k)=kb(k)-LAP
     ke(k)=kb(k)+LAP-1
     allocate(precv(kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))
     Num_data=(ke(1)-kb(1)+1)*(ke(2)-kb(2)+1)*(ke(3)-kb(3)+1)
     tag=B%Block_no*1000+Bc%f_no
     call MPI_Recv(precv,Num_data,OCFD_DATA_TYPE,Recv_from_ID,tag,MPI_COMM_WORLD,status,ierr)
     do k=kb(3),ke(3)
       do j=kb(2),ke(2)
         do i=kb(1),ke(1)
           B%p(i,j,k)=precv(i,j,k)
         enddo
       enddo
     enddo
     deallocate(precv)
   enddo
 enddo
 end subroutine p_recv_mpi


!===============================================================================
! Halo-exchange self test -- called from main when IF_Debug = 3.
!
! Verifies that the same-class ghost buffer exchange REALLY transfers data, for
! every vector component of U(1:NVAR) and for the low-speed/porous pressure p,
! on every subface of every local block.
!
! Method: overwrite the whole array (interior + ghost) with a value that is
! unique per (global block, i, j, k, component)
!
!     sent = ((blk*1000 + i+16)*1000 + j+16)*1000 + (k+16) + m/8
!
! run update_buffer_onemesh + update_p_buffer_onemesh, and compare every ghost
! cell of target block nb1 with sent(blk_src, i1,j1,k1, m), where blk_src is
! the block that OWNS the source cell (B%Block_no) and (i1,j1,k1) is obtained
! with the SAME L1/L2/L3/P mapping that the packers use.  Both sides of an
! interface know (B, nb1, i1,j1,k1), so the expected value can be computed
! locally on each rank: no cross-rank memory access and no extra communication
! is needed, which is exactly what lets this test also cover the cross-process
! subfaces.
!
! The p check is applied only to LOWSPEED/POROUS pairs, i.e. exactly the
! subfaces that p_send_mpi/p_recv_mpi are supposed to serve.
!
! The original state is saved before the test and restored afterwards, so the
! initial condition is not perturbed.  Output: one line per subface (each rank,
! its own subfaces -- "local = T/F" shows whether the subface was exchanged in
! memory or through MPI) plus an MPI_Allreduce summary with the global maximum
! error and a PASS/FAIL verdict on rank 0.
!===============================================================================
 subroutine check_buffer_exchange(nMesh)
   use Global_Var
   use const_var
   implicit none
   Type (Block_TYPE),pointer:: B,B1
   Type (BC_MSG_TYPE),pointer:: Bc
   integer,intent(in):: nMesh
   integer:: mBlock,ksub,m,i,j,k,i1,j1,k1,mb
   integer:: kb(3),ke(3),kb1(3),ke1(3),ka(3),L1,L2,L3,P(3),ks(3)
   integer:: ncel,ncel_all,nface_l,nface_c,nface_l_all,nface_c_all
   integer:: nbad,nbad_all,nfcel,ntot,ntotp,iu,ip,ierr
   real(PRE_EC):: ex,uerr,perr,eu,ep,eu_all,ep_all
   real(PRE_EC),allocatable:: Usave(:),psave(:)
   real(PRE_EC),allocatable:: Uexp(:,:,:,:),pexp(:,:,:),Ugot(:,:,:,:),pgot(:,:,:)
   integer:: status(MPI_STATUS_SIZE)
   integer,parameter:: CK_TAG=900000      ! tag base of the verification messages
   logical:: chkp
   interface
     subroutine update_buffer_onemesh(nMesh)
       use precision_EC
       integer:: nMesh
     end subroutine update_buffer_onemesh
     subroutine update_p_buffer_onemesh(nMesh)
       use precision_EC
       integer:: nMesh
     end subroutine update_p_buffer_onemesh
   end interface

!---- pass 1: save the original state, then fill it with the sentinel --------
   ntot=0; ntotp=0
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%U)) cycle
     ntot=ntot+size(B%U)
     if(associated(B%p)) ntotp=ntotp+size(B%p)
   enddo
   allocate(Usave(max(ntot,1)),psave(max(ntotp,1)))
   iu=0; ip=0
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%U)) cycle
     do k=lbound(B%U,4),ubound(B%U,4)
     do j=lbound(B%U,3),ubound(B%U,3)
     do i=lbound(B%U,2),ubound(B%U,2)
       do m=1,NVAR
         iu=iu+1; Usave(iu)=B%U(m,i,j,k)
         B%U(m,i,j,k)=sent(B%Block_no,i,j,k,m)
       enddo
     enddo; enddo; enddo
     if(associated(B%p)) then
       do k=lbound(B%p,3),ubound(B%p,3)
       do j=lbound(B%p,2),ubound(B%p,2)
       do i=lbound(B%p,1),ubound(B%p,1)
         ip=ip+1; psave(ip)=B%p(i,j,k)
         B%p(i,j,k)=sent(B%Block_no,i,j,k,0)
       enddo; enddo; enddo
     endif
   enddo

!---- the exchange under test ------------------------------------------------
   call update_buffer_onemesh(nMesh)
   call update_p_buffer_onemesh(nMesh)


!---- pass 2: verify every SAME-PROCESS subface directly --------------------
!  (for a cross-process subface the neighbour block does not exist on this rank,
!   so B1 cannot be dereferenced: those are verified in passes 3-4)
   eu=0.d0; ep=0.d0; ncel=0; nbad=0; nface_l=0; nface_c=0
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%U)) cycle
     do ksub=1,B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .ge. 0) cycle                 ! physical BC, not a halo face
       if(B_proc(Bc%nb1) .ne. my_id) cycle    ! cross-process: passes 3-4
!      source range (identical to Umessage_send_mpi / p_send_mpi)
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1
       if(Bc%face .gt. 3) kb(k)=kb(k)-LAP
       ke(k)=kb(k)+LAP-1
!      target (ghost) range
       kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
       k=mod(Bc%face1-1,3)+1
       if(Bc%face1 .le. 3) kb1(k)=kb1(k)-LAP
       ke1(k)=kb1(k)+LAP-1
!      dimension mapping / order
       L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)
       L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
       L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)
       do k=1,3
         if(P(k) .gt. 0) then
           ks(k)=kb(k)
         else
           ks(k)=ke(k)
         endif
       enddo
       mb=B_n(Bc%nb1)
       B1 => Mesh(nMesh)%Block(mb)
!      does this subface carry a pressure halo? (same filter as p_send_mpi)
       chkp = associated(B1%p) .and. &
              (B%Block_type==BLOCK_LOWSPEED .or. B%Block_type==BLOCK_POROUS) .and. &
              (Block_Type_List(Bc%nb1)==BLOCK_LOWSPEED .or. &
               Block_Type_List(Bc%nb1)==BLOCK_POROUS)
       uerr=0.d0; perr=0.d0; nfcel=0
       do k=kb1(3),ke1(3)
       do j=kb1(2),ke1(2)
       do i=kb1(1),ke1(1)
         ka(1)=i-kb1(1); ka(2)=j-kb1(2); ka(3)=k-kb1(3)
         i1=ks(1)+ka(L1)*P(1)
         j1=ks(2)+ka(L2)*P(2)
         k1=ks(3)+ka(L3)*P(3)
         nfcel=nfcel+1
         do m=1,NVAR
           ex=sent(B%Block_no,i1,j1,k1,m)
           uerr=max(uerr,abs(B1%U(m,i,j,k)-ex))
         enddo
         if(chkp) then
           ex=sent(B%Block_no,i1,j1,k1,0)
           perr=max(perr,abs(B1%p(i,j,k)-ex))
         endif
       enddo; enddo; enddo
       ncel=ncel+nfcel
       if(nfcel .gt. 0 .and. (uerr .gt. 0.d0 .or. perr .gt. 0.d0)) nbad=nbad+1
       nface_l=nface_l+1
       eu=max(eu,uerr); ep=max(ep,perr)
       write(*,'(A,I3,A,I4,A,I2,A,I4,A,I2,A,I6,A,ES11.3,A,ES11.3,A)') &
         ' [bufchk] rank',my_id,' blk',B%Block_no,' face',Bc%face, &
         ' <-> blk',Bc%nb1,' face',Bc%face1,' cells=',nfcel, &
         ' max|dU|=',uerr,' max|dp|=',perr,' local= T'
     enddo
   enddo



!---- pass 3: cross-process subfaces -- send the neighbour a buffer that holds
!     EXACTLY what its ghost layers should contain after the exchange ----------
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%U)) cycle
     do ksub=1,B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .ge. 0) cycle
       if(B_proc(Bc%nb1) .eq. my_id) cycle
!      source range (B) and target ghost range (neighbour), as in the packers
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1
       if(Bc%face .gt. 3) kb(k)=kb(k)-LAP
       ke(k)=kb(k)+LAP-1
       kb1(1)=Bc%ib1; ke1(1)=Bc%ie1-1; kb1(2)=Bc%jb1; ke1(2)=Bc%je1-1; kb1(3)=Bc%kb1; ke1(3)=Bc%ke1-1
       k=mod(Bc%face1-1,3)+1
       if(Bc%face1 .le. 3) kb1(k)=kb1(k)-LAP
       ke1(k)=kb1(k)+LAP-1
       L1=abs(Bc%L1); P(1)=sign(1,Bc%L1)
       L2=abs(Bc%L2); P(2)=sign(1,Bc%L2)
       L3=abs(Bc%L3); P(3)=sign(1,Bc%L3)
       do k=1,3
         if(P(k) .gt. 0) then
           ks(k)=kb(k)
         else
           ks(k)=ke(k)
         endif
       enddo
       chkp = (B%Block_type==BLOCK_LOWSPEED .or. B%Block_type==BLOCK_POROUS) .and. &
              (Block_Type_List(Bc%nb1)==BLOCK_LOWSPEED .or. &
               Block_Type_List(Bc%nb1)==BLOCK_POROUS)
       allocate(Uexp(NVAR,kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
       do k=kb1(3),ke1(3)
       do j=kb1(2),ke1(2)
       do i=kb1(1),ke1(1)
         ka(1)=i-kb1(1); ka(2)=j-kb1(2); ka(3)=k-kb1(3)
         i1=ks(1)+ka(L1)*P(1)
         j1=ks(2)+ka(L2)*P(2)
         k1=ks(3)+ka(L3)*P(3)
         do m=1,NVAR
           Uexp(m,i,j,k)=sent(B%Block_no,i1,j1,k1,m)
         enddo
       enddo; enddo; enddo
       call MPI_Bsend(Uexp,size(Uexp),OCFD_DATA_TYPE,B_proc(Bc%nb1), &
                      CK_TAG+Bc%nb1*1000+Bc%f_no1,MPI_COMM_WORLD,ierr)
       deallocate(Uexp)
       if(chkp) then
         allocate(pexp(kb1(1):ke1(1),kb1(2):ke1(2),kb1(3):ke1(3)))
         do k=kb1(3),ke1(3)
         do j=kb1(2),ke1(2)
         do i=kb1(1),ke1(1)
           ka(1)=i-kb1(1); ka(2)=j-kb1(2); ka(3)=k-kb1(3)
           i1=ks(1)+ka(L1)*P(1)
           j1=ks(2)+ka(L2)*P(2)
           k1=ks(3)+ka(L3)*P(3)
           pexp(i,j,k)=sent(B%Block_no,i1,j1,k1,0)
         enddo; enddo; enddo
         call MPI_Bsend(pexp,size(pexp),OCFD_DATA_TYPE,B_proc(Bc%nb1), &
                        CK_TAG+500000+Bc%nb1*1000+Bc%f_no1,MPI_COMM_WORLD,ierr)
         deallocate(pexp)
       endif
     enddo
   enddo

!---- pass 4: cross-process subfaces -- receive the expected buffer and
!     compare it with the ghost layers this rank actually holds ---------------
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%U)) cycle
     do ksub=1,B%subface
       Bc => B%bc_msg(ksub)
       if(Bc%bc .ge. 0) cycle
       if(B_proc(Bc%nb1) .eq. my_id) cycle
!      my own ghost range next to this face
       kb(1)=Bc%ib; ke(1)=Bc%ie-1; kb(2)=Bc%jb; ke(2)=Bc%je-1; kb(3)=Bc%kb; ke(3)=Bc%ke-1
       k=mod(Bc%face-1,3)+1
       if(Bc%face .le. 3) kb(k)=kb(k)-LAP
       ke(k)=kb(k)+LAP-1
       chkp = (B%Block_type==BLOCK_LOWSPEED .or. B%Block_type==BLOCK_POROUS) .and. &
              (Block_Type_List(Bc%nb1)==BLOCK_LOWSPEED .or. &
               Block_Type_List(Bc%nb1)==BLOCK_POROUS)
       allocate(Ugot(NVAR,kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))
       call MPI_Recv(Ugot,size(Ugot),OCFD_DATA_TYPE,B_proc(Bc%nb1), &
                     CK_TAG+B%Block_no*1000+Bc%f_no,MPI_COMM_WORLD,status,ierr)
       uerr=0.d0; perr=0.d0; nfcel=0
       do k=kb(3),ke(3)
       do j=kb(2),ke(2)
       do i=kb(1),ke(1)
         nfcel=nfcel+1
         do m=1,NVAR
           uerr=max(uerr,abs(B%U(m,i,j,k)-Ugot(m,i,j,k)))
         enddo
       enddo; enddo; enddo
       deallocate(Ugot)
       if(chkp) then
         allocate(pgot(kb(1):ke(1),kb(2):ke(2),kb(3):ke(3)))
         call MPI_Recv(pgot,size(pgot),OCFD_DATA_TYPE,B_proc(Bc%nb1), &
                       CK_TAG+500000+B%Block_no*1000+Bc%f_no,MPI_COMM_WORLD,status,ierr)
         do k=kb(3),ke(3)
         do j=kb(2),ke(2)
         do i=kb(1),ke(1)
           perr=max(perr,abs(B%p(i,j,k)-pgot(i,j,k)))
         enddo; enddo; enddo
         deallocate(pgot)
       endif
       ncel=ncel+nfcel
       if(nfcel .gt. 0 .and. (uerr .gt. 0.d0 .or. perr .gt. 0.d0)) nbad=nbad+1
       nface_c=nface_c+1
       eu=max(eu,uerr); ep=max(ep,perr)
       write(*,'(A,I3,A,I4,A,I2,A,I4,A,I2,A,I6,A,ES11.3,A,ES11.3,A)') &
         ' [bufchk] rank',my_id,' blk',B%Block_no,' face',Bc%face, &
         ' <- blk',Bc%nb1,' face',Bc%face1,' cells=',nfcel, &
         ' max|dU|=',uerr,' max|dp|=',perr,' local= F'
     enddo
   enddo

!---- pass 5: restore the original state ------------------------------------
   iu=0; ip=0
   do mBlock=1,Mesh(nMesh)%Num_Block
     B => Mesh(nMesh)%Block(mBlock)
     if(.not. associated(B%U)) cycle
     do k=lbound(B%U,4),ubound(B%U,4)
     do j=lbound(B%U,3),ubound(B%U,3)
     do i=lbound(B%U,2),ubound(B%U,2)
       do m=1,NVAR
         iu=iu+1; B%U(m,i,j,k)=Usave(iu)
       enddo
     enddo; enddo; enddo
     if(associated(B%p)) then
       do k=lbound(B%p,3),ubound(B%p,3)
       do j=lbound(B%p,2),ubound(B%p,2)
       do i=lbound(B%p,1),ubound(B%p,1)
         ip=ip+1; B%p(i,j,k)=psave(ip)
       enddo; enddo; enddo
     endif
   enddo
   deallocate(Usave,psave)

!---- global summary --------------------------------------------------------
   call MPI_Allreduce(eu,eu_all,1,OCFD_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(ep,ep_all,1,OCFD_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(ncel,ncel_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(nbad,nbad_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(nface_l,nface_l_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPI_Allreduce(nface_c,nface_c_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0) then
     write(*,'(A)') ' [bufchk] ===== halo-exchange self test (IF_Debug=3) ====='
     write(*,'(A,I4,A,I8)') ' [bufchk] np =',Total_proc,'  ghost cells checked =',ncel_all
     write(*,'(A,I6,A,I6)') ' [bufchk] subfaces: same-process =',nface_l_all, &
                            '  cross-process =',nface_c_all
     write(*,'(A,ES12.3,A,ES12.3)') ' [bufchk] global max|dU| =',eu_all, &
                                    '  global max|dp| =',ep_all
     if(nbad_all .eq. 0 .and. ncel_all .gt. 0) then
       write(*,'(A)') ' [bufchk] PASS: every ghost cell equals the mapped neighbour value'
     else
       write(*,'(A,I6)') ' [bufchk] FAIL: subfaces with a mismatch =',nbad_all
     endif
   endif

   contains
     real(PRE_EC) function sent(blk,ii,jj,kk,mm)
       integer,intent(in):: blk,ii,jj,kk,mm
       sent=( (dble(blk)*1000.d0+dble(ii+16))*1000.d0+dble(jj+16) )*1000.d0 &
            +dble(kk+16)+dble(mm)*0.125d0
     end function sent
 end subroutine check_buffer_exchange
