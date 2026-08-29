
!  ��ȡ������Ϣ
!  ������(0�Ž���)��ȡ���ݣ��������̽���
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
	if( Mesh_File_Format .eq. 1) then   ! ��ʽ�ļ�
     open(99,file="Mesh3d.dat")
     read(99,*) NB   ! Block number
     allocate( NI(NB),NJ(NB),NK(NB) )
     read(99,*) (NI(k), NJ(k), NK(k), k=1,NB)
	else                                ! �޸�ʽ�ļ�
     open(99,file="Mesh3d.dat",form="unformatted")
     read(99) NB                       ! �ܿ���
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
   if(B_proc(m) .eq. 0) then            ! ��Щ�����ڸ�����
      mt=B_n(m)                          ! �ÿ��ڽ����ڲ��ı��
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
     else                        ! ���ÿ����ݷ��ͳ�
	   Num_data=nx*ny*nz*3
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
!	   call MPI_Bsend(Ux,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
	   call MPI_send(Ux,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
     endif
     deallocate(Ux)
   enddo

   deallocate(NI,NJ,NK)
   close(99)
  
  else     ! �Ǹ��ڵ�
    do m=1,MP%Num_Block     ! �����̰����Ŀ�
     B=>MP%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
   	 allocate(Ux(nx,ny,nz,3))
   	 Num_data=nx*ny*nz*3
 	 tag=m
	 call MPI_Recv(Ux,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,Status,ierr )
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
  
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "read Mesh3d.dat OK"
 end subroutine read_main_mesh

!-------------------------------------------------------------------------------------
!  ��ȡ�������� ������ľ���
!  ������(0�Ž���)��ȡ���ݣ��������̽���
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

!------�����̶�ȡ����----------------------------- 
  if(my_id .eq. 0) then 
    print*, " read distance to the wall:  wall_dist.dat"
  
   open(99,file="wall_dist.dat",form="unformatted")

   do m=1,Total_block
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(dw(nx-1,ny-1,nz-1))
     read(99) (((dw(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1) 
     
	 if(B_proc(m) .eq. 0) then            ! ��Щ�����ڸ�����
      mt=B_n(m)                          ! �ÿ��ڽ����ڲ��ı��
	  B=>MP%Block(mt)
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	   B%dw(i,j,k)=dw(i,j,k)
	  enddo
	  enddo
 	  enddo
     else                        ! ���ÿ����ݷ��ͳ�
	   Num_data=(nx-1)*(ny-1)*(nz-1)
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
	   call MPI_send(dw,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
     endif
     deallocate(dw)
   enddo
   close(99)
  
  else     ! �Ǹ��ڵ�
    do m=1,MP%Num_Block     ! �����̰����Ŀ�
     B=>MP%Block(m)
	 nx=B%nx; ny=B%ny; nz=B%nz
   	 allocate(dw(nx-1,ny-1,nz-1))
	 Num_data=(nx-1)*(ny-1)*(nz-1)
 	 tag=m
	 call MPI_Recv(dw,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,Status,ierr )
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
  
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "read wall_dist.dat OK"
 end subroutine read_dw




!----------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!  ��ȡ����: d,u,v,w,T;  SA, SST �еı���
!  ������(0�Ž���)��ȡ���ݣ��������̽���
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

!------�����̶�ȡ���ݣ����͵���������----------------------------- 
  if(my_id .eq. 0) then 
    print*, " read flow data:  flow3d.dat"
 
     open(99,file="flow3d.dat",form="unformatted")
    
 	 if(NVAR1 .eq. 6) then        ! 6���Ա���
      Inquire(file="SA3d.dat",exist=Ex)
      if(Ex) then
        open(100,file="SA3d.dat",form="unformatted")
      endif
	 endif
     
	 if(NVAR1 .eq. 7) then          ! 7���Ա���
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
!------------����--------------------------------
     
	 if(B_proc(m) .eq. 0) then            ! ��Щ�����ڸ�����
      mt=B_n(m)                           ! �ÿ��ڽ����ڲ��ı��
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
     else                        ! ���ÿ����ݷ��ͳ�
	   Num_data=(nx+1)*(ny+1)*(nz+1)*NVAR1
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
	   call MPI_send(U,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
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

  else     ! �Ǹ��ڵ�
   
    do m=1,MP%Num_Block     ! �����̰����Ŀ�
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz,NVAR1))
	  Num_data=(nx+1)*(ny+1)*(nz+1)*NVAR1
 	  tag=m
	 
	 call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,Status,ierr )
      
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
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "read flow3d.dat OK"

!----------------------------------Transform data----------------
! ���������λd,u,v,w,T, ת��Ϊ�غ����
 do m=1,MP%Num_Block     ! �����̰����Ŀ�
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

!---------------����ʱ�䲽-----------------------
  if(my_id .eq. 0) then
   Inquire(file="Step_mess.dat",exist=Ex)
    if(Ex) then
     open(88,file="Step_mess.dat")
     read(88,*) Mesh(1)%Kstep, Mesh(1)%tt
    endif
  print*, "Init data OK, Kstep,tt=", Mesh(1)%Kstep, Mesh(1)%tt
 endif

call MPI_bcast(Mesh(1)%Kstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_bcast(Mesh(1)%tt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!---------------------------------------------------
 end subroutine read_flow_data








!----------------------------------------------------------------------
!  ������μ������� ��Plot3d��ʽ��, ��ϸ����flow3d.dat  

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


!   ������㲽����ʱ����Ϣ
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
     write(filename, "('flow3d-'I8.8'.dat')") Mesh(1)%Kstep
     open(99,file=filename,form="unformatted")    ! d,u,v,w,T
   endif




   if(NVAR1 .eq. 6) open(100,file="SA3d.dat",form="unformatted")    ! U6
   if(NVAR1 .eq. 7) open(101,file="SST3d.dat",form="unformatted")   ! U6,U7
  
   do m=1, Total_block   ! ȫ����
     
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz,NVAR1))

	if(B_proc(m) .eq. 0) then            ! ��Щ�����ڸ�����
      mt=B_n(m)                           ! �ÿ��ڽ����ڲ��ı��
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

    else                        ! ���ոÿ���Ϣ
	   Num_data=NVAR1*(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! �ڸÿ��еı��
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD,Status,ierr )
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
 
 else     ! ��0�ڵ�

    do m=1,MP%Num_Block     ! �����̰����Ŀ�
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
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "write flow3d.dat OK"

  end subroutine output_flow




!----------------------------------------------------------------------
!  �������ճ��ϵ��vt ��Plot3d��ʽ��, ��ϸ����vt.dat  

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
 
   do m=1, Total_block   ! ȫ����
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz))

	if(B_proc(m) .eq. 0) then            ! ��Щ�����ڸ�����
      mt=B_n(m)                           ! �ÿ��ڽ����ڲ��ı��
	  B=>MP%Block(mt)
	 
	  do k=0,nz
	  do j=0,ny
	  do i=0,nx
		 U(i,j,k)=B%mu_t(i,j,k)*Re                     ! mu_t
	  enddo
	  enddo
 	  enddo
    else                        ! ���ոÿ���Ϣ
	   Num_data=(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! �ڸÿ��еı��
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD,Status,ierr )
    endif
! write Data ....
    
	write(99) ((( U(i,j,k),i=0,nx),j=0,ny),k=0,nz)    

	deallocate(U)
   enddo
   close(99)
   close(100)
   close(101)
 
 else     ! ��0�ڵ�

    do m=1,MP%Num_Block     ! �����̰����Ŀ�
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
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "write vt.dat OK"

  end subroutine output_vt



!----------------------------------------------------------------------
!  ���������ľ���  

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
  
   do m=1, Total_block   ! ȫ����
 	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(nx-1,ny-1,nz-1))

	if(B_proc(m) .eq. 0) then            ! ��Щ�����ڸ�����
      mt=B_n(m)                          ! �ÿ��ڽ����ڲ��ı��
	  B=>MP%Block(mt)
	 
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
        U(i,j,k)=B%dw(i,j,k)
	  enddo
	  enddo
 	  enddo
    else                        ! ���ոÿ���Ϣ
	   Num_data=(nx-1)*(ny-1)*(nz-1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m)             ! �ڸÿ��еı��
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD,Status,ierr )
    endif
! write Data ....
 	write(99) ((( U(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)    
	deallocate(U)
   enddo
   close(99)
 
 else     ! ��0�ڵ�

    do m=1,MP%Num_Block     ! �����̰����Ŀ�
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
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "write wall_dist.dat OK"

  end subroutine write_dw




!----Boundary message (bc3d.inc, OpenCFD-EC Build-in format)----------------------------------------------------------
! OpenCFD-EC�ڽ���.inc�߽����Ӹ�ʽ���� Gridgen ��.inp��ʽ�����Ϸ�չ���ġ�
!  ��.inp��ʽ����һЩ������Ϣ��������������f_no,
!  ������face �Լ����ӵ������f_no1,���ӵ�������face1
!  �Լ����Ӵ���L1, L2, L3  (����L1=1��ʾ��ά�����ӿ�ĵ�1ά�����ӣ� L1=-1��ʾ�����ӿ�ĵ�1Ϊ��������).
!  ��Щ������ϢΪ��-��֮���ͨ�ţ�������MPI����ͨ�ţ��ṩ�˱����������ڼ�ͨ�Ŵ���
!--------------------------------------------------------------------------------------------------------------------
  subroutine read_inc 
   use Global_Var
   implicit none
   integer,parameter:: NC=21         ! .inc�ļ�ÿ��Ԫ��21��Ԫ��
   integer:: nx,ny,nz,NB,Nsub,m,mt,k,j
   integer:: Send_to_ID,tag,ierr,Status(MPI_Status_SIZE)
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer,pointer,dimension(:,:):: Bs 

 !  ��Gridgen .inp ��ʽת��Ϊ .inc��ʽ   
	if(my_id .eq. 0) then
	  call  convert_inp_inc 
	endif


!  read bc3d.inc, �����̶��룬����������������
   if(my_id .eq. 0) then
     print*, "read bc3d.inc (Link/boundary file)......"
     open(88,file="bc3d.inc")
     read(88,*)
     read(88,*) NB
 !   ��ȡ.inc�ļ��е�Ԫ��   
    do m=1,NB
     read(88,*) nx,ny,nz
     read(88,*)
     read(88,*) Nsub   ! m���������
 	 allocate(Bs(NC,Nsub))
	 do k=1,Nsub
	  read(88,*) (Bs(j,k),j=1,9)
	  read(88,*) (Bs(j,k),j=10,21)
!	 read(88,*) Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
!	 read(88,*) Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
     enddo
!     ���ÿ���Ϣ���ͳ�ȥ
      
     if(B_Proc(m) .eq. 0) then
!          �ÿ���0����
       	 mt=B_n(m)   ! ��0�����е��ڲ����
		 B=>Mesh(1)%Block(mt)
		 B%subface=Nsub   
	     allocate(B%bc_msg(B%subface))   ! �߽�����
         do k=1,Nsub
		 Bc=>B%bc_msg(k)
          Bc%ib=Bs(1,k); Bc%ie=Bs(2,k); Bc%jb=Bs(3,k); Bc%je=Bs(4,k)
		  Bc%kb=Bs(5,k); Bc%ke=Bs(6,k); Bc%bc=Bs(7,k); Bc%face=Bs(8,k); Bc%f_no=Bs(9,k)

          Bc%ib1=Bs(10,k); Bc%ie1=Bs(11,k); Bc%jb1=Bs(12,k); Bc%je1=Bs(13,k)
		  Bc%kb1=Bs(14,k); Bc%ke1=Bs(15,k); Bc%nb1=Bs(16,k); Bc%face1=Bs(17,k)
		  Bc%f_no1=Bs(18,k); Bc%L1=Bs(19,k); Bc%L2=Bs(20,k); Bc%L3=Bs(21,k)
         enddo
      else
!        ��nsub ��Bs ���ͳ�ȥ
	     Send_to_ID=B_proc(m)              ! ����Ŀ������ڵĽ��̺�
	     tag=B_n(m)                        ! ���
  	     call MPI_send(Nsub,1,MPI_INTEGER, Send_to_ID, tag, MPI_COMM_WORLD,ierr )        !������
  	     call MPI_send(Bs,size(Bs),MPI_INTEGER, Send_to_ID, tag, MPI_COMM_WORLD,ierr )    !����������Ϣ
     endif
	  deallocate(Bs)
    enddo
	  close(88)
   endif

! �Ǹ�����
   if(my_id .ne. 0) then
      do m=1,Mesh(1)%Num_Block
	    B=>Mesh(1)%Block(m)
	    call MPI_Recv(Nsub,1,MPI_INTEGER,0,m,MPI_COMM_WORLD,status,ierr)
 	    allocate(Bs(Nc,Nsub))
 	    call MPI_Recv(Bs,Nsub*Nc,MPI_INTEGER,0,m,MPI_COMM_WORLD,status,ierr)
	    B%subface=Nsub   
      	allocate(B%bc_msg(B%subface))   ! �߽�����
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
	 call MPI_Barrier(MPI_COMM_WORLD,ierr)
	 if(my_id .eq. 0) print*, "read bc3d.inc OK"

  end  subroutine read_inc 

