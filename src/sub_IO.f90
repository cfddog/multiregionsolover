
 !
 !
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
	if( Mesh_File_Format .eq. 1) then !
     open(99,file="Mesh3d.dat")
     read(99,*) NB   ! Block number
     allocate( NI(NB),NJ(NB),NK(NB) )
     read(99,*) (NI(k), NJ(k), NK(k), k=1,NB)
    elseif( Mesh_File_Format .eq. 2) then !
!    Mesh3d.x format: binary, each block stores x,y,z in one unformatted record
     open(99,file="Mesh3d.x",form="unformatted")
     read(99) NB
     allocate( NI(NB),NJ(NB),NK(NB) )
     read(99) (NI(k), NJ(k), NK(k), k=1,NB)
    else !
     open(99,file="Mesh3d.dat",form="unformatted")
     read(99) NB !
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
     elseif( Mesh_File_Format .eq. 2) then
!      Mesh3d.x: each block stores x,y,z in three separate records
       read(99) (((Ux(i,j,k,1),i=1,nx),j=1,ny),k=1,nz)
       read(99) (((Ux(i,j,k,2),i=1,nx),j=1,ny),k=1,nz)
       read(99) (((Ux(i,j,k,3),i=1,nx),j=1,ny),k=1,nz)
	 else
       read(99)   (((Ux(i,j,k,1),i=1,nx),j=1,ny),k=1,nz) , &
                  (((Ux(i,j,k,2),i=1,nx),j=1,ny),k=1,nz) , &
                  (((Ux(i,j,k,3),i=1,nx),j=1,ny),k=1,nz)
	 endif
   if(B_proc(m) .eq. 0) then !
      mt=B_n(m) !
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
     else !
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
  
  else !
    do m=1,MP%Num_Block !
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
   if(my_id .eq. 0) then
     if(Mesh_File_Format .eq. 2) then
       print*, "read Mesh3d.x OK"
     else
       print*, "read Mesh3d.dat OK"
     endif
   endif
 end subroutine read_main_mesh

!-------------------------------------------------------------------------------------
 !
 !
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

 !
  if(my_id .eq. 0) then 
    print*, " read distance to the wall:  wall_dist.dat"
  
   open(99,file="wall_dist.dat",form="unformatted")

   do m=1,Total_block
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(dw(nx-1,ny-1,nz-1))
     read(99) (((dw(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1) 
     
	 if(B_proc(m) .eq. 0) then !
      mt=B_n(m) !
	  B=>MP%Block(mt)
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
	   B%dw(i,j,k)=dw(i,j,k)
	  enddo
	  enddo
 	  enddo
     else !
	   Num_data=(nx-1)*(ny-1)*(nz-1)
	   Send_to_ID=B_proc(m)
	   tag=B_n(m)
	   call MPI_send(dw,Num_data,OCFD_DATA_TYPE, Send_to_ID, tag, MPI_COMM_WORLD,ierr )
     endif
     deallocate(dw)
   enddo
   close(99)
  
  else !
    do m=1,MP%Num_Block !
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
 !
 !
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

 !
  if(my_id .eq. 0) then 
    print*, " read flow data:  flow3d.dat"
 
     open(99,file="flow3d.dat",form="unformatted")
    
 	 if(NVAR1 .eq. 6) then !
      Inquire(file="SA3d.dat",exist=Ex)
      if(Ex) then
        open(100,file="SA3d.dat",form="unformatted")
      endif
	 endif
     
	 if(NVAR1 .eq. 7) then !
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
 !
     
	 if(B_proc(m) .eq. 0) then !
      mt=B_n(m) !
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
     else !
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

  else !
   
    do m=1,MP%Num_Block !
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
 !  For BLOCK_LOWSPEED / BLOCK_POROUS blocks B%U already stores primitives
 !  (rho,u,v,w,T) and must NOT be converted to conservative variables.
 do m=1,MP%Num_Block !
    B=>MP%Block(m)
    nx=B%nx; ny=B%ny; nz=B%nz
   if(B%Block_type == BLOCK_LOWSPEED .or. B%Block_type == BLOCK_POROUS) cycle
 
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

 !  Restart files carry only U(1..5); re-seed pressure and solid-frame
 !  temperature for porous blocks on every process.
  do m=1,MP%Num_Block
    B=>MP%Block(m)
    if(B%Block_type /= BLOCK_POROUS) cycle
    nx=B%nx; ny=B%ny; nz=B%nz
    do k=1-LAP,nz+LAP-1; do j=1-LAP,ny+LAP-1; do i=1-LAP,nx+LAP-1
      B%p(i,j,k)  = LS_P_out
      B%Ts(i,j,k) = Porous_T_ref
    enddo; enddo; enddo
  enddo
 !
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
 !

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


 !
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
  
   do m=1, Total_block !
     
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz,NVAR1))

	if(B_proc(m) .eq. 0) then !
      mt=B_n(m) !
	  B=>MP%Block(mt)
	 
	  if(B%Block_type == BLOCK_LOWSPEED .or. B%Block_type == BLOCK_POROUS) then
	    do k=0,nz
	    do j=0,ny
	    do i=0,nx
		  U(i,j,k,1)=B%U(1,i,j,k)        ! rho
          U(i,j,k,2)=B%U(2,i,j,k)        ! u
          U(i,j,k,3)=B%U(3,i,j,k)        ! v
          U(i,j,k,4)=B%U(4,i,j,k)        ! w
          U(i,j,k,5)=B%U(5,i,j,k)        ! T
	    enddo
	    enddo
 	    enddo
	  else
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
	  endif

    else !
	   Num_data=NVAR1*(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m) !
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
 
 else !

    do m=1,MP%Num_Block !
      B=>MP%Block(m)
	  nx=B%nx; ny=B%ny; nz=B%nz
   	  allocate(U(0:nx,0:ny,0:nz,NVAR1))
	  Num_data=(nx+1)*(ny+1)*(nz+1)*NVAR1
 	  tag=m
	 
	   if(B%Block_type == BLOCK_LOWSPEED .or. B%Block_type == BLOCK_POROUS) then
	     do k=0,nz
	     do j=0,ny
	     do i=0,nx
		   U(i,j,k,1)=B%U(1,i,j,k)
           U(i,j,k,2)=B%U(2,i,j,k)
           U(i,j,k,3)=B%U(3,i,j,k)
           U(i,j,k,4)=B%U(4,i,j,k)
           U(i,j,k,5)=B%U(5,i,j,k)
	     enddo
	     enddo
 	     enddo
	   else
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
	   endif
	   call MPI_Send(U,Num_data,OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD,ierr )
      deallocate(U)
    enddo
   
  endif
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0)  print*, "write flow3d.dat OK"

  end subroutine output_flow

!----------------------------------------------------------------------
! Output VTK format (legacy STRUCTURED_GRID) for Paraview
! Writes one file per block:
!   flow3d_block_NNN.vtk  - fluid blocks (cell-centered d,u,v,w,T,p)
!   Ts_block_NNN.vtk      - solid blocks (cell-centered temperature)
!----------------------------------------------------------------------
  subroutine output_vtk
   use Global_Var
   use const_var
   implicit none
   
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:,:,:):: U, G
   real(PRE_EC),allocatable,dimension(:,:,:):: Ts_buf
   integer:: NB,NVAR1,m,m1,nx,ny,nz,i,j,k,mt,Num_data
   integer:: Recv_from_ID,tag,ierr, status(MPI_status_size)
   integer:: npts, ncells, grid_size, flow_size
   real(PRE_EC):: d1, u1, v1, w1, p1, T1
   real(PRE_EC),parameter:: R_AIR_SI=287.0d0
   real(PRE_EC),parameter:: MU_SI0=1.716d-5, T_SI0=273.15d0, S_SI=110.4d0
   real(PRE_EC):: sc_rho, sc_u, sc_T, sc_p, a_ref, mu_inf_p
   character(len=50):: filename

   MP=>Mesh(1)
   NVAR1=MP%NVAR

!  Optional SI conversion for compressible (BLOCK_FLUID) blocks in VTK output:
!  U_ref=Ma*a_ref, rho_ref=Re*mu_SI(T_inf)/(Ma*a_ref*Lscale), T_ref=T_inf,
!  p_ref=rho_ref*U_ref^2  (real-gas reference state).
   sc_rho=1.d0; sc_u=1.d0; sc_T=1.d0; sc_p=1.d0
   if(Iflag_vtk_SI .eq. 1) then
     a_ref = sqrt(gamma*R_AIR_SI*T_inf)
     mu_inf_p = MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
     sc_u   = Ma*a_ref
     sc_rho = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
     sc_T   = T_inf
     sc_p   = sc_rho*sc_u*sc_u
   endif

   if(my_id .eq. 0) then
     print*, "write VTK files ......"
     
!    --- Flow data VTK (all blocks) ---
     if(Iflag_vtk_onefile .eq. 0) then
!      default: one structured-grid file per block
     do m=1, Total_block
       nx=bNi(m); ny=bNj(m); nz=bNk(m)
       npts = nx*ny*nz
       ncells = (nx-1)*(ny-1)*(nz-1)
       allocate(U(0:nx,0:ny,0:nz,6))
       allocate(G(nx,ny,nz,3))

!      Gather grid coords + flow data to master process
       if(B_proc(m) .eq. 0) then
         mt=B_n(m)
         B=>MP%Block(mt)
!        Grid coordinates
         do k=1,nz; do j=1,ny; do i=1,nx
           G(i,j,k,1)=B%x(i,j,k)
           G(i,j,k,2)=B%y(i,j,k)
           G(i,j,k,3)=B%z(i,j,k)
         enddo; enddo; enddo
!        Flow data (cell-centered primitive variables)
          if(Block_Type_List(m) == BLOCK_LOWSPEED .or. Block_Type_List(m) == BLOCK_POROUS) then
            do k=0,nz; do j=0,ny; do i=0,nx
              U(i,j,k,1) = B%U(1,i,j,k)      ! rho
              U(i,j,k,2) = B%U(2,i,j,k)      ! u
              U(i,j,k,3) = B%U(3,i,j,k)      ! v
              U(i,j,k,4) = B%U(4,i,j,k)      ! w
              U(i,j,k,5) = B%U(5,i,j,k)      ! T
              U(i,j,k,6) = B%p(i,j,k)        ! p
            enddo; enddo; enddo
          else
            do k=0,nz; do j=0,ny; do i=0,nx
              d1 = B%U(1,i,j,k)
              u1 = B%U(2,i,j,k)/max(d1, 1.d-20)
              v1 = B%U(3,i,j,k)/max(d1, 1.d-20)
              w1 = B%U(4,i,j,k)/max(d1, 1.d-20)
              p1 = (B%U(5,i,j,k) - 0.5d0*d1*(u1*u1+v1*v1+w1*w1)) * (gamma-1.d0)
              T1 = gamma * Ma * Ma * p1 / max(d1, 1.d-20)
              U(i,j,k,1) = d1
              U(i,j,k,2) = u1
              U(i,j,k,3) = v1
              U(i,j,k,4) = w1
              U(i,j,k,5) = T1
              U(i,j,k,6) = p1
            enddo; enddo; enddo
          endif
       else
         grid_size = nx*ny*nz*3
         flow_size = 6*(nx+1)*(ny+1)*(nz+1)
         Recv_from_ID = B_proc(m)
         tag = B_n(m)*2
         call MPI_Recv(G, grid_size, OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD, Status, ierr)
         tag = B_n(m)*2+1
         call MPI_Recv(U, flow_size, OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD, Status, ierr)
       endif

!      SI conversion: compressible (non-dimensional) block -> physical units
       if(Iflag_vtk_SI .eq. 1) then
         if(Block_Type_List(m) /= BLOCK_LOWSPEED .and. Block_Type_List(m) /= BLOCK_POROUS) then
           do k=0,nz; do j=0,ny; do i=0,nx
             U(i,j,k,1) = U(i,j,k,1)*sc_rho
             U(i,j,k,2) = U(i,j,k,2)*sc_u
             U(i,j,k,3) = U(i,j,k,3)*sc_u
             U(i,j,k,4) = U(i,j,k,4)*sc_u
             U(i,j,k,5) = U(i,j,k,5)*sc_T
             U(i,j,k,6) = U(i,j,k,6)*sc_p
           enddo; enddo; enddo
         endif
       endif

!      Write VTK file for this block
       write(filename, '("flow3d_block_",I0,".vtk")') m
       open(99, file=filename, status='replace')
       write(99, '(A)') '# vtk DataFile Version 3.0'
       write(99, '(A)') 'OpenCFD-EC output'
       write(99, '(A)') 'ASCII'
       write(99, '(A)') 'DATASET STRUCTURED_GRID'
       write(99, '(A,3I6)') 'DIMENSIONS ', nx, ny, nz
       write(99, '(A,I12,A)') 'POINTS ', npts, ' double'
!      Grid nodes
       do k=1,nz; do j=1,ny; do i=1,nx
         write(99, '(3ES25.15)') G(i,j,k,1), G(i,j,k,2), G(i,j,k,3)
       enddo; enddo; enddo
!      Cell-centered data
       write(99, '(A,I12)') 'CELL_DATA ', ncells
!      Density
       write(99, '(A)') 'SCALARS density double 1'
       write(99, '(A)') 'LOOKUP_TABLE default'
       do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
         write(99, '(ES25.15)') U(i,j,k,1)
       enddo; enddo; enddo
!      Velocity (vector)
       write(99, '(A)') 'VECTORS velocity double'
       do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
         write(99, '(3ES25.15)') U(i,j,k,2), U(i,j,k,3), U(i,j,k,4)
       enddo; enddo; enddo
!      Temperature
       write(99, '(A)') 'SCALARS temperature double 1'
       write(99, '(A)') 'LOOKUP_TABLE default'
       do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
         write(99, '(ES25.15)') U(i,j,k,5)
       enddo; enddo; enddo
!      Pressure
       write(99, '(A)') 'SCALARS pressure double 1'
       write(99, '(A)') 'LOOKUP_TABLE default'
       do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
         write(99, '(ES25.15)') U(i,j,k,6)
       enddo; enddo; enddo

       close(99)
       deallocate(U, G)
     enddo
     else
!      Iflag_vtk_onefile = 1: all flow blocks written into one unstructured file
       call output_vtk_merged_flow
     endif


!    --- Solid temperature VTK (all solid blocks) ---
     do m=1, Total_block
       nx=bNi(m); ny=bNj(m); nz=bNk(m)
       npts = nx*ny*nz
       ncells = (nx-1)*(ny-1)*(nz-1)

!      Determine if this is a solid block (check via Block_Type_List)
       if(Block_Type_List(m) /= BLOCK_SOLID .and. Block_Type_List(m) /= BLOCK_POROUS) cycle

       allocate(G(nx,ny,nz,3))
       allocate(Ts_buf(nx-1,ny-1,nz-1))

!      Gather grid + Ts data
       if(B_proc(m) .eq. 0) then
         mt=B_n(m)
         B=>MP%Block(mt)
         do k=1,nz; do j=1,ny; do i=1,nx
           G(i,j,k,1)=B%x(i,j,k)
           G(i,j,k,2)=B%y(i,j,k)
           G(i,j,k,3)=B%z(i,j,k)
         enddo; enddo; enddo
         do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
           Ts_buf(i,j,k)=B%Ts(i,j,k)
         enddo; enddo; enddo
       else
         grid_size = nx*ny*nz*3
         Recv_from_ID = B_proc(m)
         tag = B_n(m)*4
         call MPI_Recv(G, grid_size, OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD, Status, ierr)
         tag = B_n(m)*4+1
         call MPI_Recv(Ts_buf, (nx-1)*(ny-1)*(nz-1), OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD, Status, ierr)
       endif

       write(filename, '("Ts_block_",I0,".vtk")') m
       open(99, file=filename, status='replace')
       write(99, '(A)') '# vtk DataFile Version 3.0'
       write(99, '(A)') 'OpenCFD-EC solid temperature'
       write(99, '(A)') 'ASCII'
       write(99, '(A)') 'DATASET STRUCTURED_GRID'
       write(99, '(A,3I6)') 'DIMENSIONS ', nx, ny, nz
       write(99, '(A,I12,A)') 'POINTS ', npts, ' double'
       do k=1,nz; do j=1,ny; do i=1,nx
         write(99, '(3ES25.15)') G(i,j,k,1), G(i,j,k,2), G(i,j,k,3)
       enddo; enddo; enddo
       write(99, '(A,I12)') 'CELL_DATA ', ncells
       write(99, '(A)') 'SCALARS temperature double 1'
       write(99, '(A)') 'LOOKUP_TABLE default'
       do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
         write(99, '(ES25.15)') Ts_buf(i,j,k)
       enddo; enddo; enddo
       close(99)
       deallocate(G, Ts_buf)
     enddo

   else
!    --- Non-master processes: send data for local blocks ---
     do m=1, MP%Num_Block
       B=>MP%Block(m)
       nx=B%nx; ny=B%ny; nz=B%nz
       mt = m  ! local block index = B_n(global_block_index)

!      Find global block index for this local block
       do m1=1, Total_block
         if(B_proc(m1) .eq. my_id .and. B_n(m1) .eq. m) exit
       enddo
       if(m1 .gt. Total_block) cycle

!      Send grid coordinates
       allocate(G(nx,ny,nz,3))
       do k=1,nz; do j=1,ny; do i=1,nx
         G(i,j,k,1)=B%x(i,j,k)
         G(i,j,k,2)=B%y(i,j,k)
         G(i,j,k,3)=B%z(i,j,k)
       enddo; enddo; enddo
       grid_size = nx*ny*nz*3
       tag = m*2
       call MPI_Send(G, grid_size, OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD, ierr)
       deallocate(G)

!      Send flow data (cell-centered primitive variables)
       allocate(U(0:nx,0:ny,0:nz,6))
       if(Block_Type_List(m1) == BLOCK_LOWSPEED .or. Block_Type_List(m1) == BLOCK_POROUS) then
         do k=0,nz; do j=0,ny; do i=0,nx
           U(i,j,k,1) = B%U(1,i,j,k)
           U(i,j,k,2) = B%U(2,i,j,k)
           U(i,j,k,3) = B%U(3,i,j,k)
           U(i,j,k,4) = B%U(4,i,j,k)
           U(i,j,k,5) = B%U(5,i,j,k)
           U(i,j,k,6) = B%p(i,j,k)
         enddo; enddo; enddo
       else
         do k=0,nz; do j=0,ny; do i=0,nx
           d1 = B%U(1,i,j,k)
           u1 = B%U(2,i,j,k)/max(d1, 1.d-20)
           v1 = B%U(3,i,j,k)/max(d1, 1.d-20)
           w1 = B%U(4,i,j,k)/max(d1, 1.d-20)
           p1 = (B%U(5,i,j,k) - 0.5d0*d1*(u1*u1+v1*v1+w1*w1)) * (gamma-1.d0)
           T1 = gamma * Ma * Ma * p1 / max(d1, 1.d-20)
           U(i,j,k,1) = d1
           U(i,j,k,2) = u1
           U(i,j,k,3) = v1
           U(i,j,k,4) = w1
           U(i,j,k,5) = T1
           U(i,j,k,6) = p1
         enddo; enddo; enddo
       endif
       flow_size = 6*(nx+1)*(ny+1)*(nz+1)
       tag = m*2+1
       call MPI_Send(U, flow_size, OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD, ierr)
       deallocate(U)

!      Send solid Ts data if this is a solid block
       if(B%Block_type == BLOCK_SOLID .or. B%Block_type == BLOCK_POROUS) then
         allocate(G(nx,ny,nz,3))
         do k=1,nz; do j=1,ny; do i=1,nx
           G(i,j,k,1)=B%x(i,j,k)
           G(i,j,k,2)=B%y(i,j,k)
           G(i,j,k,3)=B%z(i,j,k)
         enddo; enddo; enddo
         allocate(Ts_buf(nx-1,ny-1,nz-1))
         do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
           Ts_buf(i,j,k) = B%Ts(i,j,k)
         enddo; enddo; enddo
         tag = m*4
         call MPI_Send(G, grid_size, OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD, ierr)
         tag = m*4+1
         call MPI_Send(Ts_buf, (nx-1)*(ny-1)*(nz-1), OCFD_DATA_TYPE, 0, tag, MPI_COMM_WORLD, ierr)
         deallocate(G, Ts_buf)
       endif
     enddo
   endif
   
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(my_id .eq. 0) print*, "write VTK files OK"
  end subroutine output_vtk

  subroutine output_vtk_merged_flow
   use Global_Var
   use const_var
   implicit none
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real(PRE_EC),allocatable,dimension(:,:):: xyz
   real(PRE_EC),allocatable,dimension(:):: dbuf, tbuf, pbuf
   real(PRE_EC),allocatable,dimension(:,:):: ubuf
   real(PRE_EC),allocatable,dimension(:,:,:,:):: G
   real(PRE_EC),allocatable,dimension(:,:,:,:):: U
   integer:: m, mt, nx,ny,nz, i,j,k, pt_base, cl_base, idx, nb
   integer:: tot_pts, tot_cells, Recv_from_ID, tag, ierr
   integer:: grid_size, flow_size
   integer:: status(MPI_status_size)
   integer:: ip0,ip1,ip2,ip3,ip4,ip5,ip6,ip7
   real(PRE_EC):: d1,u1,v1,w1,p1,T1
   real(PRE_EC),parameter:: R_AIR_SI=287.0d0
   real(PRE_EC),parameter:: MU_SI0=1.716d-5, T_SI0=273.15d0, S_SI=110.4d0
   real(PRE_EC):: sc_rho, sc_u, sc_T, sc_p, a_ref, mu_inf_p

!  This subroutine runs on the master process only and writes one legacy VTK
!  UNSTRUCTURED_GRID file (flow3d.vtk) containing ALL flow-type blocks
!  (fluid / low-speed / porous).  Node and cell ordering follows the global
!  block sequence, so the file layout (POINTS, CELLS, CELL_TYPES, CELL_DATA)
!  can be reproduced from the per-block structured grid topology.

   MP=>Mesh(1)

!  SI conversion scales for compressible blocks (see output_vtk for definition)
   sc_rho=1.d0; sc_u=1.d0; sc_T=1.d0; sc_p=1.d0
   if(Iflag_vtk_SI .eq. 1) then
     a_ref = sqrt(gamma*R_AIR_SI*T_inf)
     mu_inf_p = MU_SI0*sqrt((T_inf/T_SI0)**3)*(T_SI0+S_SI)/(T_inf+S_SI)
     sc_u   = Ma*a_ref
     sc_rho = Re*mu_inf_p/(Ma*a_ref*max(Lscale,1.d-30))
     sc_T   = T_inf
     sc_p   = sc_rho*sc_u*sc_u
   endif

   tot_pts=0; tot_cells=0
   do m=1, Total_block
     tot_pts  = tot_pts  + bNi(m)*bNj(m)*bNk(m)
     tot_cells= tot_cells + (bNi(m)-1)*(bNj(m)-1)*(bNk(m)-1)
   enddo
   allocate(xyz(3,tot_pts))
   allocate(dbuf(tot_cells), ubuf(3,tot_cells), tbuf(tot_cells), pbuf(tot_cells))

   pt_base=0; cl_base=0
   do m=1, Total_block
     nx=bNi(m); ny=bNj(m); nz=bNk(m)
     allocate(U(0:nx,0:ny,0:nz,6))
     allocate(G(nx,ny,nz,3))

     if(B_proc(m) .eq. 0) then
       mt=B_n(m)
       B=>MP%Block(mt)
       do k=1,nz; do j=1,ny; do i=1,nx
         G(i,j,k,1)=B%x(i,j,k)
         G(i,j,k,2)=B%y(i,j,k)
         G(i,j,k,3)=B%z(i,j,k)
       enddo; enddo; enddo
       if(Block_Type_List(m) == BLOCK_LOWSPEED .or. Block_Type_List(m) == BLOCK_POROUS) then
         do k=0,nz; do j=0,ny; do i=0,nx
           U(i,j,k,1) = B%U(1,i,j,k)
           U(i,j,k,2) = B%U(2,i,j,k)
           U(i,j,k,3) = B%U(3,i,j,k)
           U(i,j,k,4) = B%U(4,i,j,k)
           U(i,j,k,5) = B%U(5,i,j,k)
           U(i,j,k,6) = B%p(i,j,k)
         enddo; enddo; enddo
       else
         do k=0,nz; do j=0,ny; do i=0,nx
           d1 = B%U(1,i,j,k)
           u1 = B%U(2,i,j,k)/max(d1, 1.d-20)
           v1 = B%U(3,i,j,k)/max(d1, 1.d-20)
           w1 = B%U(4,i,j,k)/max(d1, 1.d-20)
           p1 = (B%U(5,i,j,k) - 0.5d0*d1*(u1*u1+v1*v1+w1*w1)) * (gamma-1.d0)
           T1 = gamma * Ma * Ma * p1 / max(d1, 1.d-20)
           U(i,j,k,1) = d1
           U(i,j,k,2) = u1
           U(i,j,k,3) = v1
           U(i,j,k,4) = w1
           U(i,j,k,5) = T1
           U(i,j,k,6) = p1
         enddo; enddo; enddo
       endif
     else
       grid_size = nx*ny*nz*3
       flow_size = 6*(nx+1)*(ny+1)*(nz+1)
       Recv_from_ID = B_proc(m)
       tag = B_n(m)*2
       call MPI_Recv(G, grid_size, OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD, status, ierr)
       tag = B_n(m)*2+1
       call MPI_Recv(U, flow_size, OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD, status, ierr)
     endif

!    SI conversion: compressible (non-dimensional) block -> physical units
     if(Iflag_vtk_SI .eq. 1) then
       if(Block_Type_List(m) /= BLOCK_LOWSPEED .and. Block_Type_List(m) /= BLOCK_POROUS) then
         do k=0,nz; do j=0,ny; do i=0,nx
           U(i,j,k,1) = U(i,j,k,1)*sc_rho
           U(i,j,k,2) = U(i,j,k,2)*sc_u
           U(i,j,k,3) = U(i,j,k,3)*sc_u
           U(i,j,k,4) = U(i,j,k,4)*sc_u
           U(i,j,k,5) = U(i,j,k,5)*sc_T
           U(i,j,k,6) = U(i,j,k,6)*sc_p
         enddo; enddo; enddo
       endif
     endif

!    store node coordinates (block order)
     idx=0
     do k=1,nz; do j=1,ny; do i=1,nx
       xyz(1,pt_base+idx+1) = G(i,j,k,1)
       xyz(2,pt_base+idx+1) = G(i,j,k,2)
       xyz(3,pt_base+idx+1) = G(i,j,k,3)
       idx = idx + 1
     enddo; enddo; enddo

!    store cell-centred data (block order)
     idx=0
     do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
       dbuf(cl_base+idx+1)        = U(i,j,k,1)
       ubuf(1,cl_base+idx+1)      = U(i,j,k,2)
       ubuf(2,cl_base+idx+1)      = U(i,j,k,3)
       ubuf(3,cl_base+idx+1)      = U(i,j,k,4)
       tbuf(cl_base+idx+1)        = U(i,j,k,5)
       pbuf(cl_base+idx+1)        = U(i,j,k,6)
       idx = idx + 1
     enddo; enddo; enddo

     pt_base = pt_base + nx*ny*nz
     cl_base = cl_base + (nx-1)*(ny-1)*(nz-1)
     deallocate(G, U)
   enddo

   open(99, file='flow3d.vtk', status='replace')
   write(99, '(A)') '# vtk DataFile Version 3.0'
   write(99, '(A)') 'OpenCFD-EC output (all blocks in one file)'
   write(99, '(A)') 'ASCII'
   write(99, '(A)') 'DATASET UNSTRUCTURED_GRID'
   write(99, '(A,I12,A)') 'POINTS ', tot_pts, ' double'
   do nb=1, tot_pts
     write(99, '(3ES25.15)') xyz(1,nb), xyz(2,nb), xyz(3,nb)
   enddo

   write(99, '(A,I12,I13)') 'CELLS ', tot_cells, tot_cells*9
   pt_base = 0
   do m=1, Total_block
     nx=bNi(m); ny=bNj(m); nz=bNk(m)
     do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
       ip0 = pt_base + (k-1)*ny*nx + (j-1)*nx + (i-1)
       ip1 = ip0 + 1
       ip2 = ip0 + 1 + nx
       ip3 = ip0 + nx
       ip4 = ip0 + ny*nx
       ip5 = ip0 + ny*nx + 1
       ip6 = ip0 + ny*nx + 1 + nx
       ip7 = ip0 + ny*nx + nx
       write(99, '(A,8I12)') '8 ', ip0, ip1, ip2, ip3, ip4, ip5, ip6, ip7
     enddo; enddo; enddo
     pt_base = pt_base + nx*ny*nz
   enddo

   write(99, '(A,I12)') 'CELL_TYPES ', tot_cells
   do nb=1, tot_cells
     write(99, '(I4)') 12
   enddo

   write(99, '(A,I12)') 'CELL_DATA ', tot_cells
   write(99, '(A)') 'SCALARS density double 1'
   write(99, '(A)') 'LOOKUP_TABLE default'
   do nb=1, tot_cells
     write(99, '(ES25.15)') dbuf(nb)
   enddo
   write(99, '(A)') 'VECTORS velocity double'
   do nb=1, tot_cells
     write(99, '(3ES25.15)') ubuf(1,nb), ubuf(2,nb), ubuf(3,nb)
   enddo
   write(99, '(A)') 'SCALARS temperature double 1'
   write(99, '(A)') 'LOOKUP_TABLE default'
   do nb=1, tot_cells
     write(99, '(ES25.15)') tbuf(nb)
   enddo
   write(99, '(A)') 'SCALARS pressure double 1'
   write(99, '(A)') 'LOOKUP_TABLE default'
   do nb=1, tot_cells
     write(99, '(ES25.15)') pbuf(nb)
   enddo
   close(99)

   deallocate(xyz, dbuf, ubuf, tbuf, pbuf)
   print*, "write merged flow3d.vtk OK: points=", tot_pts, " cells=", tot_cells
  end subroutine output_vtk_merged_flow





!----------------------------------------------------------------------
 !

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
 
   do m=1, Total_block !
	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(0:nx,0:ny,0:nz))

	if(B_proc(m) .eq. 0) then !
      mt=B_n(m) !
	  B=>MP%Block(mt)
	 
	  do k=0,nz
	  do j=0,ny
	  do i=0,nx
		 U(i,j,k)=B%mu_t(i,j,k)*Re                     ! mu_t
	  enddo
	  enddo
 	  enddo
    else !
	   Num_data=(nx+1)*(ny+1)*(nz+1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m) !
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD,Status,ierr )
    endif
! write Data ....
    
	write(99) ((( U(i,j,k),i=0,nx),j=0,ny),k=0,nz)    

	deallocate(U)
   enddo
   close(99)
   close(100)
   close(101)
 
 else !

    do m=1,MP%Num_Block !
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
 !

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
  
   do m=1, Total_block !
 	 nx=bNi(m); ny=bNj(m); nz=bNk(m)
	 allocate(U(nx-1,ny-1,nz-1))

	if(B_proc(m) .eq. 0) then !
      mt=B_n(m) !
	  B=>MP%Block(mt)
	 
	  do k=1,nz-1
	  do j=1,ny-1
	  do i=1,nx-1
        U(i,j,k)=B%dw(i,j,k)
	  enddo
	  enddo
 	  enddo
    else !
	   Num_data=(nx-1)*(ny-1)*(nz-1)
	   Recv_from_ID=B_proc(m)
	   tag=B_n(m) !
 	  call MPI_Recv(U,Num_data,OCFD_DATA_TYPE, Recv_from_ID, tag, MPI_COMM_WORLD,Status,ierr )
    endif
! write Data ....
 	write(99) ((( U(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)    
	deallocate(U)
   enddo
   close(99)
 
 else !

    do m=1,MP%Num_Block !
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
 !
 !
 !
 !
 !
!--------------------------------------------------------------------------------------------------------------------
  subroutine read_inc 
   use Global_Var
   implicit none
   integer,parameter:: NC=21 !
   integer:: nx,ny,nz,NB,Nsub,m,mt,k,j
   integer:: Send_to_ID,tag,ierr,Status(MPI_Status_SIZE)
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer,pointer,dimension(:,:):: Bs
   logical:: ext_inc

 !
	if(my_id .eq. 0) then
	  inquire(file="bc3d.inp",exist=ext_inc)
	  if(ext_inc) then
	    call  convert_inp_inc
	  else
	    print*, "bc3d.inp not found, skip conversion"
	  endif
	endif


 !
   if(my_id .eq. 0) then
     print*, "read bc3d.inc (Link/boundary file)......"
     open(88,file="bc3d.inc")
     read(88,*)
     read(88,*) NB
 !
    do m=1,NB
     read(88,*) nx,ny,nz
     read(88,*)
     read(88,*) Nsub !
 	 allocate(Bs(NC,Nsub))
	 do k=1,Nsub
	  read(88,*) (Bs(j,k),j=1,9)
	  read(88,*) (Bs(j,k),j=10,21)
!	 read(88,*) Bc%ib,Bc%ie,Bc%jb,Bc%je,Bc%kb,Bc%ke,Bc%bc,Bc%face,Bc%f_no
!	 read(88,*) Bc%ib1,Bc%ie1,Bc%jb1,Bc%je1,Bc%kb1,Bc%ke1,Bc%nb1,Bc%face1,Bc%f_no1,Bc%L1,Bc%L2,Bc%L3
     enddo
 !
      
     if(B_Proc(m) .eq. 0) then
 !
       	 mt=B_n(m) !
		 B=>Mesh(1)%Block(mt)
		 B%subface=Nsub   
	     allocate(B%bc_msg(B%subface)) !
         do k=1,Nsub
		 Bc=>B%bc_msg(k)
          Bc%ib=Bs(1,k); Bc%ie=Bs(2,k); Bc%jb=Bs(3,k); Bc%je=Bs(4,k)
		  Bc%kb=Bs(5,k); Bc%ke=Bs(6,k); Bc%bc=Bs(7,k); Bc%face=Bs(8,k); Bc%f_no=Bs(9,k)

          Bc%ib1=Bs(10,k); Bc%ie1=Bs(11,k); Bc%jb1=Bs(12,k); Bc%je1=Bs(13,k)
		  Bc%kb1=Bs(14,k); Bc%ke1=Bs(15,k); Bc%nb1=Bs(16,k); Bc%face1=Bs(17,k)
		  Bc%f_no1=Bs(18,k); Bc%L1=Bs(19,k); Bc%L2=Bs(20,k); Bc%L3=Bs(21,k)
         enddo
      else
 !
	     Send_to_ID=B_proc(m) !
	     tag=B_n(m) !
  	     call MPI_send(Nsub,1,MPI_INTEGER, Send_to_ID, tag, MPI_COMM_WORLD,ierr ) !
  	     call MPI_send(Bs,size(Bs),MPI_INTEGER, Send_to_ID, tag, MPI_COMM_WORLD,ierr ) !
     endif
	  deallocate(Bs)
    enddo
	  close(88)
   endif

 !
   if(my_id .ne. 0) then
      do m=1,Mesh(1)%Num_Block
	    B=>Mesh(1)%Block(m)
	    call MPI_Recv(Nsub,1,MPI_INTEGER,0,m,MPI_COMM_WORLD,status,ierr)
 	    allocate(Bs(Nc,Nsub))
 	    call MPI_Recv(Bs,Nsub*Nc,MPI_INTEGER,0,m,MPI_COMM_WORLD,status,ierr)
	    B%subface=Nsub   
      	allocate(B%bc_msg(B%subface)) !
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

!----------------------------------------------------------------------
! Read interface connection information (bc3d_interface.inc)
! Similar to read_inc but populates bc_msg2 instead of bc_msg
!----------------------------------------------------------------------
  subroutine read_inc_interface
   use Global_Var
   implicit none
   integer,parameter:: NC=21
   integer:: nx,ny,nz,NB,Nsub,m,mt,k,j
   integer:: Send_to_ID,tag,ierr,Status(MPI_Status_SIZE)
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer,pointer,dimension(:,:):: Bs
   logical:: ex_if
   integer:: have_interface

   have_interface=0
   if(my_id .eq. 0) then
     inquire(file="bc3d_interface.inp",exist=ex_if)
     if(ex_if) then
       have_interface=1
     endif
   endif
   call MPI_bcast(have_interface,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

   if(have_interface .eq. 0) then
     if(my_id .eq. 0) print*, "bc3d_interface.inp not found, skip interface"
     return
   endif

   if(my_id .eq. 0) then
     call convert_inp_inc_interface
   endif

   if(my_id .eq. 0) then
     print*, "read bc3d_interface.inc (Interface link file)......"
     open(88,file="bc3d_interface.inc")
     read(88,*)
     read(88,*) NB
     do m=1,NB
       read(88,*) nx,ny,nz
       read(88,*)
       read(88,*) Nsub
       allocate(Bs(NC,Nsub))
       do k=1,Nsub
         read(88,*) (Bs(j,k),j=1,9)
         read(88,*) (Bs(j,k),j=10,21)
       enddo
       if(B_Proc(m) .eq. 0) then
         mt=B_n(m)
         B=>Mesh(1)%Block(mt)
         B%subface=Nsub
         allocate(B%bc_msg2(B%subface))
         do k=1,Nsub
           Bc=>B%bc_msg2(k)
           Bc%ib=Bs(1,k); Bc%ie=Bs(2,k); Bc%jb=Bs(3,k); Bc%je=Bs(4,k)
           Bc%kb=Bs(5,k); Bc%ke=Bs(6,k); Bc%bc=Bs(7,k); Bc%face=Bs(8,k); Bc%f_no=Bs(9,k)
           Bc%ib1=Bs(10,k); Bc%ie1=Bs(11,k); Bc%jb1=Bs(12,k); Bc%je1=Bs(13,k)
           Bc%kb1=Bs(14,k); Bc%ke1=Bs(15,k); Bc%nb1=Bs(16,k); Bc%face1=Bs(17,k)
           Bc%f_no1=Bs(18,k); Bc%L1=Bs(19,k); Bc%L2=Bs(20,k); Bc%L3=Bs(21,k)
         enddo
       else
         Send_to_ID=B_proc(m)
         tag=B_n(m)
         call MPI_send(Nsub,1,MPI_INTEGER, Send_to_ID, tag, MPI_COMM_WORLD,ierr)
         call MPI_send(Bs,size(Bs),MPI_INTEGER, Send_to_ID, tag, MPI_COMM_WORLD,ierr)
       endif
       deallocate(Bs)
     enddo
     close(88)
   endif

   if(my_id .ne. 0) then
     do m=1,Mesh(1)%Num_Block
       B=>Mesh(1)%Block(m)
       call MPI_Recv(Nsub,1,MPI_INTEGER,0,m,MPI_COMM_WORLD,status,ierr)
       allocate(Bs(Nc,Nsub))
       call MPI_Recv(Bs,Nsub*Nc,MPI_INTEGER,0,m,MPI_COMM_WORLD,status,ierr)
       B%subface=Nsub
       allocate(B%bc_msg2(B%subface))
       do k=1,Nsub
         Bc=>B%bc_msg2(k)
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
   if(my_id .eq. 0) print*, "read bc3d_interface.inc OK"
   call align_interface_to_bc_msg   ! reorder bc_msg2 to match bc_msg face order
  end subroutine read_inc_interface

!----------------------------------------------------------------------
! Reorder each local block's bc_msg2 entries so that bc_msg2(k) refers to
! the same physical face as bc_msg(k).  Auto-generated Gridgen files may
! list interface faces in a different order than bc3d.inp; all boundary
! loops compare bc_msg and bc_msg2 by the same ksub index, so the orders
! must match.  Matching key: face + (ib,ie,jb,je,kb,ke).
!----------------------------------------------------------------------
  subroutine align_interface_to_bc_msg
   use Global_Var
   use const_var
   implicit none
   integer:: m, k, l, kk, nsub, nbty
   Type (Block_TYPE),pointer:: B, Bn
   TYPE (BC_MSG_TYPE),pointer:: Bc1, Bc2
   integer,allocatable:: perm(:)
   logical,allocatable:: used(:)
   integer:: val1(7), val2(7)
   TYPE (BC_MSG_TYPE),allocatable:: tmp(:)

   do m=1, Mesh(1)%Num_Block
     B => Mesh(1)%Block(m)
     if(.not. associated(B%bc_msg2)) cycle
     if(.not. associated(B%bc_msg)) cycle
     nsub = B%subface
     if(size(B%bc_msg) /= nsub .or. size(B%bc_msg2) /= nsub) then
       print*, 'align_interface: size mismatch block', m, size(B%bc_msg), size(B%bc_msg2), nsub
       cycle
     endif
     allocate(perm(nsub), used(nsub))
     used(:) = .false.
     do k=1, nsub
       Bc1 => B%bc_msg(k)
       val1 = (/Bc1%face, Bc1%ib, Bc1%ie, Bc1%jb, Bc1%je, Bc1%kb, Bc1%ke/)
       perm(k) = k
       do l=1, nsub
         if(used(l)) cycle
         Bc2 => B%bc_msg2(l)
         val2 = (/Bc2%face, Bc2%ib, Bc2%ie, Bc2%jb, Bc2%je, Bc2%kb, Bc2%ke/)
         if(all(val1 == val2)) then
           perm(k) = l; used(l) = .true.
           exit
         endif
       enddo
     enddo
!    Apply permutation: bc_msg2(k) <- old bc_msg2(perm(k))
     allocate(tmp(nsub))
     tmp(:) = B%bc_msg2(:)
     do k=1, nsub
       B%bc_msg2(k) = tmp(perm(k))
     enddo
     deallocate(tmp, perm, used)

!    bc3d.inp is authoritative for boundary TYPES: bc_msg2 only supplies
!    the interface pairing (ranges, face1, L1..L3, nb1) from the
!    auto-generated bc3d_interface file, whose own bc codes are ignored.
      do k=1, nsub
        B%bc_msg2(k)%bc = B%bc_msg(k)%bc
!       Wall-placeholder pairing faces (physical code 2 with a real interface
!       connection in bc3d_interface) are promoted to the canonical cross-class
!       interface code so every is_interface_bc consumer (boundary ghost fill,
!       solid ghost, LS/AC, couples, staggered drivers) treats them as interfaces.
!       Grids already typing the face as an interface in bc3d.inp are untouched.
        Bc2 => B%bc_msg2(k)
        if(Bc2%nb1 > 0 .and. Bc2%face1 > 0 .and. .not. is_interface_bc(Bc2%bc)) then
          nbty = -1
          do kk=1, Mesh(1)%Num_Block
            if(Mesh(1)%Block(kk)%Block_no == Bc2%nb1) then
              Bn => Mesh(1)%Block(kk); nbty = Bn%Block_type; exit
            endif
          enddo
          if(nbty >= 0) then
            if((B%Block_type == BLOCK_FLUID   .and. nbty == BLOCK_SOLID) .or. &
               (B%Block_type == BLOCK_SOLID   .and. nbty == BLOCK_FLUID)) Bc2%bc = BC_Interface_FluidSolid
            if((B%Block_type == BLOCK_LOWSPEED .and. nbty == BLOCK_SOLID) .or. &
               (B%Block_type == BLOCK_SOLID   .and. nbty == BLOCK_LOWSPEED)) Bc2%bc = BC_Interface_LowSolid
            if((B%Block_type == BLOCK_FLUID   .and. nbty == BLOCK_LOWSPEED) .or. &
               (B%Block_type == BLOCK_LOWSPEED .and. nbty == BLOCK_FLUID)) Bc2%bc = BC_Interface_FluidLow
            if((B%Block_type == BLOCK_FLUID   .and. nbty == BLOCK_POROUS) .or. &
               (B%Block_type == BLOCK_POROUS  .and. nbty == BLOCK_FLUID)) Bc2%bc = BC_Interface_FluidPorous
            if((B%Block_type == BLOCK_LOWSPEED .and. nbty == BLOCK_POROUS) .or. &
               (B%Block_type == BLOCK_POROUS  .and. nbty == BLOCK_LOWSPEED)) Bc2%bc = BC_Interface_LowPorous
            if((B%Block_type == BLOCK_SOLID   .and. nbty == BLOCK_POROUS) .or. &
               (B%Block_type == BLOCK_POROUS  .and. nbty == BLOCK_SOLID)) Bc2%bc = BC_Interface_SolidPorous
          endif
        endif
      enddo
    enddo


  end subroutine align_interface_to_bc_msg 
