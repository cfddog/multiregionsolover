! Initialization: create data structures and assign initial values
! For multi-grid, create each level of mesh based on the parent mesh information
!  Check if the mesh is suitable for multi-grid
!  Single-direction grid points = 2*K+1 supports 2-level mesh, =4*K+1 supports 3-level mesh, =8*K+1 supports 4-level mesh ...  
!---------------------------------------------------------------------------------
!------------------------------------------------------------------------------     
  subroutine init
   use Global_var
   implicit none
   integer :: i,j,k,m,nx1,ny1,nz1,Num_Block1,ksub,Kmax_grid
   real(PRE_EC),allocatable,dimension(:,:,:):: xc,yc,zc
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
 !--------------------------------------------------------------------
 ! initial of const variables
   Ralfa(1)=1.d0 ;  Ralfa(2)=3.d0/4.d0 ; Ralfa(3)=1.d0/3.d0
   Rbeta(1)=1.d0 ;  Rbeta(2)=1.d0/4.d0 ; Rbeta(3)=2.d0/3.d0
   Rgamma(1)=0.d0;  Rgamma(2)=1.d0/4.d0; Rgamma(3)=2.d0/3.d0
   Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma)
!--------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
   call partation                         ! Domain decomposition (determine which process each block belongs to)
   call read_material_in
   call bcast_material
   allocate( Mesh(Num_Mesh) )             ! Main data structure: "Mesh" (whose members are "Mesh Blocks"). Mesh(1) is the finest mesh in multi-grid, Mesh(2), Mesh(3) are coarser meshes.
   call Creat_main_Mesh                   ! Create main mesh (finest mesh in multi-grid) (from mesh file Mesh3d.dat)
   call read_main_Mesh                    ! Read main mesh
   call set_block_type                    ! Set block type from Block_Type_List (after Mesh is created)
   call read_inc    ! Read mesh connectivity information (bc3d.inc)
   call read_inc_interface
   call read_solid_bc  ! Read solid thermal BCs (solid_bc.inp)
   call Update_coordinate_buffer_onemesh(1)
   call Comput_Goemetric_var(1)
   call update_Mesh_Center(1)    ! Update ghost values of center coordinates (for periodic conditions)

   if(IF_Debug==1) call Output_mesh_debug                  ! Output mesh including ghost cells
  
    if(IF_Walldist ==  1)  then
	   call comput_dist_wall   ! Compute (or read) distance to wall
    else
	   if(my_id .eq. 0) print*, " Need not read wall_dist.dat "
	endif

   if(Num_Mesh .ge. 2) then
     call Creat_Mesh(1,2)                 ! Create mesh 2 (coarse mesh) based on mesh 1 (finest mesh) information
   endif
   if(Num_Mesh .ge. 3) then
     call Creat_Mesh(2,3)                 ! Create mesh 3 (coarsest mesh) based on mesh 2 (coarse mesh) information
   endif
   
   do m=1,Num_Mesh
     call set_BcK(m)   ! Set boundary indicators (2013-11)
   enddo 

! !!! FVM_FDM FVM_FDM !!!!         hybrid Finite-Difference/ Finite-Valume Method   
   call init_FDM
! !!!----------------------------------------------------------------------------
  end   


!--------------------------------------------------------------------------------------
!    Create data structure: finest mesh (store geometric quantities and conservative variables)
  subroutine Creat_main_Mesh
   use Global_var
   implicit none
   integer:: m,Num_Cell,ierr
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   
!   print*, "-----------------------------"
! ---------node Coordinates----------------------------------------  
!  Mesh file: PLOT3D format;   
   MP=>Mesh(1)
   MP%NVAR=NVAR   ! Number of variables on the finest mesh
   MP%Num_Block=Num_Block                      ! Number of blocks included in this MPI process
   allocate(MP%Block(Num_block))               ! ??????????�B
   call set_size_blocks                        ! ?څ?????��
   call allocate_mem_Blocks(1)                 ! ?????????????????

!  ?څ???????
   allocate(MP%Res_max(NVAR),MP%Res_rms(NVAR))   ! ???��?????????��? 
   MP%Kstep=0
   MP%tt=0.d0
   Num_Cell=0   ! ???????
    do m=1,Num_Block
    B => MP%Block(m)
    Num_Cell=Num_Cell+(B%nx-1)*(B%ny-1)*(B%nz-1)
	
	B%IF_OverLimit=0                       ! ??????????? ???????????

    enddo
    call MPI_ALLREDUCE(Num_Cell,MP%Num_Cell,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
	if (my_id .eq. 0) then
	  print*, "creat main mesh OK, Num_Cell=",MP%Num_Cell
    endif
  end   subroutine Creat_main_Mesh


!---------------------------------------------
! ?څ?????��(B%nx,B%ny,B%nz), ???????????Mesh3d.dat 
  subroutine set_size_blocks
   use Global_var
   implicit none
   integer:: m,mb
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   integer:: NB,NB1,k,ierr
!   integer,allocatable,dimension(:):: NI,NJ,NK
  
   MP=>Mesh(1)
   NB=Total_Block
   allocate(bNI(NB),bNJ(NB),bNK(NB))         ! ????????????
  
  if(my_id .eq. 0) then      ! ????????��?��????

   if( Mesh_File_Format .eq. 1) then   ! ??????
     open(99,file="Mesh3d.dat")
     read(99,*) NB1   ! Block number
     read(99,*) (bNI(k), bNJ(k), bNK(k), k=1,NB)
   elseif( Mesh_File_Format .eq. 2) then   ! Mesh3d.x binary
     open(99,file="Mesh3d.x",form="unformatted")
     read(99) NB1
     if(NB1 .ne. NB) then
       print*, "Warning !!! Block number Error !!!"
       print*, "please check 'partation.dat' ..."
       stop
     endif
     read(99) (bNI(k), bNJ(k), bNK(k), k=1,NB)
    else                                ! ???????
     open(99,file="Mesh3d.dat",form="unformatted")
     read(99) NB1                       ! ?????
      if(NB1 .ne. NB) then 
	     print*, "Warning !!! Block number Error !!!"
         print*, "please check 'partation.dat' ..."
		 stop
	  endif
	 read(99) (bNI(k), bNJ(k), bNK(k), k=1,NB)
    endif
    close(99)
   endif
   
   call MPI_bcast(bNI,NB,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
   call MPI_bcast(bNJ,NB,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
   call MPI_bcast(bNK,NB,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
   
   do m=1,MP%Num_Block   ! ??????????????
    B=> MP%Block(m)       ! ????
    mb= my_blocks(m)      ! ???
    B%block_no=mb         ! ???  
    B%mpi_id=my_id        ! ?????????
    B%nx=bNI(mb)
    B%ny=bNJ(mb)
    B%nz=bNK(mb)
    B%IFLAG_FVM_FDM=Method_FVM   ! ???????????? 
	B%IF_OverLimit=0

   enddo
!   deallocate(NI,NJ,NK)
!   print*, " define size ok ...", my_id

  end subroutine set_size_blocks

!------??????????? (????????????????????????)----------------------------------
  



! ??????????????????????????m2 
  subroutine Creat_Mesh(m1,m2)
   use Global_Var
   implicit none
   integer:: NB,NVAR1,m,m1,m2,ksub,nx,ny,nz,i,j,k,i1,j1,k1,Bsub,Num_Cell,ierr
   Type (Block_TYPE),pointer:: B1,B2
   TYPE (BC_MSG_TYPE),pointer:: Bc1,Bc2
   Type (Mesh_TYPE),pointer:: MP1,MP2
   MP1=>Mesh(m1)             ! ????????? ???????
   Mp2=>Mesh(m2)             ! ????????   ????????
   
   MP2%NVAR=5                ! ???????????????

   NB=MP1%Num_Block
   MP2%Num_Block=NB          !  ????m2??m1 ???????
   MP2%Mesh_no=m2            ! ?????
   MP2%Num_Cell=0      
   NVAR1=MP2%NVAR    ! NVAR1=5  ????????????????
   allocate(MP2%Res_max(NVAR1),MP2%Res_rms(NVAR1)) 
   allocate(MP2%Block(NB))   ! ??MP2?��?????????????�B
     Num_Cell=0
   do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)
	 B2%Block_no=B1%block_no
     nx=(B1%nx-1)/2+1        ! ??????????
	 ny=(B1%ny-1)/2+1
	 nz=(B1%nz-1)/2+1
	 B2%nx=nx
	 B2%ny=ny
	 B2%nz=nz    
	 Num_Cell=Num_Cell+(nx-1)*(ny-1)*(nz-1)          ! ???MP2???????????
     B2%IFLAG_FVM_FDM=Method_FVM   ! ???????????? 
 	 B2%IF_OverLimit=0           ! ??????????????????
  
    enddo
    call MPI_ALLREDUCE(Num_Cell,MP2%Num_Cell,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)


!    ??????????????????
!--------------------------------------------------
     call allocate_mem_Blocks(m2)           ! ???????
!------------------------------------------------
    do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)
!   ?څ????????????????????????????
     do k=1,B2%nz     
	   do j=1,B2%ny
	     do i=1,B2%nx
	       i1=2*i-1 ; j1=2*j-1 ;k1=2*k-1
	       B2%x(i,j,k)=B1%x(i1,j1,k1)         !??????????????????? ???????????????????????
           B2%y(i,j,k)=B1%y(i1,j1,k1)
		   B2%z(i,j,k)=B1%z(i1,j1,k1)
	     enddo
	   enddo
	 enddo
	 enddo
!-----------------------------------------------  
    do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)

!    ???????????
     Bsub=B1%subface        ! ??????
     B2%subface=Bsub
     allocate(B2%bc_msg(Bsub))
     do ksub=1, Bsub
	   Bc1=> B1%bc_msg(ksub)    ! ?????????????????
	   Bc2=> B2%bc_msg(ksub)    ! ????????????????
      
	   Bc2%ib=(Bc1%ib-1)/2+1   ! ?????????��???????
	   Bc2%ie=(Bc1%ie-1)/2+1   ! ?????????��???????
	   Bc2%jb=(Bc1%jb-1)/2+1
	   Bc2%je=(Bc1%je-1)/2+1
	   Bc2%kb=(Bc1%kb-1)/2+1
	   Bc2%ke=(Bc1%ke-1)/2+1
 	   
	   Bc2%ib1=(Bc1%ib1-1)/2+1   ! ?????????��???????
	   Bc2%ie1=(Bc1%ie1-1)/2+1   ! ?????????��???????
	   Bc2%jb1=(Bc1%jb1-1)/2+1
	   Bc2%je1=(Bc1%je1-1)/2+1
	   Bc2%kb1=(Bc1%kb1-1)/2+1
	   Bc2%ke1=(Bc1%ke1-1)/2+1
      
	   Bc2%bc=Bc1%bc            ! ??????? ??-1 ????��
	   Bc2%face=Bc1%face        ! ??????(1-6?????? i-,j-,k-,i+,j+,k+)
	   Bc2%f_no=Bc1%f_no        ! ?????
	   Bc2%nb1=Bc1%nb1          ! ?????
	   Bc2%face1=Bc1%face1      ! ???????????(1-6)
	   Bc2%f_no1=Bc1%f_no1      ! ????????????
	   Bc2%L1=Bc1%L1            ! ??????????
	   Bc2%L2=Bc1%L2
	   Bc2%L3=Bc1%L3
     enddo
   enddo

   call Update_coordinate_buffer_onemesh(m2)
   call Comput_Goemetric_var(m2)

   call update_Mesh_Center(m2)

   Mesh(m2)%Kstep=0
   Mesh(m2)%tt=0.d0
  
  end  subroutine Creat_Mesh

! ------------------------------------------------------------------------------

  subroutine Init_flow
    use Global_var
    implicit none
    integer:: i,j,k,m,m1,NVAR1
   logical:: restarted
    Type (Mesh_TYPE),pointer:: MP
    Type (Block_TYPE),pointer:: B	 
  
   if(Iflag_restart .ge. 0) call read_restart(restarted)
   if(.not. restarted) then
   if(Iflag_init .le. 0) then
	 call init_flow_zero                   ! ??????????????????????????? ????????????????????????
   else
	 call read_flow_data
!    restart file stores only U(1..5); re-seed pressure & frame temperature
     do m=1,Mesh(1)%Num_Block
       B => Mesh(1)%Block(m)
       if(B%Block_type /= BLOCK_POROUS) cycle
       do k=1-LAP,B%nz+LAP-1; do j=1-LAP,B%ny+LAP-1; do i=1-LAP,B%nx+LAP-1
         B%p(i,j,k)  = LS_P_out
         B%Ts(i,j,k) = Porous_T_ref
       enddo; enddo; enddo
     enddo
   endif
   endif
 
 !    n??n-1??????? ????????????????, ????LU-SGS??? ??????????????
    if(Time_Method .eq. Time_dual_LU_SGS) then             
      MP=> Mesh(1)
      NVAR1=MP%NVAR
      do m=1,MP%Num_Block
      B => MP%Block(m)                
	   do k=-1,B%nz+1
        do j=-1,B%ny+1
         do i=-1,B%nx+1
          do m1=1,NVAR1
           B%Un(m1,i,j,k)=B%U(m1,i,j,k)
		   B%Un1(m1,i,j,k)=B%U(m1,i,j,k)
		  enddo
		 enddo
		enddo
	   enddo
	  enddo		   
    endif   
  end subroutine Init_flow




!--------------------------------------------------------------------------------
! ????????????? ????????????????????????????????????????  
  subroutine init_flow_zero
   use Global_var
   implicit none
   real(PRE_EC):: d0,u0,v0,w0,p0,T0,vx,tmp,pin0
   integer:: i,j,k,m,step,nMesh,n
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B	 
!-------------------------------------------------------------------
   if(my_id .eq. 0) then      
    if( Iflag_init .eq. Init_By_FreeStream)then    ! ???????????????
      print*, "Initial by Free-stream flow ......"
    else  if( Iflag_init .eq. Init_By_Zeroflow) then
      print*, "Initial by Zero flow ......"                                          ! ?????????????
    endif
   endif

  
   MP=> Mesh(Num_Mesh)   ! Mesh(Num_Mesh) ??????????
   do m=1,MP%Num_Block
     B => MP%Block(m)                
     if(B%Block_type == BLOCK_LOWSPEED) then
!        low-speed block: primitive variables U=[rho,u,v,w,T], pressure in B%p
!        (U(2..4) store VELOCITY, not momentum: U(1) carries LS_rho)
        do k=1-LAP,B%nz+LAP-1
        do j=1-LAP,B%ny+LAP-1
        do i=1-LAP,B%nx+LAP-1
          B%U(1,i,j,k)=LS_rho
          B%U(2,i,j,k)=LS_U_in
          B%U(3,i,j,k)=LS_V_in
          B%U(4,i,j,k)=0.d0
          B%U(5,i,j,k)=LS_T_ref
          B%p(i,j,k)=LS_P_out
        enddo; enddo; enddo
        cycle
     endif
     if(B%Block_type == BLOCK_POROUS) then
!        porous block: same primitive variables U=[rho,u,v,w,Tf] as low-speed,
!        pressure in B%p, solid-frame temperature in B%Ts
        do k=1-LAP,B%nz+LAP-1
        do j=1-LAP,B%ny+LAP-1
        do i=1-LAP,B%nx+LAP-1
          B%U(1,i,j,k)=LS_rho
          B%U(2,i,j,k)=LS_rho*LS_U_in
          B%U(3,i,j,k)=LS_rho*LS_V_in
          B%U(4,i,j,k)=0.d0
          B%U(5,i,j,k)=LS_T_ref
          B%p(i,j,k)=LS_P_out
          B%Ts(i,j,k)=Porous_T_ref
        enddo; enddo; enddo
        cycle
     endif
     if(IF_TurboMachinary ==0  .and. IF_Innerflow==0 ) then    ! ??????

	  d0=1.d0
      p0=1.d0/(gamma*Ma*Ma)
     if( Iflag_init .eq. Init_By_FreeStream)then    ! ???????????????
       u0=cos(A_alfa)*cos(A_beta)
       v0=sin(A_alfa)*cos(A_beta) 
       w0=sin(A_beta)
     else  if( Iflag_init .eq. Init_By_Zeroflow) then                                          ! ?????????????
 	   u0=0.d0
       v0=0.d0
       w0=0.d0
     endif


	   do k=1-LAP,B%nz+LAP-1
       do j=1-LAP,B%ny+LAP-1
       do i=1-LAP,B%nx+LAP-1

           B%U(1,i,j,k)=d0
           B%U(2,i,j,k)=d0*u0
           B%U(3,i,j,k)=d0*v0
           B%U(4,i,j,k)=d0*w0
           B%U(5,i,j,k)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0+w0*w0)

           if(MP%NVAR .eq. 6) then
! see:       http://turbmodels.larc.nasa.gov/spalart.html
		     B%U(6,i,j,k)=5.d0                 ! ?څ?????????????5?? ??0.98c???��??

		   else if (MP%NVAR .eq. 7) then
		     B%U(6,i,j,k)=10.d0*Kt_Inf   ! ????? ????????????????10????
			 B%U(7,i,j,k)=Wt_Inf         ! ??????
           endif
       enddo
       enddo
	   enddo
     
	 else   ! ??????

    
	  if( Iflag_init .eq. Init_By_Zeroflow) then                                         
        vx=0.d0
		d0=1.d0
		p0=1.d0/(gamma*Ma*Ma)
	  else 
       if(P_outlet > 0) then       !            
	     pin0=1.d0/(gamma)       ! ??????
	     p0= P_OUTLET                  ! ??? ?????????????
 	     T0= (P_OUTLET/pin0)**((gamma-1.d0)/gamma)   ! ?????? (????=1, ?????=1,???=1/gamma)
	     d0= P_OUTLET/T0*gamma*Ma*Ma          ! ??????
         vx=sqrt(2.d0*Cp*(1.d0-T0))   ! ??????
       else
         p0=1.d0/(gamma*Ma*Ma)
		 d0=1.d0
		 vx=1.d0
	   endif
  	  endif


	  do k=1-LAP,B%nz+LAP-1
      do j=1-LAP,B%ny+LAP-1
      do i=1-LAP,B%nx+LAP-1
     
	       u0=vx
		   v0= Turbo_w*B%zc(i,j,k)   ! ?????? ???????????
		   w0= -Turbo_w*B%yc(i,j,k)   
           B%U(1,i,j,k)=d0
           B%U(2,i,j,k)=d0*u0
           B%U(3,i,j,k)=d0*v0
           B%U(4,i,j,k)=d0*w0
           B%U(5,i,j,k)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0+w0*w0)
          
            if(MP%NVAR .eq. 6) then
		     B%U(6,i,j,k)=5.d0                 ! ?څ?????????????5?? ??0.98c???��??
		    else if (MP%NVAR .eq. 7) then
		     B%U(6,i,j,k)=10.d0*Kt_Inf   ! ????? ????????????????10????
			 B%U(7,i,j,k)=Wt_Inf         ! ??????
            endif

       enddo
	   enddo
	   enddo

     endif

   enddo

   call Boundary_condition_onemesh(Num_Mesh)     ! ??????? ???څGhost Cell?????
   call update_buffer_onemesh(Num_Mesh)          ! ?????????????

!  Initial smoothing    ! ?????? (????Kstep_Init_Smooth????
   do n=1,Kstep_init_smooth
     if(my_id .eq. 0 .and. mod(n,10) .eq. 0) print*, "Initial smoothing",n
	 call smoothing_oneMesh(Num_Mesh,Smooth_2nd)     
   enddo
   
   
   
!-----------------------------------------------------------------
!------------------------------------------------------
!   ???????????
!   ???????????????????????
   do nMesh=Num_Mesh,1,-1 
     do step=1, Pre_Step_Mesh(nMesh)   
       call NS_Time_advance(nMesh)
       if(mod(step,Kstep_show) .eq. 0) call output_Res(nMesh)
     enddo
!     call output (nMesh)
     if(nMesh .gt. 1) then
       call prolong_U(nMesh,nMesh-1,1)                  ! ??nMesh?????????????????????????????; flag=1 ???U????
	   call Boundary_condition_onemesh(nMesh-1)         ! ??????? ???څGhost Cell?????
	   call update_buffer_onemesh(nMesh-1)              ! ????????????? 
       print*, " Prolong  to mesh ", nMesh-1, "   OK"           
     endif
   enddo

  end subroutine init_flow_zero 


!----------------------------------------------------


  subroutine allocate_mem_Blocks(nMesh)
   use Global_var
   implicit none
   integer:: nMesh,m,nx,ny,nz,NVAR1
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B

    MP=>Mesh(nMesh)
    NVAR1=MP%NVAR

   do m=1,MP%Num_Block
	 B=>MP%Block(m)

     nx=B%nx ; ny= B%ny  ; nz= B%nz

!    Nullify pointer components that may not be allocated later
     nullify(B%bc_msg2)
     nullify(B%solid_bc_face_no)
     nullify(B%solid_bc_Tw)
     nullify(B%solid_bc_Qw)
     nullify(B%solid_bc_htc)
     nullify(B%solid_bc_Tinf)
     B%solid_bc_nface = 0
	  	 
!   ???????   (x,y,z) ??????? (xc,yc,zc)????????????; Vol ??????????? U, Un ??????
     allocate( B%x(0:nx+1,0:ny+1,0:nz+1), B%y(0:nx+1,0:ny+1,0:nz+1),B%z(0:nx+1,0:ny+1,0:nz+1))   ! ???????
  
!  ???????? LAP????????  ???????????????????????          
   	 allocate( B%xc(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1) , &
	           B%yc(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1) , &
			   B%zc(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1)  )   ! LAP ????????   
	 
	 
	 allocate( B%U(NVAR1,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1) )   ! LAP ????????   ! bug is removed
	 allocate( B%Un(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))  ! ??? Ghost Cell
     allocate( B%Vol(nx-1,ny-1,nz-1)) 
     allocate( B%Res(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))        !  ?��?
	 allocate( B%dt(-1:nx+1,-1:ny+1,-1:nz+1))               !  ?????
	 allocate( B%deltU(5,-1:nx+1,-1:ny+1,-1:nz+1))          !  ?????U???? ?????????????????????????????5??????????
	 allocate( B%dU(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))         !  ?????U???? ??LU-SGS??????
     allocate( B%Si(nx,ny,nz), B%Sj(nx,ny,nz), B%Sk(nx,ny,nz) )  ! ?????
     allocate( B%ni1(nx,ny,nz),B%ni2(nx,ny,nz),B%ni3(nx,ny,nz), & 
               B%nj1(nx,ny,nz),B%nj2(nx,ny,nz),B%nj3(nx,ny,nz), &
               B%nk1(nx,ny,nz),B%nk2(nx,ny,nz),B%nk3(nx,ny,nz))
	 allocate( B%dw(nx-1,ny-1,nz-1))         ! ??????????
!  Jocabian?��???????????????????��???? 
     allocate(B%ix1(nx,ny,nz),B%iy1(nx,ny,nz),B%iz1(nx,ny,nz), &
	          B%jx1(nx,ny,nz),B%jy1(nx,ny,nz),B%jz1(nx,ny,nz), &
              B%kx1(nx,ny,nz),B%ky1(nx,ny,nz),B%kz1(nx,ny,nz))
     allocate(B%ix2(nx,ny,nz),B%iy2(nx,ny,nz),B%iz2(nx,ny,nz), &
	          B%jx2(nx,ny,nz),B%jy2(nx,ny,nz),B%jz2(nx,ny,nz), &
              B%kx2(nx,ny,nz),B%ky2(nx,ny,nz),B%kz2(nx,ny,nz))
     allocate(B%ix3(nx,ny,nz),B%iy3(nx,ny,nz),B%iz3(nx,ny,nz), &
	          B%jx3(nx,ny,nz),B%jy3(nx,ny,nz),B%jz3(nx,ny,nz), &
              B%kx3(nx,ny,nz),B%ky3(nx,ny,nz),B%kz3(nx,ny,nz))
     allocate(B%ix0(nx,ny,nz),B%iy0(nx,ny,nz),B%iz0(nx,ny,nz), &
	          B%jx0(nx,ny,nz),B%jy0(nx,ny,nz),B%jz0(nx,ny,nz), &
              B%kx0(nx,ny,nz),B%ky0(nx,ny,nz),B%kz0(nx,ny,nz))
!-------------------------------------------------------
      allocate(B%dtime_mesh(nx-1,ny-1,nz-1))   ! ????????? ??????????????????
               B%dtime_mesh(:,:,:)=1.d0                ! ???
	 
	 if(Time_Method .eq. Time_Dual_LU_SGS) then             ! ????LU_SGS??? n-1??????????
          allocate(B%Un1(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))
	 endif

!-------------------------------------------------------

     if(If_viscous .eq. 1) then 
      allocate(B%mu(-1:nx+1,-1:ny+1,-1:nz+1))   ! ??????????
	  allocate(B%mu_t(-1:nx+1,-1:ny+1,-1:nz+1))    ! ??????????
	  B%mu(:,:,:)=1.d0/Re
	  B%mu_t(:,:,:)=0.d0
     endif
!------??????  (??????????????????)------------------------------------
     allocate( B%Surf1(ny,nz,3),B%Surf2(nx,nz,3), B%Surf3(nx,ny,3),  &
               B%Surf4(ny,nz,3),B%Surf5(nx,nz,3), B%Surf6(nx,ny,3)   )          
	 B%Surf1(:,:,:)=0.d0; B%Surf2(:,:,:)=0.d0; B%Surf3(:,:,:)=0.d0
     B%Surf4(:,:,:)=0.d0; B%Surf5(:,:,:)=0.d0; B%Surf6(:,:,:)=0.d0
!---------------------------------------------------------------------------
! --------?????????????????????1??????????----------------------------
     B%x(:,:,:)=0.d0; B%y(:,:,:)=0.d0; B%z(:,:,:)=0.d0; B%xc(:,:,:)=0.d0; B%yc(:,:,:)=0.d0; B%zc(:,:,:)=0.d0 
     B%U(1,:,:,:)=1.d0; B%U(2,:,:,:)=0.d0; B%U(3,:,:,:)=0.d0; B%U(4,:,:,:)=0.d0; B%U(5,:,:,:)=1.d0
	 
	 if(NVAR1 .eq. 6) then
	    B%U(6,:,:,:)=1.d0
	 else if (NVAR1 .eq. 7) then
	    B%U(6,:,:,:)=0.d0
		B%U(7,:,:,:)=1.d0
	 endif
	   B%Res(:,:,:,:)=0.d0
   
    if( nMesh .ne. 1) then
       allocate( B%QF(NVAR1,-1:nx+1,-1:ny+1,-1:nz+1))         ! ??????
	    B%QF(:,:,:,:)=0.d0                                    ! ????????????0
    endif
 
 !----???????? (cell-centered, LAP ghost layers, ???���Y????)
     allocate(B%Ts(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
     allocate(B%Tsn(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
     B%Ts(:,:,:)=T_inf      ! ???????��????
     B%Tsn(:,:,:)=T_inf

!---- low-speed pressure field (cell-centered, LAP ghost layers)
     allocate(B%p(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1,1-LAP:nz+LAP-1))
     B%p(:,:,:)=LS_P_out

!----??????? (1 ??????��0????)
     allocate(B%BcI(ny-1,nz-1,2),B%BcJ(nx-1,nz-1,2),B%BcK(nx-1,ny-1,2))
     B%BcI(:,:,:)=0
         B%BcJ(:,:,:)=0
         B%BcK(:,:,:)=0
    
   enddo

  end subroutine allocate_mem_Blocks


 ! ?څ??????? (0 ??????�� 1 ????), ???? ????? ???????????????
   subroutine set_BcK(nm)  
   use Global_var
   implicit none
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
   integer :: i,j,k,m,nm,ksub
   integer:: ib,ie,jb,je,kb,ke
     do m=1,Mesh(nm)%Num_Block
       B=>Mesh(nm)%Block(m)
       B%BcI(:,:,:)=0
	   B%BcJ(:,:,:)=0
	   B%BcK(:,:,:)=0
     do  ksub=1,B%subface
     Bc=> B%bc_msg(ksub)
     ib=Bc%ib; ie=Bc%ie; jb=Bc%jb; je=Bc%je ; kb=Bc%kb; ke=Bc%ke      
     if(Bc%bc >=0 ) then   ! ??????
       if(Bc%face .eq. 1 ) then   
         B%BcI(jb:je-1,kb:ke-1,1)=1
	   else if(Bc%face .eq. 2) then
         B%BcJ(ib:ie-1,kb:ke-1,1)=1
	   else if(Bc%face .eq. 3) then
         B%BcK(ib:ie-1,jb:je-1,1)=1
       else if(Bc%face .eq. 4 ) then   
         B%BcI(jb:je-1,kb:ke-1,2)=1
	   else if(Bc%face .eq. 5) then
         B%BcJ(ib:ie-1,kb:ke-1,2)=1
	   else if(Bc%face .eq. 6) then
         B%BcK(ib:ie-1,jb:je-1,2)=1
       endif
	 endif
	enddo
    enddo
   end

!-----------------------------------------------------------------------
  subroutine bcast_material
   use Global_var
   implicit none
   integer:: ierr
   if(.not. associated(Block_Type_List)) then
     if(my_id .eq. 0) then
       allocate(Block_Type_List(Total_block))
       Block_Type_List(:)=BLOCK_FLUID
     endif
   endif
   if(.not. associated(solid_rho_list)) then
     if(my_id .eq. 0) then
       allocate(solid_rho_list(Total_block), solid_Cp_list(Total_block), solid_k_list(Total_block))
       solid_rho_list(:)=0.d0; solid_Cp_list(:)=0.d0; solid_k_list(:)=0.d0
     endif
   endif
   if(my_id .ne. 0) then
     allocate(Block_Type_List(Total_block))
     allocate(solid_rho_list(Total_block), solid_Cp_list(Total_block), solid_k_list(Total_block))
   endif
   call MPI_bcast(Block_Type_List,Total_block,MPI_Integer,0,MPI_COMM_WORLD,ierr)
   call MPI_bcast(solid_rho_list,Total_block,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   call MPI_bcast(solid_Cp_list,Total_block,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   call MPI_bcast(solid_k_list,Total_block,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   if(.not. associated(porous_eps_list)) then
     if(my_id .eq. 0) then
       allocate(porous_eps_list(Total_block), porous_dp_list(Total_block), porous_hv_list(Total_block))
       porous_eps_list(:)=0.9d0; porous_dp_list(:)=1.d-3; porous_hv_list(:)=0.d0
     endif
   endif
   if(my_id .ne. 0) then
     allocate(porous_eps_list(Total_block), porous_dp_list(Total_block), porous_hv_list(Total_block))
   endif
   call MPI_bcast(porous_eps_list,Total_block,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   call MPI_bcast(porous_dp_list,Total_block,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
   call MPI_bcast(porous_hv_list,Total_block,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
  end subroutine bcast_material

!-----------------------------------------------------------------------
  subroutine set_block_type
   use Global_var
   implicit none
   integer:: m,ib
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   MP=>Mesh(1)
   do m=1,MP%Num_Block
     B=>MP%Block(m)
     ib=B%Block_no
     if(associated(Block_Type_List)) then
       B%Block_type=Block_Type_List(ib)
     else
       B%Block_type=BLOCK_FLUID
     endif
     B%solid_rho=solid_rho_list(ib)
     B%solid_Cp=solid_Cp_list(ib)
     B%solid_k=solid_k_list(ib)
     B%porous_eps=porous_eps_list(ib)
     B%porous_dp=porous_dp_list(ib)
     B%porous_hv=porous_hv_list(ib)
   enddo
  end subroutine set_block_type
!------------------------------------------------