!-----------------------------------------------------------
! Copyright by LiXinliang
! Ver 1.2   ����������ת���������
! Ver 1.3 �ɶ�ȡflow3d.dat (tecplot��ʽ�ļ�)��flow.dat (�޸�ʽ�ļ�)
! Ver 1.4  ��ȡ ver 0.8 ���ϰ汾������ ������k,w,mu_t�ȣ�
! Ver 1.5 ��ȡver0.82���ϵ�����
! Ver 1.6 ��ȡver 0.84���ϰ汾�����ݣ�.inp�ļ���
! Ver 1.7 Plot 3D flow
! Ver 1.8 Plot Cp and Cf on the wall
! Ver 1.8a, vt �������Re
! Ver 2.3, �ɶ�ȡBC_user
! Ver 2.4a, ������read wall_dist() �е�Bug
! Ver 2.5, �����������ϵ����� 
!------------------------------------------------------------
  module Const_Variables
  implicit none
  integer,parameter:: PRE_EC=8
   integer,parameter::  BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6 
   integer,parameter::  BC_Wall_Turbo=201             
   integer,parameter::  BC_Periodic=501, BC_Extrapolate=401      
   integer,parameter::  BC_Inner=-1, BC_PeriodicL=-2, BC_PeriodicR=-3   
   integer,parameter::  BC_Zero=0     
   integer,parameter::  BC_USER_FixedInlet=901, BC_USER_Inlet_time=902      
   integer,parameter:: BC_USER_Blow_Suction_Wall=903    
end module Const_Variables

 
  module Global_Variables
  use Const_Variables
  implicit none
 ! global parameter (for all blocks) 
  real(PRE_EC),save:: Ma,Re,gamma,T_inf,p00
  integer,save:: Num_Block, Mesh_File_Format, BC_number,BC_type(100)
!-----------------------------------------------------------------------------------------
! ����������Ϣ 
  TYPE BC_MSG_TYPE
   integer::   ib,ie,jb,je,kb,ke,bc
  END TYPE BC_MSG_TYPE

!  ���ı�������ÿ������洢����Ϣ ��ȫ�ֱ�����
   TYPE Block_TYPE           !  variables for each block 
     integer :: Block_no,nx,ny,nz,subface
     real(PRE_EC),pointer,dimension(:,:,:):: xc,yc,zc  ! coordinates of cell center, ������������ 
     real(PRE_EC),pointer,dimension(:,:,:):: dc,uc,vc,wc,Tc,pc,mut,vt,kt,wt
	 real(PRE_EC),pointer,dimension(:,:,:):: d,u,v,w,T,p,mut1,vt1
	 real(PRE_EC),pointer,dimension(:,:,:):: x,y,z,dw
     TYPE(BC_MSG_TYPE),pointer,dimension(:):: bc_msg
   End TYPE Block_TYPE  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
  
  end module Global_Variables
!------------------------------------------------------------------


!-------------------------------------------
     program readflow   
     use Global_Variables
     implicit none
     integer:: If3d,m
	  character(len=50):: filename
    
	   call init
   
    print*, "Plot boundary plane ;  bounary-xxx.dat  Grid value, boundaryc-xxx.dat  Grid center value"
    print*, "xxx:  002 BC_Wall,  003BC_Symmetry, 004 BC_Farfield, 005 BC_Inflow, 006 BC_Outflow"
     
     do m=1,Bc_number
           write(filename,"('boundary-'I3.3'.dat')") Bc_type(m)
           call plot_p(filename,BC_type(m)) 

            write(filename,"('boundaryc-'I3.3'.dat')") Bc_type(m)
            call plot_center(filename,BC_type(m)) 
    enddo
    
    	   
	   print*, "if you want to plot 2D plane flow ?  1 for yes, 0 for no"
	   read(*,*) If3d
       if(If3d .eq. 1) then
	    call Plot_2D
       endif
      
	 
	   print*, "if you want to plot 3D flow  (tecplot type) ?  1 for yes, 0 for no"
	   read(*,*) If3d
	   if(If3d .eq. 1) then
	    call Plot_3D
	   endif
  	 
	   print*, "if you want to plot 3D flow (PLOT3D type) ?  1 for yes, 0 for no"
	   read(*,*) If3d
	   if(If3d .eq. 1) then
	    call Plot_3D1
	   endif
  
	 end

!=====================================================================================
     subroutine plot_3D
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition
	 open(100,file="flow3d-tec.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone i=", B%nx, " j= ", B%ny , " k= ", B%nz  
      do k=1,B%nz
      do j=1,B%ny
      do i=1,B%nx
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
      enddo
     enddo
     close(100)
    end

     subroutine plot_3D1
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition
	 open(100,file="flow3d-plot3d.dat", form="unformatted")
     do m=1,Num_Block
     B=>Block(m)
	   write(100)   (((B%d(i,j,k), i=1, B%nx), j=1, B%ny), k=1, B%nz),   &
	                       (((B%u(i,j,k), i=1, B%nx), j=1, B%ny), k=1, B%nz),   &
	                       (((B%v(i,j,k), i=1, B%nx), j=1, B%ny), k=1, B%nz),   &
	                       (((B%w(i,j,k), i=1, B%nx), j=1, B%ny), k=1, B%nz),   &
	                       (((B%T(i,j,k), i=1, B%nx), j=1, B%ny), k=1, B%nz)     
	 
	 enddo
     close(100)
    end


     subroutine plot_2D
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,i0,j0,k0,Iflag

	 print*, "If you want to plot 2D section ? (0/1) "
	 read(*,*)  Iflag
	 if (Iflag ==0) return
	 print*, "Please input i0, j0,k0"
	 read(*,*) i0,j0,k0

     print*, "... flow2d-i.dat..."
     i=i0
	 open(100,file="flow2d-i.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone j=", B%ny , " k= ", B%nz  
      do k=1,B%nz
      do j=1,B%ny
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
     enddo
     close(100)

       print*, "... flow2d-j.dat..."
      
     j=j0
  	 open(100,file="flow2d-j.dat")
	 write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone i=", B%nx,  " k= ", B%nz  
      do k=1,B%nz
      do i=1,B%nx
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
     enddo
     close(100)

     print*, " ... flow2d-k.dat..."
	 k=k0
	 open(100,file="flow2d-k.dat")
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
      write(100,*) "zone i=", B%nx, " j= ", B%ny   
      do j=1,B%ny
      do i=1,B%nx
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
     enddo
     close(100)
  
    end




   
    subroutine plot_p(filename,bCondition)
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition
     character(len=50):: filename

	
	 open(100,file=trim(filename))
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     do m=1,Num_Block
     B=>Block(m)
     do n=1,B%subface
     Bc=>B%bc_msg(n)
     if(Bc%bc .eq. bCondition) then
      write(100,*) "zone i=", (Bc%ie-Bc%ib+1), " j= ",(Bc%je-Bc%jb+1) , " k= ", Bc%ke-Bc%kb+1  
      do k=Bc%kb,Bc%ke
      do j=Bc%jb,Bc%je
      do i=Bc%ib,Bc%ie
	   write(100,"(13E18.9)") B%x(i,j,k),B%y(i,j,k),B%z(i,j,k), &
	                          B%d(i,j,k),B%u(i,j,k),B%v(i,j,k),  B%w(i,j,k),B%T(i,j,k), &
							  B%p(i,j,k),B%mut1(i,j,k)*Re,B%vt1(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
      enddo
     endif
 
     enddo
     enddo
     close(100)
    end

    

 
   
!------------------------------------------------------------------------------------

   subroutine plot_center(filename,bCondition)
     use Global_Variables
     implicit none
     Type (Block_TYPE),pointer:: B
     TYPE(BC_MSG_TYPE),pointer::Bc
     integer::  i,j,k,m,n,bCondition,ib,ie,jb,je,kb,ke
      real(PRE_EC):: Cpw,Cf,Cfx,mu0,Tsb

	 character(len=50):: filename

     open(100,file=trim(filename))
     write(100,*) "variables=x,y,z,d,u,v,w,T,p,mut,vt,Kt,Wt"
     if(BCondition==BC_Wall) then
	 open(101,file="CpCf.dat")
	 write(101,*) "variables=x,y,z,Cp,Cf,Cfx,Pw"
	 endif
  
     Tsb=110.4d0/T_inf
	
	 do m=1,Num_Block
     B=>Block(m)
     do n=1,B%subface
     Bc=>B%bc_msg(n)
     if(Bc%bc .eq. bCondition) then
       ib=Bc%ib; ie=Bc%ie-1; jb=Bc%jb; je=Bc%je-1; kb=Bc%kb; ke=Bc%ke-1
       if(Bc%ib .eq. Bc%ie) then
	     if(Bc%ib .eq. 1) then
		  ie=1
		 else
		  ib=ib-1
		 endif
		endif

       if(Bc%jb .eq. Bc%je) then
	     if(Bc%jb .eq. 1) then
		  je=1
		 else
		  jb=jb-1
		 endif
		endif

       if(Bc%kb .eq. Bc%ke) then
	     if(Bc%kb .eq. 1) then
		  ke=1
		 else
		  kb=kb-1
		 endif
		endif

			    
	  write(100,*) "zone i=", (ie-ib+1), " j= ",(je-jb+1) , " k= ", ke-kb+1  
  	  do k=kb,ke
      do j=jb,je
      do i=ib,ie
       write(100,"(13E18.9)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),B%dc(i,j,k),B%uc(i,j,k),B%vc(i,j,k), &
                     B%wc(i,j,k),B%Tc(i,j,k),B%pc(i,j,k),B%mut(i,j,k),B%vt(i,j,k),B%Kt(i,j,k),B%Wt(i,j,k)
      enddo
      enddo
      enddo
      
	  if(bCondition==BC_Wall) then
 	  write(101,*) "zone i=", (ie-ib+1), " j= ",(je-jb+1) , " k= ", ke-kb+1  
  	  do k=kb,ke
      do j=jb,je
      do i=ib,ie
       mu0=1.d0/Re*(1.d0+Tsb)*sqrt(B%Tc(i,j,k)**3)/(Tsb+B%Tc(i,j,k))
       Cpw=(B%pc(i,j,k)-p00)*2.d0
	   Cf=mu0*(sqrt(B%uc(i,j,k)**2+B%vc(i,j,k)**2+B%wc(i,j,k)**2)/B%dw(i,j,k))*2.d0
	   Cfx=mu0*(B%uc(i,j,k)/B%dw(i,j,k))*2.d0
	   write(101,"(13E18.9)") B%xc(i,j,k),B%yc(i,j,k),B%zc(i,j,k),Cpw,Cf,Cfx,B%dc(i,j,k)*B%Tc(i,j,k)
      enddo
      enddo
      enddo
      endif

     endif

     enddo
     enddo
     close(100)
	 close(101)
    end
 



!------------------------------------------------------------------------------     
! Read the message of the mesh and the initial flow;
! ��ȡ���񣬳�ʼ������Ϣ; 
! �����ڴ������
! ���㼸������
!------------------------------------------------------------------------------
   subroutine init
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,km,nx,ny,nz,Num_Block1,ksub,tmp,NB1,Iflag
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Block_TYPE),pointer:: B
   TYPE(BC_MSG_TYPE),pointer::Bc
   logical:: EX
   integer:: ib,ie,jb,je,kb,ke
!------------------------------------------------------------------
!    print*, "Please input Ma, Re, T_inf, Mesh_File_Format"
!    read(*,*) Ma, Re, T_inf,Mesh_File_Format


    call read_parameter_ec 

    print*, "Ma,Re,gamma,Tinf=", Ma,Re,gamma,T_inf
    print*, "Mesh_File_Format=", Mesh_File_Format

    p00=1.d0/(gamma*Ma*Ma)

! ---------node Coordinates----------------------------------------  
!  �����ļ���PLOT3D��ʽ��   
   print*, "read Mesh3d.dat... (PLOT3D Format)"
   if(Mesh_File_Format .eq. 0) then
    open(99,file="Mesh3d.dat",form="unformatted")
    read(99) Num_Block         ! �ܿ���
   else
    open(99,file="Mesh3d.dat")
    read(99,*) Num_Block         ! �ܿ���
   endif
   print*, "Num_Block=", Num_Block
   
   
      
    allocate(Block(Num_Block))             
    allocate(NI(Num_Block),NJ(Num_Block),NK(Num_Block) )   ! ÿ��Ĵ�С
  if(Mesh_File_Format .eq. 0) then
   read(99) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  else
   read(99,*) (NI(k), NJ(k), NK(k), k=1,Num_Block)
  endif

! ��ȡÿ����Ϣ----------------------------------------   
    do m=1,Num_Block
     B => Block(m)
     B%nx=NI(m); B%ny=NJ(m) ; B%nz=NK(m)   ! nx,ny,nz ÿ��Ĵ�С
     nx=B%nx ; ny= B%ny ; nz=B%nz
! ----------  ������ -----------------------------------------------
    allocate(B%xc(0:nx,0:ny,0:nz), B%yc(0:nx,0:ny,0:nz), B%zc(0:nx,0:ny,0:nz)) 
    allocate(B%dc(0:nx,0:ny,0:nz),B%uc(0:nx,0:ny,0:nz),B%vc(0:nx,0:ny,0:nz), &
       B%wc(0:nx,0:ny,0:nz),B%Tc(0:nx,0:ny,0:nz),B%pc(0:nx,0:ny,0:nz), &
	   B%mut(0:nx,0:ny,0:nz),B%vt(0:nx,0:ny,0:nz),B%Kt(0:nx,0:ny,0:nz),B%Wt(0:nx,0:ny,0:nz))
    allocate(B%x(nx,ny,nz),B%y(nx,ny,nz),B%z(nx,ny,nz))
    allocate(B%d(nx,ny,nz),B%u(nx,ny,nz),B%v(nx,ny,nz), &
             B%w(nx,ny,nz),B%T(nx,ny,nz),B%p(nx,ny,nz),B%mut1(nx,ny,nz),B%vt1(nx,ny,nz))
    allocate(B%dw(nx-1,ny-1,nz-1))

   B%mut(:,:,:)=0.d0
   B%vt(:,:,:)=0.d0
   B%mut1(:,:,:)=0.d0
   B%Kt(:,:,:)=0.d0
   B%Wt(:,:,:)=0.d0
   B%vt1(:,:,:)=0.d0

   print*, "m=", m
   print*, "nx,ny,nz=", nx,ny,nz
   if(Mesh_File_Format .eq. 0) then
	read(99)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   else
  	read(99,*)   (((B%x(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%y(i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
               (((B%z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
   endif
   print*, "read mesh ok "
!----------��ֵ���������ĵ������-----------------------------------
   do k=1,B%nz-1
   do j=1,B%ny-1
   do i=1,B%nx-1
     B%xc(i,j,k)=(B%x(i,j,k)+B%x(i+1,j,k)+B%x(i,j+1,k)+        &
	              B%x(i,j,k+1)+B%x(i+1,j+1,k)+B%x(i+1,j,k+1) + &
                  B%x(i,j+1,k+1)+B%x(i+1,j+1,k+1))/8.d0
     B%yc(i,j,k)=(B%y(i,j,k)+B%y(i+1,j,k)+B%y(i,j+1,k)+        &
	              B%y(i,j,k+1)+B%y(i+1,j+1,k)+B%y(i+1,j,k+1) + &
                  B%y(i,j+1,k+1)+B%y(i+1,j+1,k+1))/8.d0
     B%zc(i,j,k)=(B%z(i,j,k)+B%z(i+1,j,k)+B%z(i,j+1,k)+        &
	              B%z(i,j,k+1)+B%z(i+1,j+1,k)+B%z(i+1,j,k+1) + &
                  B%z(i,j+1,k+1)+B%z(i+1,j+1,k+1))/8.d0
   
   enddo
   enddo
   enddo
  enddo
 
  close(99)

    inquire( file="wall_dist.dat", exist=EX)
    if(EX) then
      open(100,file="wall_dist.dat",form="unformatted")
      do m=1,Num_Block
      B => Block(m)
      read(100) (((B%dw(i,j,k),i=1,B%nx-1),j=1,B%ny-1),k=1,B%nz-1)
     enddo
     close(100)
   endif
   
!---------------��flow3d.dat�ж�ȡ����Ϊ��ֵ----------------------
  print*, "input 1 or 2,  1 read flow3d.dat,  2 read flow3d_average.dat "
  read(*,*) Iflag
  
  if(Iflag == 1) then   
	open(99,file="flow3d.dat",form="unformatted")
  else
	open(99,file="flow3d_average.dat",form="unformatted")
  endif
    
	 do m=1,Num_Block
     B => Block(m)
        read(99)   (((B%dc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%uc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%vc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
                   (((B%wc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) , &
		           (((B%Tc(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) 
          do k=0, B%nz
		  do j=0, B%ny
		  do i=0, B%nx
		    B%pc(i,j,k)=p00*B%dc(i,j,k)*B%Tc(i,j,k)
		  enddo
		  enddo
		  enddo

	 enddo
	close(99)
    
	 inquire(file="SA3d.dat",exist=EX)
	 if(EX) then
     open(100,file="SA3d.dat",form="unformatted")
      do m=1,Num_Block
        B => Block(m)
        read(100)   (((B%vt(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz) 
 	  enddo
      close(100)
     endif

    inquire(file="SST3d.dat",exist=EX)
	 if(EX) then
     open(100,file="SST3d.dat",form="unformatted")
     do m=1,Num_Block
     B => Block(m)
        read(100)   (((B%Kt(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz),  (((B%Wt(i,j,k),i=0,B%nx),j=0,B%ny),k=0,B%nz)
 	 enddo
 	 close(100)
    endif
!=================================================================================
   call comput_value_in_mesh


!---------------------------------------------------- 
   BC_number=0
   BC_type=0
   
!--------------------------------------------------
    print*, "read bc3d.inp"
    open(88,file="bc3d.inp")
    read(88,*)
    read(88,*) Num_Block1
    if(Num_Block1 .ne. Num_Block) then
      print*, "Error!  Block number in bc2d.in is not equal to that in Mesh3d.dat !"
      stop
    endif
    
    do m=1,Num_Block
     B => Block(m)
     read(88,*)
     read(88,*)
     read(88,*) B%subface   !number of the subface in the Block m
!-------------------------------------------------------------------------------
     allocate(B%bc_msg(B%subface))
     do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     read(88,*) ib,ie,jb,je,kb,ke,Bc%bc
	 if(Bc%bc .lt. 0) then
  	   read(88,*)
  	  endif
  	  
 ! Find boundary type 	  
   if(Bc%bc .ge. 0) then
      if(BC_number .eq. 0) then
       BC_number=1
       BC_type(1)=Bc%bc
      else
       km=1
       do k=1,BC_number
        if(Bc%bc .eq. BC_type(k)) km=0
       enddo
      if(km .eq. 1) then
       Bc_number=Bc_number+1
       Bc_type(Bc_number)=Bc%bc
      endif
     endif
    endif



	  Bc%ib=min(abs(ib),abs(ie)) ;  Bc%ie=max(abs(ib),abs(ie))
	  Bc%jb=min(abs(jb),abs(je)) ;  Bc%je=max(abs(jb),abs(je))
	  Bc%kb=min(abs(kb),abs(ke)) ;  Bc%ke=max(abs(kb),abs(ke))
     enddo
     enddo
     close(88)
	 print*, "read bc3d.in ok ..."


 end subroutine init  

   subroutine comput_value_in_mesh
   use  Global_Variables
   implicit none
   integer:: nx,ny,nz,i,j,k,m
   Type (Block_TYPE),pointer:: B
   do m=1, NUM_BLOCK
   B=> Block(m)
   nx=B%nx; ny=B%ny; nz=B%nz
!  �����Ĵ�����������ֵ�������   
   do k=1,B%nz
   do j=1,B%ny
   do i=1,B%nx
     B%d(i,j,k)=(B%dc(i-1,j-1,k-1)+B%dc(i,j-1,k-1)+B%dc(i-1,j,k-1)+B%dc(i,j,k-1)  &
	             +B%dc(i-1,j-1,k)+B%dc(i,j-1,k)+B%dc(i-1,j,k)+B%dc(i,j,k))/8.d0
     B%u(i,j,k)=(B%uc(i-1,j-1,k-1)+B%uc(i,j-1,k-1)+B%uc(i-1,j,k-1)+B%uc(i,j,k-1)  &
	             +B%uc(i-1,j-1,k)+B%uc(i,j-1,k)+B%uc(i-1,j,k)+B%uc(i,j,k))/8.d0
     B%v(i,j,k)=(B%vc(i-1,j-1,k-1)+B%vc(i,j-1,k-1)+B%vc(i-1,j,k-1)+B%vc(i,j,k-1)  &
	             +B%vc(i-1,j-1,k)+B%vc(i,j-1,k)+B%vc(i-1,j,k)+B%vc(i,j,k))/8.d0
     B%w(i,j,k)=(B%wc(i-1,j-1,k-1)+B%wc(i,j-1,k-1)+B%wc(i-1,j,k-1)+B%wc(i,j,k-1)  &
	             +B%wc(i-1,j-1,k)+B%wc(i,j-1,k)+B%wc(i-1,j,k)+B%wc(i,j,k))/8.d0

     B%T(i,j,k)=(B%Tc(i-1,j-1,k-1)+B%Tc(i,j-1,k-1)+B%Tc(i-1,j,k-1)+B%Tc(i,j,k-1)  &
	             +B%Tc(i-1,j-1,k)+B%Tc(i,j-1,k)+B%Tc(i-1,j,k)+B%Tc(i,j,k))/8.d0

     B%p(i,j,k)=(B%pc(i-1,j-1,k-1)+B%pc(i,j-1,k-1)+B%pc(i-1,j,k-1)+B%pc(i,j,k-1)  &
	             +B%pc(i-1,j-1,k)+B%pc(i,j-1,k)+B%pc(i-1,j,k)+B%pc(i,j,k))/8.d0

     B%mut1(i,j,k)=(B%mut(i-1,j-1,k-1)+B%mut(i,j-1,k-1)+B%mut(i-1,j,k-1)+B%mut(i,j,k-1)  &
	             +B%mut(i-1,j-1,k)+B%mut(i,j-1,k)+B%mut(i-1,j,k)+B%mut(i,j,k))/8.d0

     B%vt1(i,j,k)=(B%vt(i-1,j-1,k-1)+B%vt(i,j-1,k-1)+B%vt(i-1,j,k-1)+B%vt(i,j,k-1)  &
	             +B%vt(i-1,j-1,k)+B%vt(i,j-1,k)+B%vt(i-1,j,k)+B%vt(i,j,k))/8.d0

   enddo
   enddo
   enddo
   enddo
   end



!------read parameter (Namelist type)---------------- 
  subroutine read_parameter_ec 
   use Global_Variables
   implicit none
   real(PRE_EC),parameter:: PI=3.14159265358979d0
   real(PRE_EC):: R0, a0, d0, mu0,mu1
   logical:: has_legacy, has_freestream, has_flow
   integer:: ios
!---- Extra declarations so that ANY case control.ec can be parsed (legacy
!---- single group or the new grouped namelists).  These names are not used by
!---- this tool, but they MUST be declared: a namelist group containing a name
!---- the namelist does not know is a read error.
   real(PRE_EC):: LS_rho, LS_mu, LS_k, LS_Cp, LS_T_ref, LS_U_in, LS_V_in, LS_W_in, &
       LS_Mdot_in, LS_P_in, LS_P_out, LS_T_wall, LS_U_lid, LS_alpha_p, LS_alpha_u, &
       LS_alpha_T, LS_Tol, AC_beta, AC_CFL, AC_CFLv, AC_Tol, AC_w, AC_WenoBlend, &
       AC_MomFrac, Porous_T_ref, Porous_alpha_Ts, Porous_Tol, Porous_U_in, &
       Porous_V_in, Porous_W_in, Porous_T_in, Twall_Couple_Init, Tol_Couple_Tw, &
       Tol_Couple_p, Tol_Couple_u, Solid_GS_Omega, Solid_Tol
   integer:: LS_Inlet_Type, LS_Max_Iter, LS_Scheme, LS_Algorithm, AC_Max_Iter, &
       AC_Print, AC_Flux, AC_Recon, AC_Limiter, AC_WallRecon, AC_WallP, AC_MomDiss, &
       Porous_Max_Iter, Iflag_Couple_Scheme, Kstep_Couple_Comp, Niter_Couple_Outer, &
       Porous_Chunk_Iter, Niter_Couple_Warm, Kstep_Couple_Min, Iflag_Couple_WallFlux, &
       Solid_Max_Iter, Solid_Min_Iter
   
   real(PRE_EC)::  AoA, AoS, p_outlet, t_end, Lscale, &
	    PrL, PrT, &
        dt_global,CFL,dtmax,dtmin,Time_Method, &
		w_LU, Twall,Kt_inf,Wt_inf,&
		Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, &
        Ref_S,Ref_L,Centroid(3), &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
 		Turbo_P0,Turbo_T0, Turbo_L0,Turbo_w, Turbo_Periodic_seta, &
		Periodic_dX, Periodic_dY, Periodic_dZ
  
   integer:: Kstep_save, Kstep_average, Iflag_turbulence_model,Iflag_init,If_viscous,  &
           Iflag_local_dt,  If_Residual_smoothing, If_dtime_mesh,  &
		   Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction, Pre_Step_Mesh(3), &
           Kstep_show,Kstep_smooth,Kstep_init_smooth, Bound_Scheme,  &
           Num_Mesh,IF_Debug, NUM_THREADS, Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug(4),  &
           IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, IF_Innerflow,Iflag_savefile, &
		   Iflag_vtk_onefile, Iflag_vtk_SI, Iflag_bc_check

!---- grouped namelists (new-style control.ec).  Only the two groups this tool
!---- needs are declared, but ALL of their variables must be listed: reading a
!---- group that contains a name the namelist does not know is an error.
	namelist /freestream_ec/ Ma, Re, AoA, AoS, p_outlet, gamma, PrL, PrT, &
	    T_inf, Twall, Lscale, Ref_S, Ref_L, Centroid, Cood_Y_UP, Kt_inf, Wt_inf, &
	    IF_TurboMachinary, Ref_medium_usrdef, &
	    Turbo_P0, Turbo_T0, Turbo_L0, Turbo_w, Turbo_Periodic_seta

	namelist /flow_ec/ t_end, Kstep_save, Kstep_show, Kstep_average, &
	    Kstep_smooth, Kstep_init_smooth, CFL, dt_global, dtmax, dtmin, &
	    Iflag_local_dt, Time_Method, Step_Inner_Limit, Res_Inner_Limit, w_LU, &
	    If_Residual_smoothing, If_dtime_mesh, &
	    Iflag_Scheme, Iflag_Flux, IFlag_Reconstruction, Bound_Scheme, &
	    IF_Scheme_Positivity, IFLAG_LIMIT_FLOW, Iflag_turbulence_model, &
	    If_viscous, Iflag_init, MUT_MAX, CP1_NSA, CP2_NSA, &
	    Ldmin, Ldmax, Lpmin, Lpmax, Lumax, LSAmax, &
	    Mesh_File_Format, Num_Mesh, Pre_Step_Mesh, NUM_THREADS, IF_Debug, Pdebug, &
	    Periodic_dX, Periodic_dY, Periodic_dZ, IF_Innerflow, &
	    Iflag_savefile, Iflag_vtk_onefile, Iflag_vtk_SI, Iflag_bc_check

!---- legacy single group (all existing case files)
	namelist /control_ec/ Ma, Re, AoA, AoS, p_outlet, t_end, &
	    gamma, PrL, PrT, &
	    Kstep_save, &
	    Iflag_turbulence_model,Iflag_init,If_viscous,  &
        Iflag_local_dt,dt_global,CFL,dtmax,dtmin,Time_Method, &
		If_Residual_smoothing,w_LU,If_dtime_mesh,  &
        Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction, &
		Mesh_File_Format,Kstep_show,Kstep_average,Kstep_smooth,Kstep_init_smooth,  &
        Num_Mesh,T_inf,Twall,Kt_inf,Wt_inf,IF_Debug,NUM_THREADS,  &
		Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, Bound_Scheme, &
        Pre_Step_Mesh,Ref_S,Ref_L,Centroid,Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug, &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
        IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, &
		Turbo_P0,Turbo_T0, Turbo_L0,Turbo_w, Turbo_Periodic_seta, &
		Periodic_dX, Periodic_dY, Periodic_dZ, &
		IF_Innerflow, Iflag_savefile, &
		Lscale, LS_rho, LS_mu, LS_k, LS_Cp, LS_T_ref, LS_Inlet_Type, &
		LS_U_in, LS_V_in, LS_W_in, LS_Mdot_in, LS_P_in, LS_P_out, &
		LS_T_wall, LS_U_lid, LS_alpha_p, LS_alpha_u, LS_alpha_T, &
		LS_Max_Iter, LS_Tol, LS_Scheme, LS_Algorithm, &
		AC_Max_Iter, AC_Print, AC_beta, AC_CFL, AC_CFLv, AC_Tol, AC_w, &
		AC_Flux, AC_Recon, AC_Limiter, AC_WenoBlend, AC_WallRecon, AC_WallP, &
		AC_MomDiss, AC_MomFrac, &
		Porous_T_ref, Porous_alpha_Ts, Porous_Max_Iter, Porous_Tol, &
		Porous_U_in, Porous_V_in, Porous_W_in, Porous_T_in, &
		Solid_GS_Omega, Solid_Max_Iter, Solid_Min_Iter, Solid_Tol, &
		Iflag_Couple_Scheme, Kstep_Couple_Comp, Niter_Couple_Outer, &
		Porous_Chunk_Iter, Niter_Couple_Warm, Kstep_Couple_Min, &
		Twall_Couple_Init, Tol_Couple_Tw, Tol_Couple_p, Tol_Couple_u, &
		Iflag_Couple_WallFlux


	

!---- default--------
        Ma=1.d0
        Re=100.d0
    	gamma=1.4d0
        T_inf=288.15
        Mesh_File_Format=0
        
!---------------------------------
	open(99,file="control.ec")
    call util_scan_groups(has_legacy, has_freestream, has_flow)
    print*, ' readflow3d: control.ec groups found: legacy=', has_legacy, &
            ' freestream_ec=', has_freestream, ' flow_ec=', has_flow
    if(has_legacy) then
      rewind(99); read(99,nml=control_ec,iostat=ios)
      if(ios /= 0) then
        print*, ' ERROR: cannot read namelist $control_ec in control.ec, iostat=', ios
        stop 1
      endif
    endif
    if(has_freestream) then
      rewind(99); read(99,nml=freestream_ec,iostat=ios)
      if(ios /= 0) then
        print*, ' ERROR: cannot read namelist $freestream_ec in control.ec, iostat=', ios
        stop 1
      endif
    endif
    if(has_flow) then
      rewind(99); read(99,nml=flow_ec,iostat=ios)
      if(ios /= 0) then
        print*, ' ERROR: cannot read namelist $flow_ec in control.ec, iostat=', ios
        stop 1
      endif
    endif
    close(99)
 
 !---- convert parameters ----------------------
 ! Ref_medium_usrdef==0 ʹ��Ĭ�Ͻ��� (Ma=1, �������¡���ѹ���� Re) �� ==1 ʹ���Զ������ ����Ϊ����Ma, Re�ȣ�
    if( (IF_TurboMachinary ==1 .or. IF_Innerflow ==1) .and.  Ref_medium_usrdef == 0) then   ! Ĭ�Ͽ������ʣ�����Mach���� Reynolds��
      
	  T_inf=Turbo_T0  ! �ο��¶� ���������£�
      gamma=1.4d0    ! 
	  PrL=0.7d0   ! Prandtl��
	  PrT=0.9d0
      R0= 287.06d0   ! ���������峣��R
	  a0= sqrt(gamma*R0*Turbo_T0)    ! �ο��¶��µ����� 
	  mu0=1.179d-5     ! ����ճ��ϵ�� (288.15K)  
      mu1=mu0* sqrt((Turbo_T0/288.15d0)**3)*(288.15d0+110.4d0)/(Turbo_T0+110.4d0)  ! �ο��¶��µĿ���ճ��ϵ��
      d0=Turbo_P0/(R0*Turbo_T0)
	  Re=d0*a0*Turbo_L0/mu1    ! �ο��¶��£��������˶���Reynolds��
	  Ma=1.d0     ! Mach��    ����������Ϊ�ο��ٶȣ�����ο�Mach��Ϊ1��
      Turbo_w= 2.d0*PI*Turbo_w/(a0/Turbo_L0)   ! �����ٽ��ٶ� Turbo_W��ת/�룩
      P_outlet=P_outlet/Turbo_P0    ! ��ѹ �������٣�
	endif

 !  print*, "Re, Ma, gamma, T_inf=", Re, Ma, gamma, T_inf
 !  print*, "Mesh_File_Format=", Mesh_File_Format


 !----------------------------------------------
 end
!--------------------------------------------------
 



!==============================================================================
! Lower-case a string in place (used by the control.ec group scan below).
!==============================================================================
  subroutine util_lower(s)
   implicit none
   character(len=*),intent(inout):: s
   integer:: i, ic
   do i=1,len(s)
     ic=iachar(s(i:i))
     if(ic >= 65 .and. ic <= 90) s(i:i)=achar(ic+32)   ! 'A'..'Z' -> 'a'..'z'
   enddo
  end subroutine util_lower

!==============================================================================
! Report which control.ec namelist groups are present: the historical single
! group "$control_ec" and/or the new grouped namelists "$freestream_ec" and
! "$flow_ec" (the only two this post-processing tool needs).  Text after '!' is
! ignored; both '$' and '&' delimiters are accepted; the scan is case
! insensitive.  Reading an absent namelist group is a runtime error in Fortran,
! hence the pre-scan.
!==============================================================================
  subroutine util_scan_groups(has_legacy, has_freestream, has_flow)
   implicit none
   logical,intent(out):: has_legacy, has_freestream, has_flow
   character(len=512):: line
   integer:: ios, k, nz
   has_legacy=.false.; has_freestream=.false.; has_flow=.false.
   open(96,file="control.ec",status='old')
   do
     read(96,'(A)',iostat=ios) line
     if(ios /= 0) exit
     k=index(line,'!')
     if(k > 0) line(k:)=' '
     call util_lower(line)
!    only a real group header counts (first non-blank char must be '$' or '&')
     nz=1
     do
       if(nz > len(line)) exit
       if(line(nz:nz) == ' ' .or. iachar(line(nz:nz)) == 9) then
         nz=nz+1
       else
         exit
       endif
     enddo
     if(nz > len(line)) cycle
     if(line(nz:nz) /= '$' .and. line(nz:nz) /= '&') cycle
     if(index(line,'control_ec') > 0)    has_legacy=.true.
     if(index(line,'freestream_ec') > 0) has_freestream=.true.
     if(index(line,'flow_ec') > 0)       has_flow=.true.
   enddo
   close(96)
  end subroutine util_scan_groups

