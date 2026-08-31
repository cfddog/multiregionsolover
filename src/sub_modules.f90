 !
!-----------------------------------------------------------------------------------
! Consts  
  module precision_EC
    use mpi
  implicit none
  !include "mpif.h"
!    integer,parameter:: PRE_EC=4            ! Single precision
     integer,parameter:: PRE_EC=8           ! Double Precision
!    integer,parameter:: OCFD_DATA_TYPE=MPI_REAL  
     integer,parameter:: OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION ! Double precision
  end module  precision_EC 

!------------------------------------------------------------------------------------

 !
  module const_var
   use precision_EC
   implicit none
   real(PRE_EC),parameter::  PI=3.1415926535897932d0, Lim_Zero=1.d-20 !
   integer,parameter:: LAP=4 !
   integer,parameter:: Scheme_UD1=0, Scheme_NND2=1, Scheme_UD3=2,Scheme_MUSCL2U=3,Scheme_MUSCL2C=4,   &
                       Scheme_MUSCL3=5,Scheme_OMUSCL2=6,Scheme_WENO5=7,Scheme_UD5=8, Scheme_WENO7=9
   integer,parameter:: Scheme_CD2=20, Scheme_none=-1 !
   integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_Van_Leer=5,Flux_Ausm=6
   integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2
   integer,parameter:: BLOCK_FLUID=0, BLOCK_SOLID=1, BLOCK_LOWSPEED=2, BLOCK_POROUS=3


 !
!  integer,parameter:: BC_Wall=-10, BC_Farfield=-20, BC_Periodic=-30,BC_Symmetry=-40,BC_Outlet=-22
   integer,parameter::  BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6 
   integer,parameter::  BC_Wall_Turbo=201 !
   integer,parameter::  BC_Periodic=501, BC_Extrapolate=401 !
 !

 !

 !
   integer,parameter::  BC_Inner=-1, BC_PeriodicL=-2, BC_PeriodicR=-3   
   
   integer,parameter::  BC_Zero=0 !

 !
   integer,parameter::  BC_USER_FixedInlet=901, BC_USER_Inlet_time=902 !
   integer,parameter:: BC_USER_Blow_Suction_Wall=903 !

   integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
   integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2, Turbulence_SST=3,Turbulence_NewSA=21
   integer,parameter:: Init_continue=1, Init_By_FreeStream=0, Init_By_Zeroflow=-1,  Smooth_2nd=0,Smooth_4th=1 
!   real(PRE_EC), parameter::  Density_LIMIT=1.d-4,Temperature_LIMIT=1.d-4,Pressure_LIMIT=1.d-4

   integer,parameter:: Method_FVM=0, Method_FDM=1 !
   integer,parameter:: FD_WENO5=1,FD_WENO7=2,FD_OMP6=3 !
   integer,parameter:: FD_Steger_Warming=1,FD_Van_Leer=2


  end module const_var

 !
 module Mod_Type_Def
   use precision_EC
   implicit none

    TYPE BC_MSG_TYPE !
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no !
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1 !
	 integer:: L1,L2,L3 !
   END TYPE BC_MSG_TYPE

 !
   TYPE Block_TYPE !
     integer::  Block_no,mpi_id !
	 integer::  nx,ny,nz !
	 integer::  subface !
 !
     real(PRE_EC),pointer,dimension(:,:,:):: x,y,z !
     real(PRE_EC),pointer,dimension(:,:,:):: xc,yc,zc !
     real(PRE_EC),pointer,dimension(:,:,:):: Vol,Si,Sj,Sk !
	 real(PRE_EC),pointer,dimension(:,:,:):: ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3 !
 !
	 real(PRE_EC),pointer,dimension(:,:,:):: ix1,iy1,iz1,jx1,jy1,jz1,kx1,ky1,kz1
	 real(PRE_EC),pointer,dimension(:,:,:):: ix2,iy2,iz2,jx2,jy2,jz2,kx2,ky2,kz2
	 real(PRE_EC),pointer,dimension(:,:,:):: ix3,iy3,iz3,jx3,jy3,jz3,kx3,ky3,kz3
	 real(PRE_EC),pointer,dimension(:,:,:):: ix0,iy0,iz0,jx0,jy0,jz0,kx0,ky0,kz0
 !
	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U,Un,Un1 !
     real(PRE_EC),pointer,dimension(:,:,:,:) :: Res !
     real(PRE_EC),pointer,dimension(:,:,:):: dt !
     real(PRE_EC),pointer,dimension(:,:,:,:) :: QF !
     real(PRE_EC),pointer,dimension(:,:,:,:) :: deltU !
     real(PRE_EC),pointer,dimension(:,:,:,:) :: DU !
	 real(PRE_EC),pointer,dimension(:,:,:):: dw         ! turbulent viscous ; distance to the wall  (used in SA model)
	 real(PRE_EC),pointer,dimension(:,:,:):: surf1,surf2,surf3,surf4,surf5,surf6 !
	 real(PRE_EC),pointer,dimension(:,:,:):: mu,mu_t !
     real(PRE_EC),pointer,dimension(:,:,:):: dtime_mesh !

	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U_average !
	 
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg !
     integer,pointer,dimension(:,:,:):: BcI,BcJ,BcK !
   	 
	 integer:: IFLAG_FVM_FDM !
	 integer:: IF_OverLimit !
	integer:: Block_type                            ! block attribute: 0=fluid, 1=solid, 2=low-speed, 3=porous
	 real(PRE_EC),pointer,dimension(:,:,:):: Ts,Tsn  ! solid temperature (current/previous time step, cell-centered, LAP ghost layers)
	 TYPE(BC_MSG_TYPE),pointer,dimension(:):: bc_msg2  ! interface connection info (from bc3d_interface.inp)
	 real(PRE_EC):: solid_rho, solid_Cp, solid_k     ! solid material properties (density, specific heat, thermal conductivity)
	End TYPE Block_TYPE  

 !
 !
  
   TYPE Mesh_TYPE !
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep !
     integer:: NVAR !
	 real(PRE_EC)::  tt !
	 real(PRE_EC),pointer,dimension(:)::  Res_max,Res_rms !
	 TYPE (Block_TYPE),pointer,dimension(:):: Block !

 !
 !
 !
    integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction, Bound_Scheme
   End TYPE Mesh_TYPE
  
  end module Mod_Type_Def


!-------------------------------------------------------------------------------------- 
! Global Variables:
! Ma: Mach number ; Re: Reynolds number; gamma: Specific rato (=Cp/Cv); Pr: Prandtl number; 
! AoA: Angle of Attack; p00=1/(gamma*Ma*Ma), p=p00*d*T; 
! t_end: end time; 
! Num_Block: Total Block number 
! Kstep_save: Save data every Kstep_save step
! Iflag_turbulence_model: 0 no model, 1 BL model
! Iflag_flux: type of Splitting (or Riemann solver) : 0 Steger-Warming 1 HLL 2 HLLC 3 Roe
! Iflag_local_dt:  0 global time step, 1 local time step
! Iflag_Reconstruction: Schemes 1 NND 2 3rd Upwind 3 WENO3  4 MUSCL 
!----------------------------------------------------------------------------------------
 !
!========================================================================================

  module Global_Var    
   use const_var !
   use mod_type_def !
   implicit none


!---------------------------------------------------------------------------------------------
 !
!----------------------------------------------------------------------------

   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh !
   integer,save:: Num_Mesh,NVAR , Total_block, Num_block !
   integer,pointer,dimension(:):: bNi,bNj,bNk !
   integer,save::  Kstep_save, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method, &
	  Kstep_show,If_viscous,If_Residual_smoothing,Mesh_File_Format,IF_Debug, &
	  Kstep_smooth,Kstep_init_smooth,NUM_THREADS,If_dtime_mesh, Step_Inner_Limit, &  
      Bound_Scheme, & !
      IF_Walldist, IFLAG_LIMIT_FLOW, & !
      IF_Scheme_Positivity,          & !
      Kstep_average,                 & !
      Iflag_savefile !
   integer,save:: IF_TurboMachinary , Ref_medium_usrdef !
   integer,save:: IF_InnerFlow !

   integer,save:: FD_Flux,FD_scheme !
   integer,save:: KRK=0 !
   integer,save:: Istep_average=0 !
 !

   real(PRE_EC),save:: Ma,Re,gamma,Cp,Cv,t_end,P_OUTLET,&
                       A_alfa,A_beta,PrL,PrT,T_inf,Twall,w_LU,Kt_inf,Wt_inf , Res_Inner_Limit, MUT_MAX,AoA,Aos
   real(PRE_EC),save:: Turbo_Periodic_seta, Turbo_w, Turbo_P0, Turbo_T0 , Turbo_L0 !
   real(PRE_EC),save:: Periodic_dX,Periodic_dY,Periodic_dZ
   integer,pointer,dimension(:):: Block_Type_List ! global block type list (0=fluid, 1=solid, 2=low-speed, 3=porous)
   real(PRE_EC),pointer,dimension(:):: solid_rho_list, solid_Cp_list, solid_k_list  ! global solid material property lists

 
 !
   real(PRE_EC),save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin !
   integer,save:: Pre_Step_Mesh(3) !
   integer,save:: Cood_Y_UP !
   integer,save:: Pdebug(4) !
   real(PRE_EC),save:: Ref_S,Ref_L, Centroid(3) !
   real(PRE_EC),save:: Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax !
   real(PRE_EC),save:: CP1_NSA,CP2_NSA       ! parameters in New SA model
 !-----------mpi data ----------------------------------------------------------- 
   integer:: my_id,Total_proc !
   integer,pointer,dimension(:):: B_Proc, B_n !
   integer,pointer,dimension(:):: my_Blocks !
  end module Global_Var  
!----------------------------------------------------------------------------



 !
!------------SEC part---------------------------------------------------------------------

   module FDM_data
    use precision_EC
    implicit none
    TYPE FDM_Block_TYPE !
    real(PRE_EC), pointer,dimension(:,:,:):: ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac !
    End TYPE FDM_Block_TYPE  
  
   TYPE FDM_Mesh_TYPE !
	 TYPE (FDM_Block_TYPE),pointer,dimension(:):: Block !
   End TYPE FDM_Mesh_TYPE 
  
   TYPE (FDM_Mesh_TYPE),pointer,dimension(:):: FDM_Mesh !
  
   end  module FDM_data
