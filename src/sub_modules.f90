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

!  Explicit interface (multi-solver conjugate) boundary types.  Write these
!  codes directly in bc3d.inp / bc3d_interface.inp instead of the generic
!  negative code BC_Inner.  The code identifies the pair of block types and
!  therefore selects the matching coupling routine:
!    11 compressible (BLOCK_FLUID)  - solid (BLOCK_SOLID)
!    12 compressible (BLOCK_FLUID)  - low-speed (BLOCK_LOWSPEED)
!    13 low-speed  (BLOCK_LOWSPEED) - solid (BLOCK_SOLID)
!    14 solid      (BLOCK_SOLID)    - solid (BLOCK_SOLID)
!    15 compressible block-to-block interface (handled by buffer exchange)
   integer,parameter::  BC_Interface_FluidSolid = 11
   integer,parameter::  BC_Interface_FluidLow    = 12
   integer,parameter::  BC_Interface_LowSolid    = 13
   integer,parameter::  BC_Interface_SolidSolid  = 14
   integer,parameter::  BC_Interface_FluidFluid  = 15
!  Porous-media interfaces.  Same convention as 11-15: the code identifies the
!  PAIR of block types, never the direction.  On a given face the neighbour's
!  block type (Block_Type_List(nb1)) decides which routine handles the pair and
!  from which side it is processed once.
!    16 low-speed (BLOCK_LOWSPEED) - porous (BLOCK_POROUS)
!    17 solid     (BLOCK_SOLID)    - porous (BLOCK_POROUS)
!    18 porous    (BLOCK_POROUS)   - porous (BLOCK_POROUS)
!    19 compressible (BLOCK_FLUID) - porous (BLOCK_POROUS)
   integer,parameter::  BC_Interface_LowPorous    = 16
   integer,parameter::  BC_Interface_SolidPorous  = 17
   integer,parameter::  BC_Interface_PorousPorous = 18
   integer,parameter::  BC_Interface_FluidPorous  = 19
   integer,parameter::  BC_Interface_First = 11, BC_Interface_Last = 19

!  Explicit low-speed inlet / outlet types (kept distinct from the
!  compressible BC_Inflow / BC_Outflow codes).
   integer,parameter::  BC_LS_Inlet = 21
   integer,parameter::  BC_LS_Outlet = 22



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



   contains
!    True for every "inner" boundary: the legacy negative codes
!    (BC_Inner / BC_PeriodicL / BC_PeriodicR) and the explicit interface
!    types BC_Interface_FluidSolid .. BC_Interface_FluidFluid.
     logical function is_interface_bc(bc)
       integer,intent(in):: bc
       is_interface_bc = (bc < 0) .or. &
            (bc >= BC_Interface_First .and. bc <= BC_Interface_Last)
     end function is_interface_bc

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
	 real(PRE_EC),pointer,dimension(:,:,:):: p       ! pressure (low-speed solver, cell-centered, LAP ghost layers)
	 TYPE(BC_MSG_TYPE),pointer,dimension(:):: bc_msg2  ! interface connection info (from bc3d_interface.inp)
	 real(PRE_EC):: solid_rho, solid_Cp, solid_k     ! solid material properties (density, specific heat, thermal conductivity)
!     Solid thermal boundary condition data (from solid_bc.inp)
	 integer:: solid_bc_nface                       ! number of physical faces with thermal BCs
	 integer,pointer,dimension(:):: solid_bc_face_no ! face number in bc_msg (1-based)
	 real(PRE_EC),pointer,dimension(:):: solid_bc_Tw ! wall temperature (K, dimensional; >0 isothermal, <0 heat flux)
	 real(PRE_EC),pointer,dimension(:):: solid_bc_Qw ! heat flux (W/m2, dimensional; used when Tw<0)
	 real(PRE_EC),pointer,dimension(:):: solid_bc_htc  ! convective htc (W/m2/K); Robin BC when Tw==0 and htc>0
	 real(PRE_EC),pointer,dimension(:):: solid_bc_Tinf ! ambient temperature (K) for convective (Robin) BC
!     Porous-media material properties (block-uniform in the first version).
!     The solid skeleton uses solid_rho / solid_Cp / solid_k above (read from
!     material.in for every block); the extra porous quantities come from
!     porous.inp and are copied here by set_block_type.
	 real(PRE_EC):: porous_eps, porous_dp, porous_hv ! porosity, particle diam. (m), volumetric h_fs*A_fs (W/m3/K)
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
      Iflag_savefile, Iflag_vtk_onefile, Iflag_vtk_SI, &
      Iflag_bc_check !
   integer,save:: IF_TurboMachinary , Ref_medium_usrdef !
   integer,save:: IF_InnerFlow !

   integer,save:: LS_Inlet_Type, LS_Max_Iter    ! low-speed solver: inlet type (1=velocity,2=mass flow,3=pressure), SIMPLE inner iterations
   integer,save:: LS_Scheme=1    ! low-speed convection scheme: 1=1st-order upwind, 2=2nd-order upwind, 3=MUSCL(Van Leer)
   integer,save:: LS_Algorithm=1 ! low-speed pressure-velocity coupling: 1=SIMPLE, 2=SIMPLEC, 3=AC-FV (artificial compressibility)

!---- Artificial-compressibility (AC) low-speed solver controls (LS_Algorithm=3)
   integer,save:: AC_Max_Iter=40000  ! AC pseudo-time iterations per solver call
   integer,save:: AC_Print=500       ! print residual every AC_Print iterations
   real(PRE_EC),save:: AC_beta=10.d0 ! artificial compressibility: beta = AC_beta*U_ref^2
   real(PRE_EC),save:: AC_CFL=2.d0   ! inviscid CFL for the pseudo time step
   real(PRE_EC),save:: AC_CFLv=0.5d0 ! viscous CFL limit for the pseudo time step
   real(PRE_EC),save:: AC_Tol=1.d-7  ! convergence: dimensionless res_q & res_m < AC_Tol
   real(PRE_EC),save:: AC_w=1.d0     ! LU-SGS relaxation (>=0.5)
   integer,save:: AC_Flux=1   ! AC inviscid flux: 1=Rusanov(LLF), 2=Steger-Warming, 3=AUSM+
   integer,save:: AC_Recon=1  ! AC face reconstruction: 1=2nd-order MUSCL(van Leer), 2=WENO5, 3=WENO3
   integer,save:: AC_Limiter=1 ! MUSCL limiter: 1=van Leer, 2=minmod (more robust/dissipative)
   real(PRE_EC),save:: AC_WenoBlend=0.d0 ! WENO only: blend fraction of the 1st-order Rusanov flux (0..1)

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
   real(PRE_EC),pointer,dimension(:):: porous_eps_list, porous_dp_list, porous_hv_list ! global porous material lists

 
 !
   real(PRE_EC),save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin !
   integer,save:: Pre_Step_Mesh(3) !
   integer,save:: Cood_Y_UP !
   integer,save:: Pdebug(4) !
   real(PRE_EC),save:: Ref_S,Ref_L, Centroid(3) !
   real(PRE_EC),save:: Lscale ! Length scale: 1=m, 0.001=mm, etc.
   real(PRE_EC),save:: Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax !
   real(PRE_EC),save:: CP1_NSA,CP2_NSA       ! parameters in New SA model
   real(PRE_EC),save:: LS_rho, LS_mu, LS_k, LS_Cp, LS_T_ref    ! low-speed fluid properties (SI units)
   real(PRE_EC),save:: LS_U_in, LS_V_in, LS_W_in, LS_Mdot_in, LS_P_in, LS_P_out, LS_T_wall, LS_U_lid ! low-speed BC parameters
   real(PRE_EC),save:: LS_alpha_p, LS_alpha_u, LS_alpha_T, LS_Tol  ! low-speed under-relaxation + tolerance
   real(PRE_EC),save:: Porous_T_ref, Porous_alpha_Ts, Porous_Tol    ! porous: initial/frame T ref, Ts relaxation, SIMPLE tolerance
   real(PRE_EC),save:: Porous_U_in, Porous_V_in, Porous_W_in, Porous_T_in  ! porous coolant inlet velocity [m/s] + T [K]
   integer,save:: Porous_Max_Iter                                   ! porous SIMPLE inner iterations per solver call
 !---- interface-19 staggered (segmented) coupling controls --------------------
 ! Iflag_Couple_Scheme = 0 : per-step simultaneous coupling (default)
 !                       1 : staggered segmented coupling for FLUID<->POROUS (19):
 !                           gas chunk (compressible only, code-19 face = isothermal
 !                           wall at fp_Tw) -> extract q_w -> porous chunk (hot-end
 !                           heat flux q_w, coolant supply below) -> return fp_Tw.
   integer,save:: Iflag_Couple_Scheme=0
   integer,save:: Kstep_Couple_Comp=1000    ! compressible steps per gas chunk
   integer,save:: Niter_Couple_Outer=60     ! max outer staggered iterations
   integer,save:: Porous_Chunk_Iter=10      ! porous SIMPLE solver calls per porous chunk
   integer,save:: Niter_Couple_Warm=0       ! outer iters with the FULL gas chunk (warm start)
   integer,save:: Kstep_Couple_Min=1        ! gas-chunk step floor after warm-up (auto shrink)
   real(PRE_EC),save:: Twall_Couple_Init=300.d0  ! first gas-chunk interface wall T [K]
   real(PRE_EC),save:: Tol_Couple_Tw=2.d-2       ! outer convergence: max|dT_w| [K]
   real(PRE_EC),save:: Tol_Couple_p=1.d1         ! 12 (fluid-fluid) outer conv.: max|dp_w| [Pa]
   real(PRE_EC),save:: Tol_Couple_u=1.d-2        ! 12 (fluid-fluid) outer conv.: max|du_w| [m/s]
   integer,save:: Iflag_Couple_WallFlux=0  ! CHT split: 0=couple-based, 1=isothermal-Tw/qw wall-flux
 ! per-interface-face arrays used by the staggered coupling (face-cell indexed)
   real(PRE_EC),save,allocatable,dimension(:,:):: fp_Tw, fp_Tw_old, fp_qw, fp_pw
   real(PRE_EC),save,allocatable,dimension(:,:):: fp_u
   real(PRE_EC),save,allocatable,dimension(:,:):: fp_pw_old, fp_u_old
 ! solid GS q_w override descriptor (wall-flux CHT split mode)
   integer,save:: stg_ow_block=0
   integer,save:: stg_ow_face=0
   integer,save:: stg_ow_ib=0, stg_ow_ie=0, stg_ow_jb=0, stg_ow_je=0, stg_ow_kb=0, stg_ow_ke=0
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
