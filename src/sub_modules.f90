! ȫ�ֱ������ඨ��
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

! ����
  module const_var
   use precision_EC
   implicit none
   real(PRE_EC),parameter::  PI=3.1415926535897932d0, Lim_Zero=1.d-20   ! С�ڸ�ֵ��Ϊ0
   integer,parameter:: LAP=4                          ! ���������Ŀ ��ʹ��3�׸�ʽ��ֵ��С��2����ʹ��5��WENO, ��ֵ��С��3; ���ʹ��WENO7, ���ֵ��С��4��
   integer,parameter:: Scheme_UD1=0, Scheme_NND2=1, Scheme_UD3=2,Scheme_MUSCL2U=3,Scheme_MUSCL2C=4,   &
                       Scheme_MUSCL3=5,Scheme_OMUSCL2=6,Scheme_WENO5=7,Scheme_UD5=8, Scheme_WENO7=9
   integer,parameter:: Scheme_CD2=20, Scheme_none=-1         ! ��ʹ�ã��߽磩��ʽ
   integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_Van_Leer=5,Flux_Ausm=6
   integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2


!---------------�߽�����-------------------------
!  integer,parameter:: BC_Wall=-10, BC_Farfield=-20, BC_Periodic=-30,BC_Symmetry=-40,BC_Outlet=-22
   integer,parameter::  BC_Wall=2, BC_Symmetry=3, BC_Farfield=4,BC_Inflow=5, BC_Outflow=6 
   integer,parameter::  BC_Wall_Turbo=201             ! Ҷ�ֻ�е�л�ϻ�ı߽����� �������ٶ�Ϊ0�� �������ϵ����ת�� 
   integer,parameter::  BC_Periodic=501, BC_Extrapolate=401      ! (��չ) ��Griggen .inp�ļ��Ķ����������������ע��
!                     �����Ա߽������������ó�BC_Peridodic=501, ��������BC_PeriodicL=-2, BC_PeriodicR=-3

!   ���ڱ߽簴���ڱ߽紦���� �����α�����ר�Ŵ����� Ҷ�ֻ�ģʽ�£����ڱ߽�Ҳ�����⴦��   

!  -1 �ڱ߽� ���������߽磩��  -2  �����ڱ߽磻 -3 �����ڱ߽� ���������ұ��ڴ����������꣩
   integer,parameter::  BC_Inner=-1, BC_PeriodicL=-2, BC_PeriodicR=-3   
   
   integer,parameter::  BC_Zero=0       ! �ޱ߽�����

!                 �û��Զ���ı߽�������Ҫ����� >=900  
   integer,parameter::  BC_USER_FixedInlet=901, BC_USER_Inlet_time=902       !������������� �������ʱ������
   integer,parameter:: BC_USER_Blow_Suction_Wall=903    ! �����Ŷ�����

   integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
   integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2, Turbulence_SST=3,Turbulence_NewSA=21
   integer,parameter:: Init_continue=1, Init_By_FreeStream=0, Init_By_Zeroflow=-1,  Smooth_2nd=0,Smooth_4th=1 
!   real(PRE_EC), parameter::  Density_LIMIT=1.d-4,Temperature_LIMIT=1.d-4,Pressure_LIMIT=1.d-4

   integer,parameter:: Method_FVM=0, Method_FDM=1        ! ��֡��������
   integer,parameter:: FD_WENO5=1,FD_WENO7=2,FD_OMP6=3  ! ��ַ����õ���ֵ��ʽ
   integer,parameter:: FD_Steger_Warming=1,FD_Van_Leer=2


  end module const_var

! ������ (�߽����ӣ� ����飬 "����")
 module Mod_Type_Def
   use precision_EC
   implicit none

    TYPE BC_MSG_TYPE              ! �߽�������Ϣ
 !   integer::  f_no, face, ist, iend, jst, jend, kst, kend, neighb, subface, orient   ! BXCFD .in format
     integer:: ib,ie,jb,je,kb,ke,bc,face,f_no                      ! �߽��������棩�Ķ��壬 .inp format
     integer:: ib1,ie1,jb1,je1,kb1,ke1,nb1,face1,f_no1             ! ��������
	 integer:: L1,L2,L3                     ! ����ţ�����˳��������
   END TYPE BC_MSG_TYPE

!------------------------------------�����--------------------------------------
   TYPE Block_TYPE                                 ! ���ݽṹ������� ���������α�����������������Ϣ 
     integer::  Block_no,mpi_id           ! ��ţ������Ľ��̺�
	 integer::  nx,ny,nz                   ! ������nx,ny,nz
	 integer::  subface                   ! ������
 !   ������  
     real(PRE_EC),pointer,dimension(:,:,:):: x,y,z     ! coordinates of vortex, ����ڵ�����
     real(PRE_EC),pointer,dimension(:,:,:):: xc,yc,zc  ! coordinates of cell center, ������������ 
     real(PRE_EC),pointer,dimension(:,:,:):: Vol,Si,Sj,Sk ! Volume and surface area, ������������i,j,k���������߽������� 
	 real(PRE_EC),pointer,dimension(:,:,:):: ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3  !  i,j,k��������������ķ�����
!      Jocabian �任ϵ��	 
	 real(PRE_EC),pointer,dimension(:,:,:):: ix1,iy1,iz1,jx1,jy1,jz1,kx1,ky1,kz1
	 real(PRE_EC),pointer,dimension(:,:,:):: ix2,iy2,iz2,jx2,jy2,jz2,kx2,ky2,kz2
	 real(PRE_EC),pointer,dimension(:,:,:):: ix3,iy3,iz3,jx3,jy3,jz3,kx3,ky3,kz3
	 real(PRE_EC),pointer,dimension(:,:,:):: ix0,iy0,iz0,jx0,jy0,jz0,kx0,ky0,kz0
!     ������
	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U,Un,Un1    ! �غ���� (��ʱ�䲽��ǰһ������ʱ�䲽��ֵ), conversation variables 
     real(PRE_EC),pointer,dimension(:,:,:,:) :: Res         ! �в� ����ͨ����
     real(PRE_EC),pointer,dimension(:,:,:):: dt          ! (�ֲ�)ʱ�䲽��
     real(PRE_EC),pointer,dimension(:,:,:,:) :: QF      ! ǿ�Ⱥ��� (���������д�����ʹ��)
     real(PRE_EC),pointer,dimension(:,:,:,:) :: deltU   ! �غ�����Ĳ�ֵ, dU=U(n+1)-U(n)  ��������ʹ��
     real(PRE_EC),pointer,dimension(:,:,:,:) :: DU      ! U(n+1)-U(n)    LU-SGS��ʹ��
	 real(PRE_EC),pointer,dimension(:,:,:):: dw         ! turbulent viscous ; distance to the wall  (used in SA model)
	 real(PRE_EC),pointer,dimension(:,:,:):: surf1,surf2,surf3,surf4,surf5,surf6  ! �߽紦��ͨ��������������������ʹ�ã�
	 real(PRE_EC),pointer,dimension(:,:,:):: mu,mu_t    ! ����ճ��ϵ��������ճ��ϵ��
     real(PRE_EC),pointer,dimension(:,:,:):: dtime_mesh   ! ʱ�䲽������ �������������������

	 real(PRE_EC),pointer,dimension(:,:,:,:) :: U_average ! ʱ��ƽ���� ��d,u,v,w,T��
	 
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg     ! �߽�������Ϣ 
     integer,pointer,dimension(:,:,:):: BcI,BcJ,BcK     ! �߽�ָʾ�� �������߽� or �ڱ߽磩
   	 
	 integer:: IFLAG_FVM_FDM               ! ��� or �������
	 integer:: IF_OverLimit                          ! ���������ޣ��縺�¶ȣ��� ��Ҫ���;��ȣ�1�ף�
	End TYPE Block_TYPE  

!---------------------------���� -------------------------------------------------------- 
!  (�絥������ֻ��1�ף���������񣬿����ж���) 
  
   TYPE Mesh_TYPE                     ! ���ݽṹ�����񡱣� �������α���������������Ϣ
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep          ! ������ (1��Ϊ��ϸ����2��Ϊ������ 3��Ϊ��������...)�����������������Ŀ(��������), ʱ�䲽 
     integer:: NVAR       ! ��������Ŀ����������5��+ 0��1��2���������� ��BLģ��0����SAģ��1����SSTģ��2������ ������ʹ������ģ��
	 real(PRE_EC)::  tt                   !  �ƽ���ʱ��
	 real(PRE_EC),pointer,dimension(:)::  Res_max,Res_rms                 ! ���в�������в�, �ƽ���ʱ��
	 TYPE (Block_TYPE),pointer,dimension(:):: Block       ! ������顱  �������ڡ����񡱣�

!                                                       ���Ʋ��������ڿ�����ֵ������ͨ������������ģ�͵�    
!             ��Щ���Ʋ��������ڡ����񡱣���ͬ�����񡱿��Բ��ò�ͬ�ļ��㷽��������ģ�͵ȡ�	 �����磬�������õ;��ȷ�����������ʹ������ģ��,...��
!  If_dtime_mesh     �Ƿ��������������;ֲ�ʱ�䲽�����ֲ�ʱ�䲽������Ч��	
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
! ȫ�ֱ�����������������������������
!========================================================================================

  module Global_Var    
   use const_var        ! ����
   use mod_type_def     ! �߽�����
   implicit none


!---------------------------------------------------------------------------------------------
! global variables                                       ���ӳ�����ɼ���ȫ�ֱ���
!----------------------------------------------------------------------------

   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh                          ! ������ ������
   integer,save:: Num_Mesh,NVAR , Total_block, Num_block                      ! ��������� ���������� ���������, ��mpi���̵��������  
   integer,pointer,dimension(:):: bNi,bNj,bNk                            ! �����ά����ȫ�֣��� bNi(k)Ϊ��ȫ�֣���k���nx
   integer,save::  Kstep_save, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method, &
	  Kstep_show,If_viscous,If_Residual_smoothing,Mesh_File_Format,IF_Debug, &
	  Kstep_smooth,Kstep_init_smooth,NUM_THREADS,If_dtime_mesh, Step_Inner_Limit, &  
      Bound_Scheme, &                         ! �߽��ʽ
      IF_Walldist, IFLAG_LIMIT_FLOW, &             ! �Ƿ���Ҫ��ȡ ���������; �Ƿ���Ҫ����ѹ��������
      IF_Scheme_Positivity,          &    ! ����ֵ������ѹ�����ܶ��Ƿ�Ǹ�������ʹ��1��ӭ�磻
      Kstep_average,                 &          ! ʱ��ͳ�ƵĲ�������� 0 ��ͳ��
      Iflag_savefile                       ! 0 ���浽flow3d.dat, 1 ���浽flow3d-xxxxxxx.dat
   integer,save:: IF_TurboMachinary , Ref_medium_usrdef   !  �Ƿ�����Ҷ�ֻ�е�����ģʽ��0 or 1���� �Ƿ�ʹ���û��Զ������ ��Ĭ�Ͽ������ʣ�
   integer,save:: IF_InnerFlow   ! ����ģʽ (��ڱ߽�������Ҷ�ֻ�ģʽ����)

   integer,save:: FD_Flux,FD_scheme   ! ��Ƕ��ַ����õ�ͨ����ʽ����ֵ��ʽ
   integer,save:: KRK=0            ! Runge-Kutta�����е� �Ӳ�
   integer,save:: Istep_average=0  ! ͳ�Ʋ���
 ! global parameter (for all Meshes )                     ��������, ��ȫ�塰���񡱶�����

   real(PRE_EC),save:: Ma,Re,gamma,Cp,Cv,t_end,P_OUTLET,&
                       A_alfa,A_beta,PrL,PrT,T_inf,Twall,w_LU,Kt_inf,Wt_inf , Res_Inner_Limit, MUT_MAX,AoA,Aos
   real(PRE_EC),save:: Turbo_Periodic_seta, Turbo_w, Turbo_P0, Turbo_T0 , Turbo_L0   ! �����ȣ��ǣ��� ת�٣���ѹ������, �ο�����
   real(PRE_EC),save:: Periodic_dX,Periodic_dY,Periodic_dZ

 
 ! ȫ�ֿ��Ʋ�����������ֵ������ͨ������������ģ�͵� ����Щֻ����ϸ������Ч��
   real(PRE_EC),save :: Ralfa(3), Rbeta(3) , Rgamma(3),dt_global,CFL,dtmax,dtmin    ! RK�����еĳ�������ʱ�䲽���йص���
   integer,save:: Pre_Step_Mesh(3)                                                 ! ������ֵʱ��������Ԥ��������  
   integer,save:: Cood_Y_UP             ! 1 Y�ᴹֱ���ϣ� 0 Z�ᴹֱ����
   integer,save:: Pdebug(4)             ! debug ʹ�ã� ���ĳ���ĳһ���ֵ
   real(PRE_EC),save:: Ref_S,Ref_L, Centroid(3)                           ! ������������ʹ�ã��ο����, �ο�����, ��������
   real(PRE_EC),save:: Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax                ! ���ܶȡ�ѹ�����ٶȡ���SAģ���б������������� (��С�����ֵ)
   real(PRE_EC),save:: CP1_NSA,CP2_NSA       ! parameters in New SA model
 !-----------mpi data ----------------------------------------------------------- 
   integer:: my_id,Total_proc                   ! my_id (�����̺�), Total_proc �ܽ�����
   integer,pointer,dimension(:):: B_Proc, B_n    ! B_proc(m) m�����ڵĽ��̺�; B_n(m) m�����ڽ����е��ڲ����
   integer,pointer,dimension(:):: my_Blocks       ! �����̰����Ŀ��
  end module Global_Var  
!----------------------------------------------------------------------------



! ���-���������Ϸ���ʹ��
!------------SEC part---------------------------------------------------------------------

   module FDM_data
    use precision_EC
    implicit none
    TYPE FDM_Block_TYPE                              ! ���ݽṹ������� ���������α�����������������Ϣ 
    real(PRE_EC), pointer,dimension(:,:,:):: ix,iy,iz,jx,jy,jz,kx,ky,kz,Jac   ! Jocabian�任ϵ��
    End TYPE FDM_Block_TYPE  
  
   TYPE FDM_Mesh_TYPE                     ! ���ݽṹ�����񡱣� �������α���������������Ϣ
	 TYPE (FDM_Block_TYPE),pointer,dimension(:):: Block       ! ������顱  �������ڡ����񡱣�
   End TYPE FDM_Mesh_TYPE 
  
   TYPE (FDM_Mesh_TYPE),pointer,dimension(:):: FDM_Mesh       ! ������ ������
  
   end  module FDM_data
