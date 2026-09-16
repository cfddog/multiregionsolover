$control_ec
!=======================================================================
!  2D circular cylinder, Re = rho*U*D/mu = 1*1*1/0.025 = 40
!  UPPER HALF only (theta = 0..pi) + symmetry plane y = 0 (BC 3),
!  body-fitted polar O-grid, single block, artificial-compressibility
!  (AC) low-speed solver:  LS_Algorithm = 3   (BLOCK_LOWSPEED)
!   -> non-dimensionalisation: D = 1, rho = 1, U = 1, mu = 1/Re
!  BCs (bc3d.inp): wall(2) on i=1 | outlet(6)/inlet(5) on r=R_out |
!                  symmetry(3) on j=1,j=ny (y=0) and k=1,nz (2D)
!=======================================================================
  Iflag_init=0
  Ma=0.1d0
  Re=40.d0
  AoA=0.d0
  If_viscous=0            ! inviscid (Euler) - BC_Wall then acts as a slip wall
  Iflag_turbulence_model=0
  Num_Mesh=1
  Mesh_File_Format=2
  Kstep_save=1
  Kstep_show=1
  Kstep_average=0
  Kstep_smooth=-1
  Kstep_init_smooth=0
  t_end=1.d0
  dt_global=1.d0
  CFL=1.d0
  Time_Method=0
  T_inf=288.15d0
  Twall=-1.d0
  Iflag_Scheme=5
  Iflag_Flux=5
  IFlag_Reconstruction=0
  Ref_S=1.d0
  Ref_L=1.d0
  Lscale=1.d0
  Cood_Y_UP=1
  NUM_THREADS=1
  IF_Debug=0
  IFLAG_LIMIT_FLOW=0
  Pdebug=1,1,1,1
  Ldmin=1.d-6
  Ldmax=1000.d0
  Lpmin=1.d-6
  Lpmax=1000.d0
  Lumax=1000.d0
  LSAmax=1000.d0
  CP1_NSA=0.2d0
  CP2_NSA=100.d0
  Periodic_dX=0.d0
  Periodic_dY=0.d0
  Periodic_dZ=0.d0
  IF_TurboMachinary=0
  Ref_medium_usrdef=0
  IF_Scheme_Positivity=1
  Turbo_P0=101330.d0
  Turbo_T0=288.15d0
  Turbo_L0=1.d0
  Turbo_w=0.d0
  Turbo_Periodic_seta=0.d0
  IF_Innerflow=0
  Iflag_savefile=0
  Iflag_vtk_onefile=0
  Iflag_vtk_SI=0
  Iflag_bc_check=0
  gamma=1.4d0
  PrL=0.7d0
  PrT=0.9d0
!---------------- low-speed (incompressible) parameters ----------------
  LS_rho=1.0d0
  LS_mu=0.025d0
  LS_k=0.0d0
  LS_Cp=1.0d0
  LS_T_ref=300.d0
  LS_Inlet_Type=1
  LS_U_in=1.0d0
  LS_V_in=0.d0
  LS_W_in=0.d0
  LS_Mdot_in=0.d0
  LS_P_in=0.d0
  LS_P_out=0.d0
  LS_T_wall=-1.d0
  LS_U_lid=0.d0
  LS_alpha_p=0.3d0
  LS_alpha_u=0.7d0
  LS_alpha_T=0.7d0
  LS_Max_Iter=5000
  LS_Tol=1.d-8
  LS_Scheme=3
  LS_Algorithm=3
!---------------- AC (artificial compressibility) solver ---------------
  AC_beta=1.d0
  AC_CFL=20.d0
  AC_CFLv=0.5d0
  AC_Max_Iter=300000
  AC_Print=2000
  AC_Tol=1.d-8
  AC_w=1.0d0
  AC_Flux=3
  AC_Recon=1
  AC_Limiter=1
  AC_WenoBlend=0.d0
$end
