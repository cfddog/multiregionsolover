!  -----------------------????????????????????----------------------
  subroutine read_parameter
   use Global_var
   implicit none
   logical ext1
!-----------------------------------------------------------------------------
   call set_default_parameter ! ???��???????

   if(my_id .eq. 0) then
     inquire(file="control.ec",exist=ext1)
      if(ext1) then           ! ??????control.ec (Namelist ??????????) 
       call read_parameter_ec
      else
       print*, "Can not find 'control.ec', stop !"
	   stop
	  endif

 
   endif

   call bcast_para      ! ???????????
   
   call set_const_para
!--------------------------------------------------------------------

  end subroutine read_parameter
   
!--------------------------------------------------------------------    
! ???��?????????  
  subroutine set_default_parameter 
   use Global_var
   implicit none
    Ma=1.d0        ! Mach number     
    Re= 1000.d0    ! Reynolds number
    gamma=1.4d0    ! 

	AOA=0.d0       ! Angle of attack    
	AOS=0.d0       ! Angle of Slide
	P_outlet=-1.d0   ! Outlet pressure (<0 extrapolation)  
	t_end=100.d0     ! End time (non-dimensional)
	Kstep_save=1000  ! Save data per xxx steps
    Iflag_turbulence_model=0   ! turbulence model (0 none, 1 BL, 2 SA, 3 SST)
	Iflag_init=0  ! 0 ???????????????????????? 1  ???? -1 ??0 ???????????
    If_viscous=1  ! 0 ????? 1 ???
    Iflag_local_dt=1   ! 0 ????????  1 ????????
    dt_global=0.01     ! Global time step
    CFL=1.d0           ! CFL number 
    dtmax=10.d0       !Limit of maximum time step 
	dtmin=1.d-9       !Limit of minimum time step
	Time_Method=0     ! Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0, Time_dual_LU_SGS=-1
    If_Residual_smoothing=0   ! 0 Do not need smoothing
    w_LU=1.d0         ! factor in LU-SGS
    If_dtime_mesh=1   ! decrease time step when grid quality is not good
    Iflag_Scheme= 5    ! 0 UD1;  1 NND2 ; 2 UD3 ; 3 MUSCL2U ;  4 MUSCL2C; 5 MUSCL3; 6 OMUSCL2; 7 WENO5 ; 8 UD5; 9 WENO7 
    Iflag_Flux=5      ! 1 Steger_Warming; 2 HLL, 3 HLLC,  4 Roe, 5  Van_Leer, 6 Ausm
    IFlag_Reconstruction=0   !  0 Original, 1 Conservative, 2 Characteristic
    Mesh_File_Format=0    ! 0 unformatted, 1 formatted
    Kstep_show=1          ! show per xxx steps
    Kstep_average=0       ! do not average
    Kstep_smooth=-1       ! Smoothing flow per xxx steps (<0 donot smooth)
    Kstep_init_smooth=0   ! Smoothing flow at initial time
    Num_Mesh=1            ! Number of mesh (1 single-grid), 2,3 multi-grid
    T_inf=288.15d0        ! Reference temperature  (in viscous coefficient)
    Twall=-1.d0           ! Wall temperature (in K degree) (<0 adabitic)
    Kt_inf=1.d-5          ! initial of Kt for SST model
    Wt_inf=0.01d0         ! Initial of Wt for SST model
    IF_Debug=0            ! 1 Debug mode
    NUM_THREADS=1         ! Threads for OpenMP
    Step_Inner_Limit=20   ! Inner time advance limit (Dual-time)
    Res_Inner_Limit=1.d-10 ! Inner Resdial limit (Dual-time)
    MUT_MAX=-1.d0          ! Limit for vt/vs (<0 no limit)
	Bound_Scheme= Scheme_MUSCL2C         ! Boundary scheme (Default: MUSCL2C)
    Pre_Step_Mesh(1:3)=0   ! Pre step for multi-grid
    Ref_S=1.d0             ! Ref. area
	Ref_L=1.d0             ! Ref. length
	Lscale=1.d0            ! Length scale: 1=m, 0.001=mm, etc.
	Centroid(1:3)=0.d0  ! Centroid coordinate  
    Cood_Y_UP=1   !       ???Y???????
	  IFLAG_LIMIT_FLOW=1         ! ??????????????????????????څ?1 ?? 2016-10-21
	Pdebug(1:4)=1
    PrL=0.7d0          ! Linear Prandtl number
	PrT=0.9d0          ! Turbulent Prandtl number
    Ldmin=1.d-6
	Ldmax=1000.d0
	Lpmin=1.d-6
	Lpmax=1000.d0
	Lumax=1000.d0
	LSAmax=1000.d0
	CP1_NSA=0.2d0   ! for New SA
	CP2_NSA=100.d0
	Periodic_dX=0.d0   ! ??????????????
	Periodic_dY=0.d0 
	Periodic_dZ=0.d0

    Iflag_savefile=0       ! ???��??flow3d.dat
    Iflag_vtk_onefile=0    ! 0=one vtk per block (flow3d_block_*.vtk), 1=single merged flow3d.vtk
    Iflag_vtk_SI=0         ! 0=compressible blocks written non-dimensional, 1=convert to SI units in vtk
!----for Turbomachinary solver------------
    IF_TurboMachinary=0    ! ??????????????
	Ref_medium_usrdef=0    ! ???????????? ???????????
    IF_Scheme_Positivity=1     ! ?????????????????????????????????1????��
    Turbo_P0= 101330.d0    ! ??? ??????1?????????
	Turbo_T0= 288.15d0     ! ???? ?????288.15K)
    Turbo_L0= 1.d0         ! ?��????? ??????1m)
    Turbo_w=0.d0           ! ???  ( ?/?? ?? ???0)
    Turbo_Periodic_seta=0.d0   ! ?????????????
!-------------------------------------
    IF_InnerFlow=0   ! ??????

!---- Low-speed (incompressible) solver parameters (SI units) -------------
    LS_rho=1.2d0        ! density [kg/m3]
    LS_mu=1.8d-5        ! dynamic viscosity [Pa.s]
    LS_k=0.025d0        ! thermal conductivity [W/(m.K)]
    LS_Cp=1005.d0       ! specific heat [J/(kg.K)]
    LS_T_ref=288.15d0   ! reference/inlet temperature [K]
    LS_Inlet_Type=1     ! 1=velocity inlet, 2=mass-flow inlet, 3=pressure inlet
    LS_U_in=1.d0; LS_V_in=0.d0; LS_W_in=0.d0   ! inlet velocity [m/s]
    LS_Mdot_in=0.d0     ! inlet mass flow rate [kg/s] (type 2)
    LS_P_in=101325.d0   ! inlet total pressure [Pa] (type 3)
    LS_P_out=101325.d0  ! outlet static pressure [Pa]
    LS_T_wall=-1.d0     ! wall temperature [K] (<0 adiabatic)
    LS_U_lid=0.d0       ! lid-driven cavity: tangential wall velocity on j+ face [m/s]
    LS_alpha_p=0.3d0; LS_alpha_u=0.7d0; LS_alpha_T=0.7d0   ! under-relaxation
    LS_Max_Iter=5000    ! SIMPLE inner iterations per solver call
    LS_Tol=1.d-8        ! SIMPLE convergence tolerance (pressure correction residual)
    LS_Scheme=1         ! convection scheme: 1=1st-order upwind, 2=2nd-order upwind, 3=MUSCL(Van Leer)
   LS_Algorithm=1      ! pressure-velocity coupling: 1=SIMPLE, 2=SIMPLEC, 3=AC-FV

!---- AC (artificial compressibility) solver defaults (LS_Algorithm=3) --------
    AC_Max_Iter=40000   ! pseudo-time iterations per solver call
    AC_Print=500        ! print interval
    AC_beta=10.d0       ! beta = AC_beta*U_ref^2
    AC_CFL=2.d0         ! inviscid CFL
    AC_CFLv=0.5d0       ! viscous CFL
    AC_Tol=1.d-7        ! convergence tolerance (dimensionless residual)
    AC_w=1.d0           ! LU-SGS relaxation

!---- Porous-media solver parameters (SI units) --------------------------------
    Porous_T_ref=288.15d0   ! initial/reference solid-frame temperature [K]
    Porous_alpha_Ts=0.7d0   ! under-relaxation for the Ts (solid-frame) solve
    Porous_Max_Iter=5000    ! SIMPLE inner iterations per solver call
    Porous_Tol=1.d-8        ! SIMPLE convergence tolerance (pressure correction residual)

!---- interface-19 staggered (segmented) coupling controls --------------------
    Iflag_Couple_Scheme=0   ! 0 = per-step coupling; 1 = staggered segmented (gas chunk / porous chunk)
    Kstep_Couple_Comp=1000  ! compressible steps per gas chunk
    Niter_Couple_Outer=60   ! max outer staggered iterations
    Porous_Chunk_Iter=10    ! porous SIMPLE solver calls per porous chunk
    Niter_Couple_Warm=0     ! outer iters with the FULL gas chunk (warm start, auto-shrink afterwards)
    Kstep_Couple_Min=1      ! gas-chunk step floor after warm-up (auto shrink)
    Twall_Couple_Init=300.d0 ! first gas-chunk interface wall temperature [K]
    Tol_Couple_Tw=2.d-2     ! outer convergence: max|dT_w| [K]
    Tol_Couple_p=1.d1       ! 12 (fluid-fluid): outer convergence max|dp_w| [Pa]
    Tol_Couple_u=1.d-2      ! 12 (fluid-fluid): outer convergence max|du_w| [m/s]
    Iflag_Couple_WallFlux=0 ! 11/13 CHT split: 0=couple-based, 1=isothermal-Tw/qw wall-flux
end

!------read parameter (Namelist type)---------------- 
  subroutine read_parameter_ec 
   use Global_var
   implicit none
   real(PRE_EC):: R0, a0, d0, mu0,mu1

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
        Pre_Step_Mesh,Ref_S,Ref_L,Lscale,Centroid,Cood_Y_UP,IFLAG_LIMIT_FLOW,Pdebug, &
		Ldmin,Ldmax,Lpmin,Lpmax,Lumax,LSAmax,CP1_NSA,CP2_NSA, &
        IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, &
		Turbo_P0,Turbo_T0, Turbo_L0,Turbo_w, Turbo_Periodic_seta, &
		Periodic_dX, Periodic_dY, Periodic_dZ, &
		IF_Innerflow, Iflag_savefile, Iflag_vtk_onefile, Iflag_vtk_SI, &
		LS_rho, LS_mu, LS_k, LS_Cp, LS_T_ref, LS_Inlet_Type, &
		LS_U_in, LS_V_in, LS_W_in, LS_Mdot_in, LS_P_in, LS_P_out, &
		LS_T_wall, LS_U_lid, LS_alpha_p, LS_alpha_u, LS_alpha_T, &
		LS_Max_Iter, LS_Tol, LS_Scheme, LS_Algorithm, &
		AC_Max_Iter, AC_Print, AC_beta, AC_CFL, AC_CFLv, AC_Tol, AC_w, &
		Porous_T_ref, Porous_alpha_Ts, Porous_Max_Iter, Porous_Tol, &
		Iflag_Couple_Scheme, Kstep_Couple_Comp, Niter_Couple_Outer, &
		Porous_Chunk_Iter, Niter_Couple_Warm, Kstep_Couple_Min, &
		Twall_Couple_Init, Tol_Couple_Tw, Tol_Couple_p, Tol_Couple_u, Iflag_Couple_WallFlux


	open(99,file="control.ec")
	read(99,nml=control_ec)
    close(99)
 
 !---- convert parameters ----------------------
 ! Ref_medium_usrdef==0 ????????? (Ma=1, ???????????????? Re) ?? ==1 ??????????? ?????????Ma, Re???
 ! ??????????????????????

    if( (IF_TurboMachinary ==1 .or. IF_Innerflow==1 )    &  
	  .and.  Ref_medium_usrdef == 0) then   ! ???????????????Mach???? Reynolds??
      
	  T_inf=Turbo_T0  ! ?��???? ???????????
      gamma=1.4d0    ! 
	  PrL=0.7d0   ! Prandtl??
	  PrT=0.9d0
      R0= 287.06d0   ! ?????????�i??R
	  a0= sqrt(gamma*R0*Turbo_T0)    ! ?��??????????? 
	  mu0=1.179d-5     ! ?????????? (288.15K)  
      mu1=mu0* sqrt((Turbo_T0/288.15d0)**3)*(288.15d0+110.4d0)/(Turbo_T0+110.4d0)  ! ?��????????????????
      d0=Turbo_P0/(R0*Turbo_T0)
	  Re=d0*a0*Turbo_L0/mu1    ! ?��??????????????????Reynolds??
	  Ma=1.d0     ! Mach??    ????????????��?????????��?Mach???1??
      Turbo_w= 2.d0*PI*Turbo_w/(a0/Turbo_L0)   ! ?????????? Turbo_W???/??
!	  P_outlet=P_outlet/Turbo_P0    ! ???,  Bug !!
	  P_outlet=P_outlet/(d0*a0*a0)    ! ??? ?????? ???????
	endif
      Turbo_Periodic_seta=Turbo_Periodic_seta*PI/180.d0                 ! Turbo_Periodic_seta ???


 !---output paramters----------------------------
    open(99,file="output_para.out")
	write(99,*) "-------OpenCFD-EC (Ver 1.1t), (c) Li Xinliang, lixl@imech.ac.cn--"
	write(99,*) "-------------------------------------------"
	write(99,*) "Ma=", Ma,  "  Re= ", Re , "gamma=", gamma
	write(99,*) "PrL=", PrL, "PrT=",PrT

	write(99,*) "A_alfa=", A_alfa, " A_beta=", A_beta 
	write(99,*) "P_outlet=", P_outlet
	write(99,*) "t_end=", t_end, " Iflag_local_dt=", Iflag_local_dt
	write(99,*) "dt_global=",dt_global 
	write(99,*) "Time_Method=", Time_Method,  " CFL= ", CFL 
	write(99,*) "dtmax=", dtmax, " dtmin=", dtmin 
    write(99,*) "Iflag_init=", Iflag_init," If_viscous=", If_viscous
	write(99,*) "Iflag_turbulence_model=",Iflag_turbulence_model
    write(99,*) "Iflag_Scheme= ",Iflag_Scheme, " Iflag_Flux=", Iflag_Flux
	write(99,*) "IFlag_Reconstruction=", IFlag_Reconstruction, "Bound_Scheme= ", Bound_Scheme
    write(99,*) "Kstep_save= ",Kstep_save,  " Kstep_show=", Kstep_show, "Kstep_average=",Kstep_average
	write(99,*) "Kstep_smooth= ", Kstep_smooth, " Kstep_init_smooth=",Kstep_init_smooth
    write(99,*) "If_Residual_smoothing=", If_Residual_smoothing, " If_dtime_mesh=",If_dtime_mesh
	write(99,*) "w_LU=",w_LU
    write(99,*) "Mesh_File_Format=",Mesh_File_Format, " Num_Mesh=",Num_Mesh
    write(99,*) "T_inf=", T_inf, " Twall=", Twall
	write(99,*) "Kt_inf=",Kt_inf,  " Wt_inf=",Wt_inf
    write(99,*) "Step_Inner_Limit=",Step_Inner_Limit, " MUT_MAX=", MUT_MAX
	write(99,*) "Pre_Step_Mesh(:)=", Pre_Step_Mesh(1:Num_Mesh)
    write(99,*) "IF_Debug=",IF_Debug
	write(99,*) "Ref_S=", Ref_S, "Ref_L=", Ref_L
	write(99,*) "Lscale=", Lscale
	write(99,*) "Centroid=", Centroid(1:3)
	write(99,*) "Cood_Y_UP=",Cood_Y_UP
    write(99,*) "Periodic_dX, dY, dZ=", Periodic_dX,Periodic_dY,Periodic_dZ
	write(99,*) "IFLAG_LIMIT_FLOW=",IFLAG_LIMIT_FLOW
	write(99,*) "Ldmin,Ldmax,Lpmin,Lpmax,Lumax=", Ldmin,Ldmax,Lpmin,Lpmax,Lumax
	write(99,*) "CP1_NSA,CP2_NSA=",CP1_NSA,CP2_NSA
    write(99,*) "IF_TurboMachinary=",  IF_TurboMachinary
	write(99,*) "Turbo_Periodic_seta=", Turbo_Periodic_seta
    write(99,*) "Ref_medium_usrdef=", Ref_medium_usrdef
	write(99,*) "Turbo_w= ", Turbo_w, "Turbo_L0=", Turbo_L0
	write(99,*) "a0 (Reference velocity)=", a0
	write(99,*) "d0, T0, p0 (Reference values)=", d0,Turbo_T0, Turbo_P0
	write(99,*) "IF_Scheme_Positivity =", IF_Scheme_Positivity
	write(99,*) "Iflag_savefile=", Iflag_savefile
	write(99,*) "NUM_THREADS=",NUM_THREADS
    write(99,*) "--------------------------------------------"

    close(99)

 !----------------------------------------------
 end
!--------------------------------------------------
 

!------------------------------------------------------------------
  
  
   subroutine bcast_para
   use Global_var
   implicit none
    integer:: Ipara(100),ierr
    real(PRE_EC):: rpara(100)
    Ipara=0
	rpara=0.d0
!----
	rpara(1)=Ma 
	rpara(2)=Re
	rpara(3)=AoA
	rpara(4)=AoS
	rpara(5)=p_outlet
	rpara(6)=t_end
	rpara(7)=dt_global
	rpara(8)=CFL
	rpara(9)=dtmax
	rpara(10)=dtmin
	rpara(11)=w_LU
	rpara(12)=T_inf
	rpara(13)=Twall
	rpara(14)=Kt_inf
	rpara(15)=Wt_inf
    rpara(16)=Res_Inner_Limit
    rpara(17)=MUT_MAX
    rpara(18)=Ref_S
	rpara(19)=Ref_L
	rpara(20:22)=Centroid(1:3)
    rpara(39)=Lscale
	rpara(23)=Ldmin
	rpara(24)=Ldmax
	rpara(25)=Lpmin
	rpara(26)=Lpmax
	rpara(27)=Lumax
	rpara(28)=LSAmax
	rpara(29)=CP1_NSA
	rpara(30)=CP2_NSA
	rpara(31)=gamma
	rpara(32)=PrL
	rpara(33)=PrT
    rpara(34)=Turbo_w
    rpara(35)=Turbo_Periodic_seta
    rpara(36)=Periodic_dX
	rpara(37)=Periodic_dY
	rpara(38)=Periodic_dZ

    rpara(40)=LS_rho
    rpara(41)=LS_mu
    rpara(42)=LS_k
    rpara(43)=LS_Cp
    rpara(44)=LS_T_ref
    rpara(45)=LS_U_in
    rpara(46)=LS_V_in
    rpara(47)=LS_W_in
    rpara(48)=LS_Mdot_in
    rpara(49)=LS_P_in
    rpara(50)=LS_P_out
    rpara(51)=LS_T_wall
    rpara(52)=LS_U_lid
    rpara(53)=LS_alpha_p
    rpara(54)=LS_alpha_u
    rpara(55)=LS_alpha_T
    rpara(56)=LS_Tol
    rpara(57)=Porous_T_ref
    rpara(58)=Porous_alpha_Ts
    rpara(59)=Porous_Tol
    rpara(60)=Twall_Couple_Init
    rpara(61)=Tol_Couple_Tw
    rpara(67)=Tol_Couple_p
    rpara(68)=Tol_Couple_u
    rpara(62)=AC_beta
    rpara(63)=AC_CFL
    rpara(64)=AC_CFLv
    rpara(65)=AC_Tol
    rpara(66)=AC_w




    Ipara(1)=Kstep_save
	Ipara(2)=Iflag_turbulence_model
	Ipara(3)=Iflag_init
	Ipara(4)=If_viscous
    Ipara(5)=Iflag_local_dt
	Ipara(6)=Time_Method
	Ipara(7)=If_Residual_smoothing
	Ipara(8)=If_dtime_mesh
    Ipara(9)=Iflag_Scheme
	Ipara(10)=Iflag_Flux
	Ipara(11)=IFlag_Reconstruction
	Ipara(12)=Mesh_File_Format
	Ipara(13)=Kstep_show
	Ipara(14)=Kstep_smooth
	Ipara(15)=Kstep_init_smooth
    Ipara(16)=Num_Mesh
	Ipara(17)=IF_Debug
	Ipara(18)=NUM_THREADS
    Ipara(19)=Step_Inner_Limit
    Ipara(20)=Bound_Scheme
    Ipara(21)=Cood_Y_UP
	Ipara(22)=IFLAG_LIMIT_Flow
    Ipara(23:26)=Pdebug(1:4)
    Ipara(27)=IF_TurboMachinary
	Ipara(28)=IF_Scheme_Positivity
   	Ipara(29)=IF_Innerflow
   	Ipara(30)=Kstep_average
    Ipara(31)=Iflag_savefile
    Ipara(32)=LS_Inlet_Type
    Ipara(33)=LS_Max_Iter
    Ipara(34)=LS_Scheme
    Ipara(35)=LS_Algorithm
    Ipara(36)=Porous_Max_Iter
    Ipara(37)=Iflag_vtk_onefile
    Ipara(38)=Iflag_vtk_SI
    Ipara(39)=Iflag_Couple_Scheme
    Ipara(40)=Kstep_Couple_Comp
    Ipara(41)=Niter_Couple_Outer
    Ipara(42)=Porous_Chunk_Iter
    Ipara(43)=Niter_Couple_Warm
    Ipara(44)=Kstep_Couple_Min
    Ipara(45)=AC_Max_Iter
    Ipara(46)=AC_Print
    Ipara(47)=Iflag_Couple_WallFlux

	 call MPI_bcast(rpara,100,OCFD_DATA_TYPE,0,  MPI_COMM_WORLD,ierr)
	 call MPI_bcast(Ipara,100,MPI_Integer,0,  MPI_COMM_WORLD,ierr)

	Ma=rpara(1) 
	Re=rpara(2)
	AoA=rpara(3)
	AoS=rpara(4)
	p_outlet=rpara(5)
	t_end=rpara(6)
	dt_global=rpara(7)
	CFL=rpara(8)
	dtmax=rpara(9)
	dtmin=rpara(10)
	w_LU=rpara(11)
	T_inf=rpara(12)
	Twall=rpara(13)
	Kt_inf=rpara(14)
	Wt_inf=rpara(15)
    Res_Inner_Limit=rpara(16)
    MUT_MAX=rpara(17)
    Ref_S=rpara(18)
	Ref_L=rpara(19)
	Centroid(1:3)=rpara(20:22)
    Lscale=rpara(39)
	Ldmin=rpara(23)
	Ldmax=rpara(24)
	Lpmin=rpara(25)
	Lpmax=rpara(26)
	Lumax=rpara(27)
	LSAmax=rpara(28)
	CP1_NSA=rpara(29)
	CP2_NSA=rpara(30)
	gamma=rpara(31)
	PrL=rpara(32)
	PrT=rpara(33)
    Turbo_w=rpara(34)
    Turbo_Periodic_seta=rpara(35)
    Periodic_dX=rpara(36)
	Periodic_dY=rpara(37)
	Periodic_dZ=rpara(38)
    LS_rho=rpara(40)
    LS_mu=rpara(41)
    LS_k=rpara(42)
    LS_Cp=rpara(43)
    LS_T_ref=rpara(44)
    LS_U_in=rpara(45)
    LS_V_in=rpara(46)
    LS_W_in=rpara(47)
    LS_Mdot_in=rpara(48)
    LS_P_in=rpara(49)
    LS_P_out=rpara(50)
    LS_T_wall=rpara(51)
    LS_U_lid=rpara(52)
    LS_alpha_p=rpara(53)
    LS_alpha_u=rpara(54)
    LS_alpha_T=rpara(55)
    LS_Tol=rpara(56)
    Porous_T_ref=rpara(57)
    Porous_alpha_Ts=rpara(58)
    Porous_Tol=rpara(59)
    Twall_Couple_Init=rpara(60)
    Tol_Couple_Tw=rpara(61)
    AC_beta=rpara(62)
    AC_CFL=rpara(63)
    AC_CFLv=rpara(64)
    AC_Tol=rpara(65)
    AC_w=rpara(66)
    Tol_Couple_p=rpara(67)
    Tol_Couple_u=rpara(68)



    Kstep_save=Ipara(1)
	Iflag_turbulence_model=Ipara(2)
	Iflag_init=Ipara(3)
	If_viscous=Ipara(4)
    Iflag_local_dt=Ipara(5)
	Time_Method=Ipara(6)
	If_Residual_smoothing=Ipara(7)
	If_dtime_mesh=Ipara(8)
    Iflag_Scheme=Ipara(9)
	Iflag_Flux=Ipara(10)
	IFlag_Reconstruction=Ipara(11)
	Mesh_File_Format=Ipara(12)
	Kstep_show=Ipara(13)
	Kstep_smooth=Ipara(14)
	Kstep_init_smooth=Ipara(15)
    Num_Mesh=Ipara(16)
	IF_Debug=Ipara(17)
	NUM_THREADS=Ipara(18)
    Step_Inner_Limit=Ipara(19)
    Bound_Scheme=Ipara(20)
    Cood_Y_UP=Ipara(21)
	IFLAG_LIMIT_FLOW=Ipara(22)
    Pdebug(1:4)=Ipara(23:26)
    IF_TurboMachinary=Ipara(27)
	IF_Scheme_Positivity=Ipara(28)
   	IF_Innerflow=Ipara(29)
   	Kstep_average=Ipara(30)
    Iflag_savefile=Ipara(31)
    LS_Inlet_Type=Ipara(32)
    LS_Max_Iter=Ipara(33)
    LS_Scheme=Ipara(34)
    LS_Algorithm=Ipara(35)
    Porous_Max_Iter=Ipara(36)
    Iflag_vtk_onefile=Ipara(37)
    Iflag_vtk_SI=Ipara(38)
    Iflag_Couple_Scheme=Ipara(39)
    Kstep_Couple_Comp=Ipara(40)
    Niter_Couple_Outer=Ipara(41)
    Porous_Chunk_Iter=Ipara(42)
    Niter_Couple_Warm=Ipara(43)
    Kstep_Couple_Min=Ipara(44)
    AC_Max_Iter=Ipara(45)
    AC_Print=Ipara(46)
    Iflag_Couple_WallFlux=Ipara(47)


    call MPI_bcast(Pre_Step_Mesh,Num_Mesh,MPI_Integer,0,  MPI_COMM_WORLD,ierr)
    end  subroutine bcast_para
  
  
      
! ?څ????
   subroutine set_const_para 
   use Global_var
   implicit none

   Twall=Twall/T_inf                         ! wall temperature
   Lpmin=Lpmin/(gamma*Ma*Ma)                 ! rato of free-stream pressure
   Lpmax=Lpmax/(gamma*Ma*Ma)



   if(Bound_Scheme== Scheme_none)  Bound_Scheme=Iflag_Scheme    ! ?��?????????????????????   
   if(If_viscous .eq. 0) Iflag_turbulence_model=Turbulence_NONE !    ????????????????????????
    
   if(Iflag_turbulence_model .eq. Turbulence_SA .or. & 
      Iflag_turbulence_model .eq. Turbulence_NewSA) then
     NVAR=6                       ! SA???????6??????
   else if (Iflag_turbulence_model .eq. Turbulence_SST) then
     NVAR=7                       ! SST ???????7??????
   else
     NVAR=5                       ! 5??????
   endif

  if(Iflag_turbulence_model .eq. 0) then
   IF_Walldist=0                   ! ?????????
  else
   IF_Walldist=1
  endif
 
!--------------------------------------------------------------------------------
  AoA=AoA*PI/180.d0   ! Angle of attack
  AoS=AoS*PI/180.d0   ! Angle of Slide
 if(Cood_Y_UP ==1) then   ! Y ??????? or Z??????? 
   A_alfa=AoA             ! Y??????? A_alfa?????
   A_beta=AoS
 else
   A_alfa=AoS  
   A_beta=AoA
 endif
   
  
   Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma) 
   Cp=Cv*gamma 
!--------------------------------------------------------------------

!  ????????????
!   if(Iflag_turbulence_model .eq. Turbulence_SA .or. Iflag_turbulence_model .eq. Turbulence_SST) then
!     inquire(file="wall_dist.dat",exist=file_exist)
!	 if( .not. file_exist) then
!	   call  comput_wall_dist  
!     endif
!   endif
    
	if(Num_Mesh .ne. 1) then
	   if(Time_Method .eq. Time_LU_SGS .or. Time_Method .eq. Time_Dual_LU_SGS) then
	    print*, "In this version, LU_SGS (or Dual_time_LU_SGS) DO NOT Support Multi-Grid !!!"
		stop
	   endif
	endif
  end

  subroutine read_material_in
   use Global_var
   implicit none
   integer:: m,bt
   logical:: ex
   if(my_id .ne. 0) return
   inquire(file="material.in",exist=ex)
   if(.not. ex) then
     print*, "material.in not found, all blocks set to FLUID (0)"
     return
   endif
   open(99,file="material.in")
   read(99,*) Total_block
   allocate(Block_Type_List(Total_block))
   read(99,*) (Block_Type_List(m), m=1,Total_block)
   allocate(solid_rho_list(Total_block), solid_Cp_list(Total_block), solid_k_list(Total_block))
   do m=1,Total_block
     read(99,*,end=100) solid_rho_list(m), solid_Cp_list(m), solid_k_list(m)
     cycle
100  solid_rho_list(m)=0.d0; solid_Cp_list(m)=0.d0; solid_k_list(m)=0.d0
   enddo
   close(99)
   print*, "Read material.in OK, Total_block=", Total_block
   do m=1,Total_block
     print*, "  Block", m, ": Block_type=", Block_Type_List(m), &
             " rho=", solid_rho_list(m), " Cp=", solid_Cp_list(m), " k=", solid_k_list(m)
   enddo
    call read_porous_prop
  end subroutine read_material_in
!----------------------------------------------------------------------
! Read porous material properties from porous.inp (block-uniform first version)
! Format:
!   Total_Porous_Blocks
!   Block_no  eps  dp(m)  hv(W/m3/K)
!   ...
! eps in (0,1], dp>0. hv<=0 disables the inter-phase heat exchange (frozen frame).
! Blocks not listed keep defaults eps=0.9, dp=1e-3 m, hv=0.
! Called on root inside read_material_in; broadcast follows in bcast_material.
!----------------------------------------------------------------------
  subroutine read_porous_prop
    use Global_var
    use const_var
    implicit none
    integer:: m, n, i
    logical:: ex
    allocate(porous_eps_list(Total_block), porous_dp_list(Total_block), &
             porous_hv_list(Total_block))
    porous_eps_list(:)=0.9d0
    porous_dp_list(:)=1.d-3
    porous_hv_list(:)=0.d0
    inquire(file="porous.inp", exist=ex)
    if(.not. ex) then
      print*, "porous.inp not found: porous blocks use defaults (eps=0.9, dp=1e-3 m, hv=0)"
      return
    endif
    open(97, file="porous.inp")
    read(97,*) n
    do i=1, n
      read(97,*,end=200) m, porous_eps_list(m), porous_dp_list(m), porous_hv_list(m)
      if(m .lt. 1 .or. m .gt. Total_block) then
        print*, "  porous.inp: block index out of range", m
        cycle
      endif
      if(porous_eps_list(m) .le. 0.d0 .or. porous_eps_list(m) .gt. 1.d0) then
        print*, "  porous.inp block", m, ": eps must be in (0,1], reset to 0.9"
        porous_eps_list(m)=0.9d0
      endif
      if(porous_dp_list(m) .le. 0.d0) then
        print*, "  porous.inp block", m, ": dp must be > 0, reset to 1e-3"
        porous_dp_list(m)=1.d-3
      endif
      cycle
200   continue
    enddo
    close(97)
    print*, "Read porous.inp OK"
    do m=1, Total_block
      if(Block_Type_List(m) .eq. BLOCK_POROUS) then
        print*, "  Porous Block", m, ": eps=", porous_eps_list(m), &
                " dp=", porous_dp_list(m), " hv=", porous_hv_list(m), &
                " rho_s=", solid_rho_list(m), " k_s=", solid_k_list(m)
      endif
    enddo
  end subroutine read_porous_prop


!----------------------------------------------------------------------
! Read solid thermal boundary conditions from solid_bc.inp
! Format:
!   Total_Solid_Blocks
!   Block_no
!   Num_physical_faces
!   face_no  Tw(K)  Qw(W/m2)
!   ...
! Tw>0: isothermal wall at Tw(K); Tw<0: heat flux wall with Qw(W/m2)
! Called after Mesh is allocated (Mesh(1)%Block pointers exist)
  subroutine read_solid_bc
   use Global_var
   use mod_type_def
   implicit none
   integer:: m, n, i, j, nf, unit, ios
   logical:: ex
   integer,allocatable:: nface_global(:)
   integer,allocatable:: face_global(:,:)
   real(PRE_EC),allocatable:: Tw_global(:,:), Qw_global(:,:)
   real(PRE_EC),allocatable:: htc_global(:,:), Tinf_global(:,:)
   real(PRE_EC):: f1, f2, f3
   character(len=256):: line
   Type (Block_TYPE),pointer:: B

   if(my_id .eq. 0) then
     inquire(file="solid_bc.inp", exist=ex)
     if(.not. ex) then
       print*, "solid_bc.inp not found, no solid thermal BCs"
       return
     endif
     open(newunit=unit, file="solid_bc.inp")
     read(unit,*) n
     allocate(nface_global(Total_block))
     nface_global = 0
     allocate(face_global(6, Total_block))
     allocate(Tw_global(6, Total_block))
     allocate(Qw_global(6, Total_block))
     allocate(htc_global(6, Total_block))
     allocate(Tinf_global(6, Total_block))
     face_global = 0; Tw_global = 0.d0; Qw_global = 0.d0
     htc_global = 0.d0; Tinf_global = 0.d0
     do i = 1, n
       read(unit,*) m
       read(unit,*) nf
       if(m .ge. 1 .and. m .le. Total_block) then
         nface_global(m) = nf
         do j = 1, nf
           read(unit,'(A)') line
           read(line,*,iostat=ios) face_global(j,m), Tw_global(j,m), Qw_global(j,m)
           if(ios .eq. 0) then
             read(line,*,iostat=ios) f1, f2, f3, htc_global(j,m), Tinf_global(j,m)
           endif
         enddo
       endif
     enddo
     close(unit)
     print*, "Read solid_bc.inp OK"
   endif

!  Broadcast to all processes
   if(.not. allocated(nface_global)) allocate(nface_global(Total_block))
   if(.not. allocated(face_global)) allocate(face_global(6, Total_block))
   if(.not. allocated(Tw_global)) allocate(Tw_global(6, Total_block))
   if(.not. allocated(Qw_global)) allocate(Qw_global(6, Total_block))
   if(.not. allocated(htc_global)) allocate(htc_global(6, Total_block))
   if(.not. allocated(Tinf_global)) allocate(Tinf_global(6, Total_block))
   call MPI_Bcast(nface_global, Total_block, MPI_INTEGER, 0, MPI_COMM_WORLD, i)
   call MPI_Bcast(face_global, 6*Total_block, MPI_INTEGER, 0, MPI_COMM_WORLD, i)
   call MPI_Bcast(Tw_global, 6*Total_block, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, i)
   call MPI_Bcast(Qw_global, 6*Total_block, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, i)
   call MPI_Bcast(htc_global, 6*Total_block, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, i)
   call MPI_Bcast(Tinf_global, 6*Total_block, OCFD_DATA_TYPE, 0, MPI_COMM_WORLD, i)

!  Set per-block data
   do m = 1, Total_block
     nf = nface_global(m)
     if(nf > 0) then
       B => Mesh(1)%Block(m)
       B%solid_bc_nface = nf
       allocate(B%solid_bc_face_no(nf))
       allocate(B%solid_bc_Tw(nf))
       allocate(B%solid_bc_Qw(nf))
       allocate(B%solid_bc_htc(nf))
       allocate(B%solid_bc_Tinf(nf))
       do j = 1, nf
         B%solid_bc_face_no(j) = face_global(j,m)
         B%solid_bc_Tw(j) = Tw_global(j,m)
         B%solid_bc_Qw(j) = Qw_global(j,m)
         B%solid_bc_htc(j) = htc_global(j,m)
         B%solid_bc_Tinf(j) = Tinf_global(j,m)
       enddo
       if(my_id .eq. 0) then
         print*, "  Block", m, " (solid):", nf, " thermal BC faces"
         do j = 1, nf
           print*, "    face", B%solid_bc_face_no(j), " Tw=", B%solid_bc_Tw(j), " Qw=", B%solid_bc_Qw(j)
         enddo
       endif
     endif
   enddo

   deallocate(nface_global, face_global, Tw_global, Qw_global, htc_global, Tinf_global)
  end subroutine read_solid_bc
