!=============================================================================
!  极限测试参照算例 (A) : run_case1_solid 的短调度副本
!    Block1 = 不锈钢固体 (BLOCK_SOLID), Block2 = Ma=3 空气
!    接口码 11 (solid j- <-> gas i=max), 交错 CHT 耦合
!  仅压缩调度, 物理设置与 run_case1_solid 完全一致:
!    Kstep_Couple_Comp=100, Niter_Couple_Outer=2, Niter_Couple_Warm=10
!    Kstep_save=100 -> Kstep=100/200 各写一次 flow3d_node.dat
!    Iflag_restart=-1 (不读不写重启, 保证与对照算例同一初场)
!  对照算例: limits/run_por_eps0 (多孔块 eps=0)
!=============================================================================
$freestream_ec
  Ma=3.0d0, Re=5000.d0, T_inf=800.d0, Lscale=0.001d0
$end

$flow_ec
  t_end=1.0d6, Kstep_save=100, Kstep_show=100,
  dt_global=1.0d0, dtmax=1000.d0,
  Mesh_File_Format=2, IFLAG_LIMIT_FLOW=0, Bound_Scheme=-1,
  Iflag_vtk_onefile=1, Iflag_vtk_SI=1,
  Iflag_restart=-1, Kstep_restart=0, Iflag_flow_node=1
$end

$lowspeed_ec
  !  本算例无低速块; 取值与源算例一致
  LS_rho=1.177d0, LS_mu=1.846d-5, LS_k=0.0262d0, LS_T_ref=300.d0,
  LS_U_in=0.d0, LS_P_in=0.d0, LS_P_out=0.d0,
  LS_alpha_p=0.05d0, LS_alpha_u=0.25d0,
  LS_Max_Iter=2500, LS_Tol=5.d-6, LS_Scheme=3, LS_Algorithm=3
$end

$ac_ec
  AC_Max_Iter=5000
$end

$porous_ec
  !  本算例无多孔块; 取值与源算例一致
  Porous_T_ref=300.d0, Porous_T_in=300.d0
$end

$solid_ec
  !  短调度: 压缩固体 GS 预算 (源算例用默认 1.7/20000/5/1e-9)
  Solid_Max_Iter=3000, Solid_Tol=1.d-8
$end

$couple_ec
  Iflag_Couple_Scheme=1, Kstep_Couple_Comp=100, Kstep_Couple_Min=100,
  Niter_Couple_Warm=10, Niter_Couple_Outer=2,
  Tol_Couple_Tw=1.0d-2
$end
