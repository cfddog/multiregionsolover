!=============================================================================
!  极限测试 (ε→0 逼近, 不打补丁) : Block1 = 多孔壁 ε=1.0e-6, Block2 = Ma=3 空气
!    porous.inp : eps=1e-6, dp=1e-3 m, hv=2e6 W/(m3 K)
!      -> 阻力 c_drag 中的 eps^2*mu/K 项: K 被下限 1e-30 截住 =>
!         阻力 ~1e12*mu 级别 => 速度被压到 ~1e-9 m/s (事实上冻结)
!      -> k_s,eff = (1-eps)*k_s ≈ k_s => 骨架为纯导热 (等效固体)
!  对照算例: limits/run_ref_case1 (固体 BLOCK_SOLID; 块的 j+ 面已对齐为
!            bc2 + solid_bc.inp 面 5 = 300 K)
!  作用: 在不改代码的前提下, 验证 "eps→0 时多孔算例退化为固体算例".
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
  !  Block1 冷却剂空气 300 K / 1 atm; 入口 G=2 kg/(m2 s) -> 1.699 m/s (与源算例相同)
  LS_rho=1.177d0, LS_mu=1.846d-5, LS_k=0.0262d0, LS_T_ref=300.d0,
  LS_U_in=0.d0, LS_V_in=1.699d0,
  LS_P_in=0.d0, LS_P_out=0.d0,
  LS_alpha_p=0.2d0, LS_alpha_u=0.5d0,
  LS_Max_Iter=150, LS_Tol=1.d-7, LS_Algorithm=3
$end

$ac_ec
  !  与 ref_case1/ref_case2 一致的 AC 预算 (便于逐项对比)
  AC_Max_Iter=5000, AC_Print=2500, AC_beta=10.d0, AC_CFL=2.d0, AC_Tol=1.d-6
$end

$porous_ec
  Porous_T_ref=300.d0, Porous_alpha_Ts=1.0d0, Porous_Max_Iter=3000,
  Porous_V_in=1.699d0, Porous_T_in=300.d0
$end

$solid_ec
  Solid_Max_Iter=3000, Solid_Tol=1.d-8
$end

$couple_ec
  Iflag_Couple_Scheme=1, Kstep_Couple_Comp=100, Kstep_Couple_Min=100,
  Niter_Couple_Warm=10, Niter_Couple_Outer=2, Porous_Chunk_Iter=2
$end
