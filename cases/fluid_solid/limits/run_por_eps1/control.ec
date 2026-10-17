!=============================================================================
!  极限测试 (ε→1) : Block1 = 多孔壁 ε=1.0, Block2 = Ma=3 空气, 接口码 19
!    porous.inp : eps=1.0, dp=1e-3 m, hv=2e6 W/(m3 K)
!      -> K = dp^2*eps^3/(150(1-eps)^2) -> ∞, 阻力 -> 0 (纯 N-S)
!      -> k_f,eff = eps*LS_k = LS_k;  k_s,eff = (1-eps)*k_s = 0
!      -> hv 有限大: Ts 被拉到 Tf (无骨架 => 壁温 = 流体温度)
!  对照算例: limits/run_ref_case2 (低速 BLOCK_LOWSPEED, 同几何/同入口)
!  说明: 接口码 19 的气体侧是"无滑移等温壁 + 可选发汗"(多孔壁模型), 而码 12
!        的气体侧是"流体-流体连续界面", 二者界面模型不同, 因此本算例的
!        对照目标是"块内流场(动量/连续方程)一致性 + 差值集中于界面一排".
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
  !  Block1 冷却剂空气 300 K / 1 atm; 质量入口 G=2 kg/(m2 s) -> 1.699 m/s (与 ref_case2 相同)
  LS_rho=1.177d0, LS_mu=1.846d-5, LS_k=0.0262d0, LS_T_ref=300.d0,
  LS_U_in=0.d0, LS_V_in=1.699d0,
  LS_P_in=0.d0, LS_P_out=0.d0,
  LS_alpha_p=0.2d0, LS_alpha_u=0.5d0,
  LS_Max_Iter=150, LS_Tol=1.d-7, LS_Algorithm=3
$end

$ac_ec
  !  与 ref_case2 完全相同的 AC 设置 (eps=1 时多孔块方程 = 低速 AC 方程)
  AC_Max_Iter=5000, AC_Print=2500, AC_beta=10.d0, AC_CFL=2.d0, AC_Tol=1.d-6
$end

$porous_ec
  !  eps=1 极限: 骨架消失 (见 porous.inp); hv 保留以给出适定的壁温
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
