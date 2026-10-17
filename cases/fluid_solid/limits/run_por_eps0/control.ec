!=============================================================================
!  极限测试 (精确 ε=0, 需打补丁) : Block1 = 多孔壁 ε=0 (无孔隙骨架), Block2 = Ma=3
!    porous.inp : eps=0.0, dp=1e-3 m, hv=2e6 W/(m3 K)
!      -> k_s,eff = (1-eps)*k_s = k_s => 纯导热 (等效固体)
!      -> 需要补丁: (1) eps 允许 0; (2) eps<=1e-8 时冻结速度 (U(2:4)=0)
!   几何/边界与本目录 bc3d.inp 完全对齐 case1:
!      块1 j- = 19 (多孔界面), 块1 j+ = 2 (等温壁, solid_bc.inp 面 5 = 300 K)
!  对照算例: limits/run_ref_case1 (固体 BLOCK_SOLID, 接口码 11)
!  说明: 码 19 的骨架侧用"热流 q_w"边界, 码 11 的固体侧用"界面温度 T_i"
!        边界, 二者界面的离散形式不同, 只能达到物理等效 (同一收敛解) 而非
!        逐位相同; 本算例的对照目标是 Ts 场与 case1 的 T 场在容差内一致.
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
  !  块1 j+ 是等温壁 (无冷却剂入口) => 入口速度取 0
  LS_rho=1.177d0, LS_mu=1.846d-5, LS_k=0.0262d0, LS_T_ref=300.d0,
  LS_U_in=0.d0, LS_V_in=0.d0,
  LS_P_in=0.d0, LS_P_out=0.d0,
  LS_alpha_p=0.2d0, LS_alpha_u=0.5d0,
  LS_Max_Iter=150, LS_Tol=1.d-7, LS_Algorithm=3
$end

$ac_ec
  !  eps=0 时速度被冻结, AC 实际只做几步; 仍与参照算例保持同一设置
  AC_Max_Iter=5000, AC_Print=2500, AC_beta=10.d0, AC_CFL=2.d0, AC_Tol=1.d-6
$end

$porous_ec
  !  eps=0 极限: hv>0 使 Tf ≡ Ts (单温度), 骨架方程退化为纯导热
  Porous_T_ref=300.d0, Porous_alpha_Ts=1.0d0, Porous_Max_Iter=3000,
  Porous_V_in=0.d0, Porous_T_in=300.d0
$end

$solid_ec
  !  与 ref_case1 相同的固体 GS 预算
  Solid_Max_Iter=3000, Solid_Tol=1.d-8
$end

$couple_ec
  Iflag_Couple_Scheme=1, Kstep_Couple_Comp=100, Kstep_Couple_Min=100,
  Niter_Couple_Warm=10, Niter_Couple_Outer=2, Porous_Chunk_Iter=2
$end
