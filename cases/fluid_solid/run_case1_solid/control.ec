!=============================================================================
!  Case 1 : Block1 = stainless-steel solid (conjugate wall), Block2 = Ma=3 air
!           interface code 11 (solid j=1 <-> gas i=max), staggered CHT coupling
!
!  control.ec 已改为 7 组 namelist（2026-09-23）：只写与代码默认值不同的项，
!  每组/每变量都可省略；逐块逐面数据仍在 material.in / solid_bc.inp / porous.inp。
!  旧式单组 $control_ec 仍兼容，但本文件已迁移（原文件备份为 control.ec.legacy）。
!  固体 GS 控制 $solid_ec 本算例取默认 1.7 / 20000 / 5 / 1e-9（=原硬编码值）。
!=============================================================================
$freestream_ec
  Ma=3.0d0, Re=5000.d0, T_inf=800.d0, Lscale=0.001d0
  !  AoA=0, Twall=-1(绝热), p_outlet=-1(外插), Ref_S/Ref_L=1, Cood_Y_UP=1 取默认
$end

$flow_ec
  t_end=1.0d6, Kstep_save=500, Kstep_show=50,
  dt_global=1.0d0, dtmax=1000.d0,
  Mesh_File_Format=2, IFLAG_LIMIT_FLOW=0, Bound_Scheme=-1,
  Iflag_vtk_onefile=1, Iflag_vtk_SI=1,
  !---- 重启与节点流场输出（新增）----
  Iflag_restart=0,        ! 存在 field_restart.dat 即自动续算；-1=完全关闭
  Kstep_restart=0,        ! 重启写出间隔(<=0 时 = Kstep_save)
  Iflag_flow_node=1       ! 写节点/SI 流场 flow3d_node.dat（配 Mesh3d.x 显示）
$end

$lowspeed_ec
  !  本算例无低速块；以下取值与原 control.ec 完全一致（便于逐项比对）
  LS_rho=1.177d0, LS_mu=1.846d-5, LS_k=0.0262d0, LS_T_ref=300.d0,
  LS_U_in=0.d0, LS_P_in=0.d0, LS_P_out=0.d0,
  LS_alpha_p=0.05d0, LS_alpha_u=0.25d0,
  LS_Max_Iter=2500, LS_Tol=5.d-6, LS_Scheme=3, LS_Algorithm=3
$end

$ac_ec
  AC_Max_Iter=5000
$end

$porous_ec
  !  本算例无多孔块；取值与原 control.ec 一致
  Porous_T_ref=300.d0, Porous_T_in=300.d0
$end

$couple_ec
  Iflag_Couple_Scheme=1, Niter_Couple_Warm=2, Niter_Couple_Outer=12,
  Tol_Couple_Tw=1.0d-2
$end
