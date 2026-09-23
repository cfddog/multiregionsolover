!=============================================================================
!  Case 2 : Block1 = low-speed air (BLOCK_LOWSPEED, coolant side),
!           Block2 = Ma=3 air; interface code 12 (LS j=1 <-> gas i=max).
!           Block1 j=max face = coolant mass inlet, 2 kg/(m^2 s) at 300 K.
!           The low-speed solver is constant-density (LS_rho=1.177 kg/m3), so
!           the mass flux is imposed as an equivalent normal velocity
!             LS_V_in = G/LS_rho = 2/1.177 = 1.699 m/s   (in +y)
!           (the LS_Inlet_Type=2 mass path assumes an i-face inlet in mesh
!            units, so the velocity inlet is used here.)
!
!  control.ec 已改为 7 组 namelist（2026-09-23）：只写与默认值不同的项。
!  原文件备份为 control.ec.legacy。$solid_ec 用默认 1.7/20000/5/1e-9。
!=============================================================================
$freestream_ec
  Ma=3.0d0, Re=5000.d0, T_inf=800.d0, Lscale=0.001d0
  !  AoA=0, Twall=-1(绝热), p_outlet=-1, Ref_S/Ref_L=1 取默认
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
  !  Block1 冷却剂空气 300 K / 1 atm；质量入口 G=2 用等效速度 1.699 m/s
  LS_rho=1.177d0, LS_mu=1.846d-5, LS_k=0.0262d0, LS_T_ref=300.d0,
  LS_U_in=0.d0, LS_V_in=1.699d0,
  LS_P_in=0.d0, LS_P_out=0.d0,
  LS_alpha_p=0.2d0, LS_alpha_u=0.5d0,
  LS_Max_Iter=150, LS_Tol=1.d-7, LS_Algorithm=3
$end

$ac_ec
  !  AC 伪时间求解器（LS_Algorithm=3）；本算例标定：必须 AC_CFL=2 / AC_beta=10
  AC_Max_Iter=5000
$end

$porous_ec
  !  本算例无多孔块；取值与原 control.ec 一致
  Porous_T_ref=300.d0, Porous_T_in=300.d0
$end

$couple_ec
  Iflag_Couple_Scheme=1, Niter_Couple_Warm=2, Niter_Couple_Outer=12,
  !  流-流(码12)外层判据（本算例整体放宽）
  Tol_Couple_Tw=2.0d0, Tol_Couple_p=1.0d2, Tol_Couple_u=1.0d-1
$end
