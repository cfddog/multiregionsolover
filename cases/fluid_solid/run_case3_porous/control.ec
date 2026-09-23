!=============================================================================
!  Case 3 : Block1 = stainless-steel porous wall (eps=0.30), Block2 = Ma=3 air
!           interface code 19 (porous j=1 <-> gas i=max), staggered porous
!           coupling (gas chunk at isothermal wall T_w -> q_w -> porous chunk
!           with coolant injection -> back to T_w).
!           Block1 j=max face = coolant mass inlet, 2 kg/(m^2 s) at 300 K.
!           The porous solver is constant-density (LS_rho=1.177 kg/m3), so
!           G=2 kg/(m^2 s) is imposed as
!             LS_V_in = Porous_V_in = G/LS_rho = 1.699 m/s   (in +y)
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
  !  AC 伪时间求解器（多孔块，LS_Algorithm=3）；预算已压缩：5000×3 ≈170 s/轮
  AC_Max_Iter=5000, AC_Print=1000, AC_beta=100.d0, AC_CFL=20.d0, AC_Tol=1.d-6
$end

$porous_ec
  !  不锈钢多孔骨架 eps=0.30（eps/dp/hv 见 porous.inp）
  Porous_T_ref=300.d0, Porous_alpha_Ts=1.0d0, Porous_Max_Iter=30000,
  Porous_V_in=1.699d0, Porous_T_in=300.d0
$end

$couple_ec
  Iflag_Couple_Scheme=1, Niter_Couple_Warm=2, Niter_Couple_Outer=12,
  Porous_Chunk_Iter=3
$end
