# Worklog — 简要记录（最新在最后）

> 本文件只记“简要条目”。**详细开发历史**见 `docs/工作日志.md`（830+ 行，权威）；
> 规则要求每次任务开始时只读本文件**最近 5 条**。

---

## 2026-09-18 — cases/fluid_solid 三算例（11/12/19，gas i+ ↔ 对侧 j-）分段交错

**目标**：在 `cases/fluid_solid` 网格上做三个分段交错耦合算例（先流体 1000 步 → 再跑另一侧 →
逐步减半间隔 → 最后 `Kstep_Couple_Min=1`）：①固体↔可压(11) ②低速↔可压(12) ③多孔↔可压(19)；
后两者的 Block1 j=max 为冷却剂入口 G=2 kg/(m²·s)、300 K；Block2 = Ma=3/Re=5000/T∞=800K、壁面绝热。

**难点**：既有耦合例程硬编码 `gas j- ↔ 对侧 j+`，而本网格是 `gas i+ (4) ↔ 对侧 j- (2)`
（实测 `L=(2,-1,3)`，切向反向）。为不改网格、保持“Block1 j=1 为对接面”，选择**扩展代码**。

**改动**（`src/`，向后兼容，未提交）：
- `sub_multi_region.f90`：`couple_highlow_fluid_face` / `couple_compressible_porous_blowing_face`
  新增 `(4,2)` 分支；`run_staggered_fluid_solid/highlow/fluid_porous` 支持该拓扑（porous 驱动
  按 `pair42` 重构：壁 ghost、q_w 提取、多孔 j- ghost、`set_porous_Ts_flux` 加方向参数）；
  新增发汗壁 `wall_bound_blowing`（`v_w=Bp%U(3,gi,jc_l,k)`）与界面诊断 `dbg_dump_interfaces`。
- `opencfd_ec3d_v1.16a.f90`：`IF_Debug==1` 时调用界面 dump。

**产物**：`cases/fluid_solid/{run_case1_solid,run_case2_lowspeed,run_case3_porous}/`；
`docs/进度_fluid_solid三算例.md`；`.clinerules` + `memory-bank/`。

**验证**：编译 EXIT=0；三例短调度试跑 EXIT=0、无 NaN；交付版参数（`Kstep_Couple_Comp=1000`）
首轮调度正确打印。case1 `max|dT_w|` 499.9→0.094 K、`T_w≈799.9 K`、`max|q_w|`≈2.1e3 W/m²；
case3 `max|q_w|` 5.3e6→3.9e6 W/m²、`T_w` 309→330 K。

**遗留/下一步**：三例 README；case2 的 AC 标定（`res_q~1e-2`）；case3 外层加速（多孔储热）；
清理产物；文档同步；提交。**必须单进程**（12/19 未实现跨进程）。

---

## 2026-09-18（第二轮）— 三算例现场测试 + 设置修正

**目标**：确认 `run_case1_solid / run_case2_lowspeed / run_case3_porous` 在**交付版 control.ec**
下能跑、设置无问题（不要求跑完）。

**做法**：`/tmp/smoke_all.sh` 依次 `make` + `mpirun -np 1`（每例 420 s 限时），日志各目录
`run_smoke.log`；另在 `/tmp/fs_t2a|t2b|t3` 做配置对比试验（`/tmp/tune2.sh`、`/tmp/tune3.sh`）。

**结果**：
- case1 ✅ **12 轮全跑完**，EXIT=0，**99 s**，无 NaN；调度 1000,1000,500,250,125,62,31,15,7,…
  正确；`max|dT_w|` 499.9→0.0996 K（平台）、`T_w` 799.90→799.20 K、`max|q_w|` 2064→1672 W/m²。
  `Tol_Couple_Tw=1e-2` 未达因**固体未达稳态**（T_w 每轮缓降 ~0.1 K）——物理瞬态，非设置错误。
- case2 ✅ 正常推进/无 NaN（420 s 限时 EXIT=124）；界面度量 `max|dT|` 321→20.2→17.5→9.2 K、
  `max|dp|`0.2–0.9 Pa、`max|du|`29→3.0 m/s。LS-AC 每轮跑满 5000 迭代、`res_q` 平台 ~1.4e-2。
- case3 ✅ 正常推进/无 NaN（437 s 限时 EXIT=124）；首轮气体段 `max|q_w|`=1.09e6 W/m²、
  多孔 AC `res_q`=5.6e-5、LTNE `dTs`→9.9e-4；但**原设置 `AC_Max_Iter=20000 × Porous_Chunk_Iter=10`
  首轮多孔段 437 s 跑不完**（≈65 min/轮）。

**设置修正（已改交付文件）**：
1. case3 `control.ec`：`AC_Max_Iter 20000→5000`、`AC_Print 5000→1000`、`AC_CFL 10→20`、
   `Porous_Chunk_Iter 10→3`。试验：**2 轮完成/无 NaN/≈170 s/轮**，`res_q`=6.4e-5（与原相当）。
   复验 `run_check_new.log`：`pair42=T`、1000 步气体段、多孔 AC 3×5001 迭代 res_q 9.7e-5→6.4e-5、无 NaN。
2. case2 **不改**（保留 `AC_CFL=2, AC_beta=10, AC_Max_Iter=5000`）：试 `AC_CFL=20/beta=100`
   虽把 `res_q` 降到 3.3e-3 但**第 3 轮 NaN**；试 `LS_Algorithm=1`(SIMPLE) **第 1 步即发散**
   （res_p→6e31, u→-1e6）。⇒ 必须 AC + CFL=2，`res_q` 平台为已知限制。
3. case1 **不改**：如需外层早停可把 `Tol_Couple_Tw` 放宽到 ~0.2 K。

**遗留**：三例 README；清理算例目录运行产物；文档同步；提交。单机预估：case1 ~100 s、
case2 ~20 min、case3 ~40 min（12 轮，单进程）。

---

## 2026-09-23 — control.ec 拆分为 7 组 namelist（兼容旧单组，零行为变化）

**目标**：`$control_ec` 一个组塞 128 个变量不便维护（实测 70~75% 的赋值行与默认值相同）；
按物理/功能域拆成多组，且**不改任何现有算例文件**。

**盘点**：`$control_ec` = **128 个变量，全部有代码默认值**；代表算例文件行数/赋值数/与默认值相同数 =
case1 121/112/75、case2 128/114/75、case3 129/113/69、channel_ac 93/91/65、lidcavity 89/87/66。

**改动**（`src/`、`util/`）：
- `sub_modules.f90`：新增 `Solid_GS_Omega/Solid_Max_Iter/Solid_Min_Iter/Solid_Tol`（原硬编码于
  `solid_solver_one_block`，默认值全等 ⇒ 行为不变）。
- `sub_read_parameter.f90`：1 组 → **7 组**（`freestream 24 / flow 49 / lowspeed 21 / ac 15 /
  solid 4 / porous 8 / couple 11` = 128+4）；旧组 `$control_ec` 完整保留（+`Solid_*`）；
  新读取逻辑 = **扫描组头**（首个非空字符须为 `$`/`&`；`!` 后为注释）→ 旧组 → 7 组可选读，
  **新组覆盖旧组**；组存在但解析失败 ⇒ 报错 + `stop 1`；`output_para.out` 末尾回显各组来源；
  新增 `scan_control_ec_groups/str_to_lower/check_nml_ios`。
- `sub_multi_region.f90`：固体 GS 常量改读全局量。
- `util/readflow3d-ver2.4a/2.5.f90`：同步双读；并补齐 `Lscale` + `LS_*/AC_*/Porous_*/couple/Solid_*`
  声明（**历史遗留 bug：这两个工具原先连 `Lscale` 都缺，读现代算例必报 iostat=5010**）。
- 修一个既有潜在 bug：`read_parameter_ec` 的 `a0/d0/mu0/mu1` 非叶轮机械模式下未初始化即写入
  `output_para.out`（garbage、逐次不同）→ 显式清零。
- 新增 `control.ec.template`、`docs/control.ec-说明.md`。

**验证**：`make` EXIT=0 无告警；`cases/solid_1d`、`cases/channel_ac` 旧↔新二进制 A/B
`output_para.out` **除新增 4 行回显外逐字节一致**；新格式取值/新组覆盖旧组/严格报错(exit=1) 全通过；
模板 7 组同时存在可解析且 `Solid_Max_Iter=30` 生效（日志 GS iterations 31）；两个 util 工具
（legacy + 新格式）EXIT=0。

**遗留/下一步**：`cases/**/control.ec` 未迁移（第1档决定，靠双读兼容）；可选
`util/split_control_ec.py`；随后做**重启文件 `field_restart.dat`** 与**节点/SI 流场输出**
（新开关放 `$flow_ec`）。

---

## 2026-09-23（第二轮）— 重启文件 `field_restart.dat` + 节点/SI 流场输出

**目标**：(1) 流场变量 + `Ts` + **缓冲区** 写入重启文件，含步数/时间、定期保存、存在即自动续算；
(2) `flow3d.dat` → 网格**节点**、**SI 单位**、**unformatted**，可与 `Mesh3d.x` 组合显示。

**改动**：
- 新增 `src/sub_restart.f90`：`write_restart/read_restart/restart_step_output`。文件含头记录
  （版本/块数/NVAR/Kstep/LAP + tt/Ma/Re/gamma/T_inf/Lscale/dt_global）与逐块 `nx,ny,nz` +
  **U(1:NVAR,1-LAP:nx+LAP-1,…)**（含 4 层鬼点缓冲）+ `Ts` + `Tsn` + `p`；读入后刷新块间缓冲。
- `src/sub_IO.f90`：`output_flow_node`（+`flow_node_one_block`）：8 格心平均→节点、SI、unformatted、
  逐块一条记录 `d,u,v,w,T`、块序同 `Mesh3d.x`；固体块 `d=u=v=w=0, T=Ts[K]`。
- 开关 `Iflag_restart(0/1/-1)`、`Kstep_restart(<=0→Kstep_save)`、`Iflag_flow_node`：
  `$flow_ec` + 旧组 + `Ipara(55..57)` 广播 + `output_para.out` 回显 + 模板同步。
- 挂接：`sub_init.f90::Init_flow`（重启优先）、主循环、三个交错驱动（`tt/Kstep` 恢复 + 保存点调用 +
  结束写出；补丁脚本 3+3+3，修掉 `contains` 位置问题）；`makefile` +`sub_restart.f90`；
  `.gitignore` 忽略两个输出文件。

**关键 bug 修复**：`output_flow_node`/`read_restart` 的全局块循环在**非 0 号进程**上也执行并访问
只在本进程有效的 `B_n(m)` ⇒ `-np>1` 越界/挂死；改为「主机走全局块、从机只走自有块」
（与 `output_flow`/`read_flow_data` 同构）。同时**必须**把新开关加入 `bcast_para`，否则各 rank
取值不一致会在 I/O 处死锁。

**验证**：`make` EXIT=0；`solid_1d` 续算往返（Kstep 11→20、无 NaN，`restart from ... Kstep= 11`）；
`flow3d_node.dat` 维度与 `Mesh3d.x` 完全一致（4131 / 13534 / 40602 点）、固体块 T 为 K 量级；
**`-np 2` 2 块**（`run_case1_solid` 改 `Iflag_Couple_Scheme=0`：固体块与可压块分处两个 rank）
写出并续算成功，SI 数值合理（Ma=3 → ρ≈0.10–0.13 kg/m³、u≈788–1701 m/s、T≈800–1174 K）；
`-np 2` 空闲 rank 不死锁；交错 `Iflag_Couple_Scheme=1` 的 case1 在 Kstep=50/100/150 写出两文件、无 NaN。
（注：`channel_ac_2blk` 等低速 AC 算例在 `-np 2` 下**旧二进制也失败**，属既有问题，与本改动无关。）

**遗留**：未改任何 `cases/**/control.ec`（`Iflag_flow_node` 默认 0；需要时加一行即可）；
交错界面量 `fp_*` 仍从初值重播（可选后续写入重启文件以做到零跳变）。

---

## 2026-09-24 — fluid_solid 三算例 control.ec 迁移到 7 组写法 + 重启实测

**目标**：把考核 case `fluid_solid/run_case1_solid|run_case2_lowspeed|run_case3_porous` 的
`control.ec` 改写成新 7 组写法，并实测**重启文件可用**。

**做法**：
- 先用脚本算出三个文件里**与代码默认值不同的项**（case1 29+13、case2 29+14、case3 35+13 项），
  按组重写 `control.ec`（只留非默认项），原文件备份为同目录 `control.ec.legacy`；
  `$flow_ec` 里显式写 `Iflag_restart=0, Kstep_restart=0, Iflag_flow_node=1`；
  `$solid_ec` 省略（默认 = 原硬编码 1.7/20000/5/1e-9）。
- 挂接的三处硬骨头：case2 `Tol_Couple_Tw=2.0/Tol_Couple_p=100/Tol_Couple_u=0.1`、
  case3 `AC_CFL=20/AC_beta=100/AC_Print=1000/Porous_Chunk_Iter=3` 等**必须保留**的“非默认项”全部保留。

**验证（`mpirun -np 1`）**：
- **等价性 A/B**：旧 `control.ec.legacy` 与新 `control.ec` 各短跑，`output_para.out` 过滤掉
  5 行回显后**逐行完全一致（42 行）** —— 三例全部 ✅
- **重启往返**：case1 写 Kstep=2500/2995（run1 整段跑完 exit=0）、case2 写 2000/2500、
  case3 写 500/1000；run2 均打印 `restart from … Kstep=…` 与
  `restart: continue staggered run at Kstep=…`，NaN=0 ✅
- **节点文件**：三例 `flow3d_node.dat` 均 2 记录、`(67,101,2)+(101,201,2)`、点数 13534/40602
  与 `Mesh3d.x` 一致；固体块 `d=u=v=w=0`、T∈[300,799.9] K；可压块 SI 合理
  （case1 ρ=0.041–0.115 kg/m³、u=38.9–1706 m/s、T=793–2106 K）✅

**遗留**：交错驱动续算是“状态恢复 + 外层轮数重新计”，run2 仍跑满 12 轮（设计行为）；
其余算例目录仍为旧写法（靠双读兼容），可按需用同一方法逐个迁移。

---

## 2026-09-24 — 提交（commit 29fc24f）

**提交**：`29fc24f` on `main`（44 files, +3772/-441），消息
`feat(io): 重启文件 field_restart.dat + 节点/SI 流场输出; refactor(control.ec): 7 组 namelist`。

**内容**：control.ec 7 组拆分（旧组双读兼容）+ `src/sub_restart.f90`（含 LAP=4 缓冲的重启读写、
自动续算）+ `output_flow_node`（节点/SI/unformatted）+ 交错驱动续算 + 三算例 `control.ec` 迁移
（备份 `.legacy`）+ `util/readflow3d-ver2.4a/2.5` 双读与声明补齐 + 文档/记忆库；
并含 2026-09-18 会话遗留的未提交改动（码 11/12/19 的 `gas i+ ↔ 对侧 j-` 拓扑等）。

**提交前验证**：`make clean && make` **EXIT=0**（新/改文件零告警；29 条告警全来自未触碰的
`sub_Finite_Difference1/2.f90`）；旧↔新二进制 `output_para.out` 除新增回显行外**逐字节一致**；
`solid_1d` 续算往返 Kstep 11→20；`-np 2` 双块写出并续算；交错 case1/2/3 重启往返 **NaN=0**；
`flow3d_node.dat` 点数与 `Mesh3d.x` 完全一致。

**未提交/遗留**：`.vscode/settings.json`（主题，有意保留）；各算例运行产物（`flow3d.dat`、
`*.vtk`、日志）与 `sample_code/`、`workplan.md`、`src/*.log` 等仍为 untracked；
`origin/main` 落后 1 个提交（未 push）。

---

<!-- 新条目请追加在下面（格式：## YYYY-MM-DD — 标题 / 目标 / 改动 / 验证 / 遗留） -->

