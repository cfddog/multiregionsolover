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
已 `git push origin main`（`adaa43b..6e15bcd`，远端与本地同步）。

---

## 2026-09-24（补）— flow3d_node.dat 补上固体/骨架温度 Ts（5→6 变量）

**问题（用户指出）**：`src/sub_IO.f90::flow_node_one_block` 只填一个温度槽位 ——
固体块填 `B%Ts`、多孔块填 `U(5)`（流体温度）⇒ **多孔块骨架温度 `B%Ts` 没进节点文件**
（固体块虽"有"温度但占用的是 `T` 槽位，也没有独立的 `Ts`）。LTNE 温差因此不可见。

**改动**：每条记录 `d,u,v,w,T` → **`d,u,v,w,T,Ts`**（6 变量，逐块同长，块序仍同 `Mesh3d.x`）：
- `BLOCK_SOLID`：`T = Ts = B%Ts`（`T` 槽位保持 `flow3d.vtk` 约定不变）
- `BLOCK_POROUS`：`T = U(5)`（流体温度 Tf）、`Ts = B%Ts`（骨架）← 新增信息
- `BLOCK_LOWSPEED` / `BLOCK_FLOW`：`T = Tf`、`Ts = T`（无固相，镜像以统一布局）
同步 `output_flow_node` 的 `Num_data/allocate/write`；新增 `util/check_flow3d_node.py`。

**验证**：`make` EXIT=0；`run_case1_solid` 固体块 `T=Ts=300–800 K`、可压块 `T=Ts`；
**`porous_ltne_1d`（q=2e6 W/m²、hv=2e6 W/(m³K)）→ Tf=300.7–1039 K、Ts=407.7–1174 K、
max|T−Ts|=157.7 K**（修复前该值不在文件内）；`run_case3_porous` 多孔块两温度分列；
三例 `STRUCTURE OK`、NaN=0。

**兼容性**：`flow3d_node.dat` 是本次新增文件，尚无外部读者；读取方需按 6 变量解析
（`util/check_flow3d_node.py` 已按 6 实现）。

---

<!-- 新条目请追加在下面（格式：## YYYY-MM-DD — 标题 / 目标 / 改动 / 验证 / 遗留） -->

---

## 2026-09-25 — 修复集群 run1 启动崩溃：control.ec 同一文件连两个 unit

**现象**（集群 `/work/home/lijunyang/sundong/PorousTest/code/`，`run_case1_solid/./OPENCFDEC`）：
`At line 976 of file sub_read_parameter.f90 (unit = 98)` →
`Fortran runtime error: File already opened in another unit`；
backtrace：`scan_control_ec_groups_ ← read_parameter_ec_ ← read_parameter_ ← MAIN__`。

**根因**：`read_parameter_ec` 第 293 行已把 `control.ec` 连到 **unit 99**（为读 7 组 namelist），
而 `scan_control_ec_groups` 第 976 行又 `open(98,file="control.ec",status='old')` 打开同一文件
⇒ **同一文件同时连到两个 unit**，违反 Fortran 标准（F2018 §12.5.6）。本地 gfortran 13.3.0
用最小程序实测同一文件开两个 unit 返回 `iostat=0`（新版 libgfortran 放宽了检查）⇒
2026-09-23/24 的本地验证全部通过；集群的旧 libgfortran 按标准致命报错，且该 open **未带 iostat**
⇒ 运行时直接 `Error termination`（不是可捕获的 ios）。**与算例/参数/并行无关，纯代码可移植性缺陷。**

**改动**（仅 `src/sub_read_parameter.f90`，1 file / +20 / −5）：
- `scan_control_ec_groups(unit, has_legacy, ...)`：新增 `integer,intent(in):: unit`，
  内部 `rewind(unit)` + `read(unit,'(A)',iostat=ios)`；**删除 `open(98,...)` 与 `close(98)`**
  （不 close：调用者随后还要用同一 unit 逐组 `rewind(99); read(99,nml=...)`）。
- 调用点 → `call scan_control_ec_groups(99, ...)`。
- 注释（中英）写明“为什么不能自己 open”：这是可移植性约束，不是风格问题。

**验证**：
- `cd src && make` **EXIT=0**（只重编 `sub_read_parameter.o` 并链接；make 输出**零 warning/error**）。
- 静态审计：`grep` 全 `src/*.f90` ⇒ `control.ec` 只剩 **1 处 open**（第 293 行 unit 99）；unit 98 不再使用。
- 功能等价：把 `cases/fluid_solid/run_case1_solid` 输入 + 新二进制放到 `/tmp/val_case1`，
  `mpirun -np 1` **整段跑完 12 轮**（`max|dT_w|` 499.9 → 0.0995 K，与历史一致），写出
  `field_restart.dat`(Kstep=2995) 与 `flow3d_node.dat`(6 变量)，**NaN=0**；
  `output_para.out` 与改动前二进制产物**逐字节一致**（组回显仍为
  `legacy=F freestream=T flow=T lowspeed=T ac=T solid=F porous=T couple=T`）。
- 全仓库 open 审计：`sub_convert_inp.f90`(88/96/99/87 成对 close)、`sub_IO.f90`、`porous.inp@97`、
  `solid_bc.inp`（`newunit=`）等均无“同文件双 unit”；`open(99,file="output_para.out")` 属
  “同一 unit 换文件”，标准允许。

**遗留**：① 两个 util 工具有同样缺陷（`readflow3d-ver2.5.f90`：99@707 + 96@791；
`readflow3d-ver2.4a.f90`：99@681 + 96@765），按用户要求**本轮不修**，下次按同法修；
② 集群侧需同步 `src/sub_read_parameter.f90` 并重新 `make` 才生效；
③ 本次改动**尚未提交**。


---

## 2026-09-25（补）— FAQ：case1 为何只算到 2995 步 / 如何继续（续算实测）

**问**：`fluid_solid/run_case1` 只算到 Kstep≈2950–2995 就结束，不能往下算？

**答（设计行为，非卡死/发散）**：
- case1 `$couple_ec` 里 `Iflag_Couple_Scheme=1` ⇒ 主程序在 `opencfd_ec3d_v1.16a.f90:169-174`
  调 `run_staggered_multiregion(1)` 后 **`mpi_finalize; stop`**，**不进入** `do while(tt<t_end)`
  ⇒ `t_end=1e6` 在交错模式下**无效**。
- 气体步数 = 调度表（`sub_multi_region.f90:2567-2574`）：前 `Niter_Couple_Warm` 轮满
  `Kstep_Couple_Comp`，之后每次减半（下限 `Kstep_Couple_Min`）。case1 取
  `Warm=2, Outer=12`，`Comp/Min` 用默认 **1000/1** ⇒ 1000,1000,500,250,125,62,31,15,7,3,1,1
  = **2995 步**（日志 `outer iter 12 : flow chunk steps = 1`、`flow chunk done, Kstep= 2995`）。
- 终止判定 `sub_multi_region.f90:2673-2684`：`it>=2 且 max|dT_w|<Tol_Couple_Tw(=1e-2)` 才提前收敛；
  case1 `max|dT_w|` 平台 ≈0.0995 K（固体未达稳态）⇒ 永不满足，跑满 12 轮打印
  `reached Niter_Couple_Outer= 12 (not fully converged)`，随后写 `flow3d.dat/vtk/field_restart.dat`
  → **EXIT=0**。

**继续往下算的两条路**：
1. **续算（推荐，不改参数）**：再跑一次同一命令；`Iflag_restart=0` ⇒ 检测到 `field_restart.dat`
   自动恢复并累加 Kstep。**实测（`/tmp/val_case1`，新二进制）**：
   `restart from 'field_restart.dat': Kstep= 2995` → `continue staggered run at Kstep= 2995`
   → 3995 → 4995 → 5495 → … → **5990**（+2995 步/次），NaN=0，结束又写 Kstep=5990 重启文件。
2. **改调度**：`$couple_ec` 里加 `Kstep_Couple_Comp/Niter_Couple_Warm/Kstep_Couple_Min`
   并调大 `Niter_Couple_Outer`（默认 60）；建议 `Kstep_Couple_Min` 别用 1（末段几十轮都是 1 步，
   几乎不推进物理时间）。收敛判据在此算例达不到 ⇒ **外层轮数是唯一的时间推进控制**（见 constraints §5）。

**“不能继续”的三种典型原因**（排查顺序）：① 集群仍是**修复前二进制** ⇒ 又报
`File already opened in another unit`（需同步 `src/sub_read_parameter.f90` + `make`）；
② run1 未写到重启文件就被作业时限/调度 kill，或脚本清理目录 ⇒ `read_restart: not found`，从头开始；
③ 只跑 12 轮 —— 即上述设计行为。


---

## 2026-09-25（第二轮）— 重启后按耦合状态续算：继续交错 / 自动切逐步强耦合

**需求（用户）**：重启后不要重新播一遍交错调度（暖机 + 每步减半），而应基于上次状态判断；
**若已满足强耦合，就应改用强耦合（逐步同时耦合）继续算**。

**实现（5 源文件 + 3 文档；`make` EXIT=0、零告警）**：
- `sub_modules.f90`：`Iflag_Couple_Restart`（默认 0）+ 耦合状态全局量
  （`Couple_State_Mode/Conv/Pair/Iter/nf/nk/Twmax/Tol/Found`、`Couple_Tw/qw/u/pw_save`）。
- `sub_restart.f90`：`iver` 1→2，块数据后追加 coupling-state trailer
  （`mode/conv/pair/iter/nf,nk` + `max|dT_w|/tol` + 界面量 `T_w,q_w,u,p_w`）；
  读侧 `iostat` 容错（缺 trailer ⇒ 状态未知 ⇒ 旧行为）、0 号进程广播；
  新增 `set_couple_state(mode,conv,pair,iter,twmax,tol)`（从全局 `fp_*` 取数组、按形状匹配）
  与 `set_couple_tight_state()`。
- `sub_multi_region.f90`：三个交错驱动（11/13、12、19）入口处若 trailer 带同 pair/同形状的界面量
  ⇒ 直接恢复 `fp_Tw/fp_qw/fp_u/fp_pw`（**跨重启零跳变**，替代 `Twall_Couple_Init` 重播）；
  每轮算出判据后 `set_couple_state(1, merge(...), pair, it, metric, Tol)` 登记。
- `opencfd_ec3d_v1.16a.f90`：分派处按状态决定——`0` 自动 / `1` 强制强耦合 / `2` 强制交错 / `-1` 忽略；
  切强耦合时 `Iflag_Couple_Scheme=0`、`Solid_Max_Iter=min(现值,20)`、`Solid_Tol=max(现值,1e-6)`
  （**只动扫掠上限与容差**，松弛因子/物理不变），打印实际值；主循环退出前补写一次重启文件。
- `sub_read_parameter.f90`：6 步流程（默认值 / `$couple_ec` / 旧组 `$control_ec` / `bcast_para`
  `Ipara(58)` / `output_para.out` 回显）。
- 文档：本文件、`control.ec.template`、`docs/control.ec-说明.md §8`、`docs/重启与节点流场输出说明.md §8`。

**验证（`mpirun -np 1`，全部 EXIT=0、NaN=0）**：
- ct1 run1（case1 原设置）：12 轮，`max|dT_w|` 0.099545→0.0995036（与改动前逐轮一致），写重启 Kstep=2995。
- ct1 run2：`coupling state mode=1 satisfied=0 pair=11 66x1` → **界面 T_w/q_w 已恢复（零跳变）**
  → 继续交错 12 轮、EXIT=0。
- ct2 run1（`Comp=100/Warm=1/Outer=4/Min=10/Tol=1000`）：it=2 判为满足（0.0986 K）、写重启 Kstep=150。
- ct2 run2（`t_end=170`）：`satisfied=1` → **切逐步强耦合**（`Solid_Max_Iter=20/Tol=1e-6`）→
  150→155→…→170 逐步推进 → 退出前补写重启（Kstep=170）。
- ct3a（状态=强耦合，`t_end=190`）：识别 `mode=0` → 保持强耦合续跑 170→190；再跑一次 190→210 ✔。
- ct3b/ct3c（`Iflag_Couple_Restart=2 / -1`）：均走交错驱动（旧行为）；
  ct3d（`=1`）：强制切强耦合 ✔。
- 等价性：`output_para.out` 与改动前相比**仅多 1 行**（`couple restart state: Iflag_Couple_Restart= …`）。

**重要发现（写进文档）**：逐步强耦合下 `solid_solver_one_block` 每步都要把固体解到
`Solid_Tol=1e-9`（实测 1400 次扫掠）⇒ case1 ≈**15–20 s/步**；故切强耦合时自动 cap 预算
（20 次 / 1e-6）。另外 case1 `max|dT_w|` 平台 0.0995 K > `Tol_Couple_Tw=1e-2`，自动判定永不触发，
需放宽 `Tol_Couple_Tw` 或显式 `Iflag_Couple_Restart=1`。

**遗留**：本次改动**尚未提交**；集群需同步 `src/` 并重新 `make`。


---

## 2026-09-25（第三轮）— 建立 LaTeX 程序说明手册 v1.0（docs/程序说明/）

**需求（用户）**：用 LaTeX（子系统已装 texlive）编写程序说明文档，详细介绍固体、高速、低速、
多孔介质及其它新增功能，以及**如何控制与选用控制参数**；并以 fluid_solid 三个耦合算例
为例介绍算例设置；**后续新功能要持续更新该文档**。

**交付**：`docs/程序说明/`（`main.tex` + `chap01..11` + `91/92/93` 三个附录 + `build.sh` + `README.md`），
编译产物 PDF `docs/OpenCFD-EC-1.16a-程序说明.pdf`（**40 页**）。
工程要点：`\documentclass[fontset=fandol]{ctexart}` + **xelatex 两遍**；日志 `build.log`；
`build.sh [clean]`；`.gitignore` 已忽略 `*.aux/.log/.toc/.out`。

**章节**：1 概述（块类型/接口码/物理码/面编号/单位与无量纲化/目录/上手）；
2 输入文件与算例结构（`control.ec` 7 组规则 + 网格/bc3d/bc3d\_interface/material/porous/solid\_bc 格式与示例）；
3 可压高速（时间推进/空间格式/湍流/来流参考/输出节奏 + 选参建议）；
4 低速（变量约定/SIMPLE·SIMPLEC·AC-FV/`LS_*`·`AC_*` 参数表/**AC 标定经验表**/质量入口限制与等效速度换算）；
5 多孔（DBF+Ergun+LTNE 方程、`porous.inp`、预算 `AC_Max_Iter×Porous_Chunk_Iter`、参数选取、已知缺口）；
6 固体（FVM 非正交修正、SOR 参数、`Tw/Qw` 边界、**CHT 界面温度公式**、两种 CHT 模式、成本）；
7 跨区域耦合（同步 vs 交错、调度表公式与判据、性能对比、L 表机制、发汗壁、FAQ）；
8 输出与重启（文件总览、`field_restart.dat` 布局 `iver=2`、**重启耦合状态 `Iflag_Couple_Restart`**、
`flow3d_node.dat`、MPI 约定）；9 **三个 fluid_solid 算例详解**（拓扑/L 描述符、共同设置、
case1/2/3 逐项参数与预期日志、时间预算、续算实测 2995→5990、调参清单）；
10 构建/运行/并行与集群注意事项；11 验证/回归/已知问题；
附录 A **参数总表**、附录 B 编号与格式速查、附录 C **更新记录**。

**验证**：`./build.sh` 两次 xelatex 均 EXIT=0（无 error、无 undefined reference）；
`pdfinfo` 40 页 A4；`pdftotext` 抽查中文与目录正常。
编译期修掉的问题：`\code{}` 内裸 `_`(6 处)、`^`(2 处)、数学命令(text mode)1 处；
一次 `insert_line` 把 A.1 表的剩余行截断导致 longtable 不闭合（"input stack size exceeded"）——
现已修复（附录表结构完整）。

**维护约定（已写入 `.clinerules` 与 `constraints.md §0`）**：新功能必须 ① 更新对应章节
② 同步附录 A（新参数）③ 在附录 C 追加记录 ④ 重跑 `./build.sh` 确认无错误。

**遗留**：本次改动**尚未提交**。


---

## 2026-09-25（收尾）— 提交并推送：代码/功能 + LaTeX 程序说明

**动作**：把工作树里的改动整理为**两个提交**并 `push origin main`（本地与远端一致）：

1. `fix(io)+feat(couple)`：`control.ec` 单连接修复（`scan_control_ec_groups` 不再二次 open）
   $+$ 重启耦合状态（`iver=2` trailer、`Iflag_Couple_Restart`、交错续算界面量零跳变、
   自动切逐步强耦合并 cap 固体预算、退出补写重启）$+$ `control.ec.template` 与两份 docs 同步
   （8 files, +420/−13）。
2. `docs(manual)+chore`：`docs/程序说明/`（LaTeX 手册 17 个源文件 $+$ 交付 PDF 40 页）
   $+$ `.clinerules`/`constraints.md §0` 的“新功能必须更新手册”规则 $+$ `.gitignore`
   （忽略 LaTeX 中间文件与 `main.pdf`）$+$ `memory-bank/` 记录（25 files, +2212/−6）。

**收尾清理**：删除早前 `build.sh` 路径写错时落在**仓库根目录**的重复 PDF（未跟踪文件）。

**未提交（有意保留）**：`.vscode/settings.json`（主题设置）；`cases/**` 下的运行产物（本就 untracked）。

**提交前复核**：`make` EXIT=0 零告警；case1 12 轮等价（`output_para.out` 与改前逐字节一致）；
`ct1/ct2/ct3a-d` 全部 EXIT=0 / NaN=0；续算 Kstep 2995$\to$5990；
`docs/程序说明 && ./build.sh` 0 error / 0 undefined reference。

**记住（下次接手）**：提交号不进记忆库（避免陈旧）；新功能按 `.clinerules` 更新手册
（章节 $+$ 附录 A $+$ 附录 C $+$ 重新 `build.sh`）。



---

## 2026-09-25（第四轮）— flow3d_node.dat 补 PLOT3D 解(Q)/函数型文件头（块数 + 维数 + 变量数）

**问题（用户提出）**：`flow3d_node.dat` 按 PLOT3D 解(Q)/**函数型文件**读法加载时报格式错误 ——
文件里**只有裸数据**（`write(99) ((((U(m1,i,j,k),i=1,nx),j=1,ny),k=1,nz),m1=1,6)`），
既没有**块维数**也没有**变量数**，任何标准读端（Tecplot/ParaView）都无法识别，
只能靠外部 `Mesh3d.x` 去猜。

**方案（用户选定 A）**：加 PLOT3D 解/函数文件标准头 3 条记录，数据段保持原样。

**改动（2 个文件）**：
- `src/sub_IO.f90::output_flow_node`：
  ① 新 `integer,parameter:: NVAR_NODE=6`；
  ② 打开文件后（仅 0 号进程）依次 `write(99) Total_block` /
     `write(99) (bNi(m),bNj(m),bNk(m), m=1,Total_block)` / `write(99) NVAR_NODE`；
  ③ 上方注释块写明布局、变量含义与「旧文件无头」的兼容说明。
  数据记录（变量序 `d,u,v,w,T,Ts`、i 最快→与 PLOT3D Q 顺序一致）**一字未改**。
- `util/check_flow3d_node.py`：先按记录读入，自动判定「有头（记录1 为 4 字节整数 = 块数）/
  无头」；有头时校验 dims 与 `Mesh3d.x` 一致、`nvars=6`，随后按 `nvars*ni*nj*nk` 校验每块长度；
  两种布局都打印 `d/|u|/T/Ts` 范围与 `max|T-Ts|`。

**验证（`mpirun -np 1`）**：
- `cases/solid_1d`（1 块 51×9×9，`Iflag_flow_node=1, Kstep_save=1, t_end=0.3`）：EXIT=0；
  `xxd` 确认头 = `01 00 00 00`(nblock) / `33,09,09`(dims) / `06`(nvars)，数据记录 198288 B；
  checker：`format=PLOT3D-Q header (nblock=1, dims, nvars=6)` + `STRUCTURE OK`（T=Ts=294–353 K）。
- `cases/channel_ac_4blk`（4 块）：EXIT=0；7 条记录 = 3 头 + 4 数据，维数记录 48 B（4×3 ints）；
  checker `STRUCTURE OK`。
- 旧格式兼容：把新文件的 3 条头记录去掉后 checker 判为 `legacy (no header)` 且 `STRUCTURE OK`。
- **等价性 A/B**：用 `git show HEAD:src/sub_IO.f90` 在 `/tmp/old_src` 单独编译**改动前二进制**，
  同一算例两跑：`cmp` 结果「旧二进制输出 == 新文件去掉 3 条头记录」**逐字节一致**
  （solid_1d 198296 B、channel_ac_4blk 330656 B）⇒ 数据段零变化。

**文档/手册**：`docs/重启与节点流场输出说明.md §3`（新布局 + 替换掉旧的 `5*ni*nj*nk` 示例为可用
python 片段 + 指向 checker）；`docs/程序能力与算例考核总结.md §2.8`；手册 `chap08`（产物表 +
节点小节新增 lstlisting 布局与兼容说明）、`93_changelog.tex` 追加一行；
`cd docs/程序说明 && ./build.sh` → **40 页、0 error、0 undefined**。

**注意（手手相传）**：集群副本要同步 `src/sub_IO.f90` 并 `make`；旧 `flow3d_node.dat` 不必删
（checker 兼容），但要用 PLOT3D 读法直接加载必须**重跑一次**生成带头文件。
本轮改动**尚未提交**。

---

## 2026-09-26 — 交错外层停机判据 `Iflag_Couple_Stop`（t_end 或全局残差达标）

**需求（用户提出）**：「以 `t_end` 为结束标志，或者以全局所有的残差都小于残差标准」。
交错耦合（`Iflag_Couple_Scheme=1`）原先只有两个出口：界面量判据（`Tol_Couple_Tw`，码 12 另加
`p/u`）或跑满 `Niter_Couple_Outer`，`t_end` 完全无效（主程序调用驱动后直接 `stop`）。

**方案（用户选定 A：新开关，默认 0 保留旧行为）**
1. `src/sub_modules.f90`：`Iflag_Couple_Stop`（默认 0）+ 每轮残差登记量
   `Couple_Res_AC_q / Couple_Res_AC_m / Couple_Res_Solid`（轮首由驱动清零，求解器取 `max`）。
2. `src/sub_multi_region.f90`
   - 新 `couple_stop_test(nMesh, tt_reached, res_ok, r_comp, r_ls, r_por, r_sol)`：
     ① 统计 **BLOCK_FLUID** 胞心 `B%Res` 的 RMS（`m=1..NVAR` 取 max，`MPI_ALLREDUCE` SUM + 胞数）
     ⇒ 与 `Res_Inner_Limit` 比；② 存在 BLOCK_LOWSPEED/POROUS 时用 AC 登记量 ⇒ 与 `AC_Tol` 比
     （`LS_Algorithm≠3` 直接判为不收敛，不误判）；③ 存在 BLOCK_SOLID 时用 GS 登记量 ⇒ 与 `Solid_Tol` 比；
     全类块都达标才 `res_ok=T`（`-1.` 表示该类块不存在，仅用于日志）；
   - 三个交错驱动（`run_staggered_fluid_solid` 11/13、`run_staggered_highlow` 12、
     `run_staggered_fluid_porous` 19）：轮首清零登记量；气体段按剩余时间截断
     `nstep=min(nstep,⌊(t_end-tt)/dt_global⌋+1)`；轮末（`set_couple_state` 之后、原判据之前）调用停机测试，
     满足 `tt_reached .or. res_ok` 即 `exit outer`；退出提示区分「到 `t_end`」/「残差达标」
     （新增 `stop_time` 局部量，避免把“到点”误报成 not fully converged）；
   - `solid_solver_one_block` 末尾登记 GS `res`。
3. `src/sub_lowspeed_ac.f90` / `src/sub_porous_ac.f90`：AC 收敛后登记 `rq/rm`。
4. `src/sub_read_parameter.f90`：6 步流程齐全（默认值注释 / `$couple_ec` / 旧 `$control_ec` /
   `Ipara(59)` + 反解 / `output_para.out` 回显 print 残差标准）；`control.ec.template` 同步注释。
5. `util/readflow3d-ver2.4a/2.5.f90`：旧组 namelist 补 `Iflag_Couple_Stop`（顺手补上历史遗漏的
   `Iflag_Couple_Restart`），两者 `gfortran -c` 编译通过。

**验证（`mpirun -np 1`，全部 EXIT=0；临时目录 `/tmp/stgA*/B/C/E/F`，日志已回存算例目录）**
| 用例 | 配置 | 结果 |
|---|---|---|
| 码 19 `t_end` 停机 (`run_stoptest_tend.log`) | `t_end=25, Comp=10, Warm=99` | 第 3 轮气体段自动截断 **10→6**，`tt=26 ≥ 25` ⇒ `stopped: physical time reached t_end`；`ended at tt=26 ... outer iterations=3` |
| 码 19 残差停机 (`run_stoptest_res.log`) | `t_end=1e6`、`AC_Tol=Res_Inner_Limit=1e30` | 第 1 轮 `all_res_ok= T` ⇒ `stopped: all global residuals below standard` |
| 码 11 实测 (`run_stoptest_case1.log`) | `Solid_Max_Iter=2000` | `solid=7.09e-10`（与 GS 打印一致）、`comp=4.14` ⇒ `all_res_ok=F` ⇒ 跑满轮数（新判据不干扰旧行为） |
| 旧行为 (`run_stoptest_legacy.log`) | `Iflag_Couple_Stop=0` | 日志**无** `stop test` 行，与旧版一致 |
| 旧单组 `$control_ec` (`run_stoptest_legacy_group.log`) | 用 `control.ec.legacy` + 新增变量 | `legacy= T`、`Iflag_Couple_Stop=1` 生效、第 1 轮停机，EXIT=0 |
| 编译 | `cd src && make` | EXIT=0 |

**文档**：手册 `chap07` 新增 **§7.2 交错外层的停机判据**（两条件 + 残差口径表表 7.x + 示例 +
“SIMPLE 路径恒不达标”说明）、`chap07` 模式表与参数表各加一行、`91_appendix_params.tex`（组计数
12→13 + `Iflag_Couple_Stop` 行）、`93_changelog.tex`（2026-09-26 / 1.2）；`./build.sh` 重编译成功
（`main.toc` 出现 `subsection.7.2`，0 error、0 undefined）；`docs/control.ec-说明.md` 新增 §9、
组表改 13；`docs/程序能力与算例考核总结.md` §2.5 与 §4 组表同步。

**注意（手手相传）**：① 交错耦合仍**必须单进程**；② 该模式下 `t_end` 只在
`Iflag_Couple_Stop=1` 时生效，默认 0 与旧版逐字节同行为；③ 可压块的“残差标准”沿用
`Res_Inner_Limit`（默认 1e-10，双时间步内迭代判据），若要真正实现“稳态收敛即停”，
需按需放宽它或后续引入专用外层残差标准；④ 本轮与第四轮（`flow3d_node.dat` 文件头）
改动**均未提交**。



---

## 2026-09-26（第六轮）— 区域残差监视输出：固体/低速/多孔输出精简 + max/平均 rms/位置

**需求（用户提出）**：「固体，低速和多孔介质减少屏幕输出，输出区域最大 rms，平均 rms 以及
最大 rms 的位置」。用户选定**方案甲**：**不新增 namelist 参数**，硬编码
「每 100 次迭代/扫掠 1 行过程监视 + 1 行结束汇总」。默认输出即精简（删掉全部 DBG/DEBUG 行）。

**设计**：把"残差"从单个 `max` 扩成三个量，口径与可压块既有的 `Mesh%Res_max/Res_rms`
（`sub_time_advance.f90::comput_max_Res_onemesh`，`rms=sqrt(ΣR²/Num_Cell)`）保持一致：
`max` = 区域内逐胞残差幅值最大值；`rms` = 区域平均（均方根）；`at (i,j,k)` = `max` 所在胞
+ 胞心物理坐标 `[m]`（`%xc/%yc/%zc × Lscale`）。实现放在 `src/sub_debug.f90`（诊断例程文件、
已在 makefile 中 ⇒ 不改 makefile）里的新模块 `region_stat_mod`：

```fortran
module region_stat_mod          ! 常量 RES_PRINT_STEP = 100（硬编码显示常量）
  subroutine region_res_stat(nMesh, mBlock, what, kk, res_max, imax, jmax, kmax, res_rms, a1, a2)
```
- 只由 0 号进程打印 1 行：
  ` Region res [solid GS |dTs|] block 1 (BLOCK_SOLID) step= 1400 : max= … rms= … at (i,j,k)=(…) xyz=[…] m`
- `kk>0` → 过程监视行 `step=`；`kk<0` → 结束汇总行 `final=`（`|kk|` = 总迭代/扫掠数）。
- `what` 选残差名与行尾旧字段：1 固体 GS `|dTs|`；2/3 低速/多孔 AC 流场（附 `res_q/res_m`）；
  4/5 低速/多孔 SIMPLE `div(F)`（附 `res_p`）；6 AC 能量 `|dT|`；7/8 LTNE 流体/骨架（附 `dT(window)`）。
- 字符字段全部在例程内生成（调用方只传标量），因此**隐式接口也安全**。

**改动（6 个源文件）**
1. `src/sub_debug.f90`：如上新增模块（+ 说明注释）。
2. `src/sub_multi_region.f90::solid_solver_one_block`：GS 内层循环同时累加 `rss_r += |ΔTs|²`
   并跟踪 argmax（原来只有 `res=max(...)`）；每 100 扫掠 1 行、结束 1 行；**删除**
   `GS iter`（前 5 次）与 `DBG: Ts(...)` 两处调试输出；`Couple_Res_Solid` 登记逻辑不变。
3. `src/sub_lowspeed.f90` / `src/sub_porous.f90`：`lowspeed_/porous_pressure_correction` 追加
   出参 `imax,jmax,kmax,res_p_rms`（在既有 `res_p` 的 `|div F|` 循环里顺带累加，无新增数组）；
   驱动每 100 次外迭代 1 行 + 结束 1 行；**删除** `DBG iter/maxU2/div max`、每 50 步 `DEBUG T iter`、
   `SIMPLE iter`（含前 5 次）与末尾 `SIMPLE iterations: res_p` 行；NaN 检测与 stop 保留。
4. `src/sub_lowspeed_ac.f90`：新外部例程 `ac_res_stat`（与 `ac_res_norm` 同一套逐胞量
   `vq/vm/1e-3*rr`，额外返回合并幅值的 max/位置/rms 与 `rq/rm`）；AC 结束打 1 行汇总，
   删掉**重复的** `iter==1` 打印（合并后的 `mod(iter,AC_Print)` 行已覆盖 iter=1）；
   `ac_lowspeed_energy`（AC 温度方程）在 50 步窗口统计 `|ΔT|` 的 max/rms/位置，删掉每 1000 步进度行，
   结束打 1 行。
5. `src/sub_porous_ac.f90`：AC 结束打 1 行汇总、删 `iter==1` 行；`ac_porous_ltne` 在 100 步窗口
   同时统计流体/骨架两组量（argmax + Σr²），删掉每 2000 步进度行，结束各打 1 行。
6. 各调用例程加 `use region_stat_mod`（**注意**：漏加会在链接期报 `undefined reference`，
   本轮在 `ac_lowspeed_energy` 上踩过一次）。

**验证（`mpirun -np 1`，临时目录 `/tmp/rt1/*`，全部 EXIT=0）**
| 用例 | 覆盖路径 | 结果 |
|---|---|---|
| `cases/solid_1d`（t_end=0.3d0） | 固体 GS | 每 100 扫掠 1 行 + `final= 20001`；`max` 与旧 `final res` 一致（3.9e-7） |
| `cases/channel_ac_4blk` | AC 低速 4 块 | 4 条 `low-speed AC flow`，行尾 `res_q/res_m` 与 `AC iter` 自洽 |
| `cases/porous_ltne_1d`（Porous_Max_Iter↓） | 多孔 SIMPLE | `porous SIMPLE div(F)` 过程行 + `final=`（res_p 单调降至 1.8e-11） |
| `cases/porous_ltne_1d/run_ac` | AC 多孔 + LTNE | `porous AC flow` + `LTNE fluid/frame` 三行，位置在末端/入口附近 |
| `cases/porous_plug_ac/eps06`（AC_Max_Iter=3000） | AC 多孔（hv=0） | 仅有 `porous AC flow` 汇总行（无 LTNE，符合预期） |
| `fluid_solid/run_case1_solid`（码 11，2 轮，Comp=5） | 耦合内固体 | 轮内固体汇总行正常；`solid chunk done`/退出信息不变 |
| `fluid_solid/run_case2_lowspeed`（码 12，1 轮） | 耦合内 AC + AC 能量 | 低速汇总行 + `AC energy` 行 + `highlow interface metrics` 正常 |
| `fluid_solid/run_case3_porous`（码 19，1 轮） | 耦合内多孔 AC + LTNE | `porous AC flow` + LTNE 两行；`porous chunk done … max|dT_w|=4.38 K` |
| `cd src && make` | 编译 | EXIT=0（无告警） |

**文档/手册**：`chap06`（固体日志小节改为统一格式，含 `lstlisting` 样例 + 三量说明）、
`chap04`/`chap05` 各新增「残差监视输出（2026-09-26 起）」小节、`92_appendix_codes.tex`
日志关键字表新增 `Region res ...` 行（并替换已删除的 `Solid solver Block ... GS iterations:` 行、
补 `The Max/R.M.S Residuals are` 行）；`93_changelog.tex` 新增 **1.3** 行；
`cd docs/程序说明 && ./build.sh` → **45 页、0 error**；`docs/程序能力与算例考核总结.md` 新增 §2.9。
（首次编译报 `Missing $ inserted`：`\code{AC iter ... res_q= ...}` 内下划线未转义 ⇒ 已改为 `res\_q`。）

**注意（手手相传）**：① 打印节奏只由 `region_stat_mod::RES_PRINT_STEP=100` 决定（硬编码，
非 namelist；要更少输出改这一处，或改 decade 式）；② 长稳态算例会出现每 100 步 1 行的趋势表
（`solid_1d` 20000 扫掠 ⇒ ~200 行）；③ **收敛判据与监视节奏重合时不重复打印**：
固体/SIMPLE 驱动用局部 `conv_ok`，即将退出时跳过该 `step=` 行（实测 case1 第 2 轮
只在 `final= 1400` 打一行，不再出现 `step= 1400` + `final= 1400` 两行）；
④ 耦合场景不变：交错仍必须单进程；
⑤ 本轮与第四、五轮改动**均未提交**。


---

## 2026-10-17 — 第四～七轮改动分组提交 + 多孔 ε 极限“负结果”澄清 + 文档补全

**目标**：把第四～六轮（节点头、停机判据、区域残差）、第七轮（ε 极限 + limits 算例）的
未提交内容分组提交，并**如实**记录哪些极限性质已验证、哪些没有。

**提交**：`c974e81` feat(output)（3 条 PLOT3D 头 + 冷启动确定性 + 两个 util 脚本）、
`ebe7053` feat(solver)（`region_stat_mod` + `Iflag_Couple_Stop` + ε 极限 + 模板/读参/util）、
`a29a6ee` cases(limits)（5 算例 + `run_limits.sh` + README）。

**关键事实**：① HEAD 的 `output_flow_node` 原本**无头**，工作树曾写**2 条**过渡布局
（块数 + 每块 `ni,nj,nk,nvar`），现为**标准 3 条**；`check_flow3d_node.py` 自动判别两种布局
（打印 `2-record variant` / `3-record PLOT3D-Q header`），2026-09-26 的 limits 文件仍可读。
② `Iflag_restart<0` 时 `restarted` 未初始化（偶发从未初始化场起步/NaN），已先置 `.false.`。
③ limits 实测（`--blocks 1`）：ε=1 vs 低速 max|Δu|=17.64、max|Δv|=28.02、max|ΔT|=96.34 K，
**排除界面首层后指标不变**；ε=0 vs 固体 max|ΔT|=498.6 K；多孔块 T 停在 300 K、Ts 最多 +3 K
⇒ ε 极限目前只是“分支不崩 + ε=0 速度逐点为 0”的**回归快照**，不是定量等价证据。

**文档**：`chap05`/`chap09` 改写口径；`chap11` 新增 §11.4 实测表（负结果 + 原因 + 后续）、
§11.5「2026-10-17 验证记录」；`93_changelog.tex` 追加 1.4 / 1.4.1；`docs/程序能力与算例考核总结.md`
§2.4/§5/§6-15；`cases/fluid_solid/limits/README.md`（实测表、分析、后续）。`./build.sh` 两次通过
（0 error、第二遍 0 undefined）。

**验证**：`cd src && make -j4` EXIT=0（无 warning）；`solid_1d`（1 块）/`channel_ac_4blk`（4 块）
`STRUCTURE OK`；`bash -n run_limits.sh`、`python3 -m py_compile` 两个脚本通过。

**注意**：`workplan.md`（仓库根，**GBK 编码**的用户笔记）与 `.vscode/settings.json` 有意未提交；
6 个新 `Mesh3d.x` 单独提交。集群副本（`/work/home/lijunyang/sundong/PorousTest/code/`）
仍需同步 + `make`，旧 libgfortran 的双 unit 检查本地测不出来。

## 🗓 2026-10-17 · 程序更名 MultiPhysMixSolver（数值惰性）

**改动**：`git mv src/opencfd_ec3d_v1.16a.f90 src/multiphysmixsolver.f90`；`src/makefile`
（`EXE = MultiPhysMixSolver`，链接后 `ln -sf $(EXE) opencfd-ec1.16a.out`，`clean` 两者都删）；
启动横幅 5 行 + 文件头注释；`sub_read_parameter.f90:383` 的 `output_para.out` 版本行；
22 个 `cases/**` 脚本/README；`.vscode/{tasks,launch}.json`；手册 `main/chap01/chap09/chap10/93_changelog`
（v1.5）与 `build.sh`；`docs/程序能力与算例考核总结.md`、`docs/重启与节点流场输出说明.md`。

**有意不改**：`flow3d.dat` / `Ts_block_*.dat` / `bc3d.inc` 的数据头字符串仍是 `OpenCFD-EC`
（后处理兼容）；`field_restart.dat` 为纯数值记录，与文件名无关；`docs/工作日志.md`、
`docs/进度_fluid_solid三算例.md`、`memory-bank/{worklog 旧条目,progress 旧条目,archive_*}` 保留历史原名。

**验证**：串行 `make clean && make` EXIT=0；以 `git show HEAD:` 复原更名前源码编译旧二进制，
与新版在 `channel_ac_2blk`（2 块 + AC，5 步）同输入对跑：`flow3d.dat` / `field_restart.dat` /
`Residual.dat` / `force.log` / 两块 `*.vtk` / `Step_mess.dat` **逐字节相同**，
`output_para.out` 仅第 1 行版本串不同（预期）；`docs/程序说明/build.sh` 两遍 0 error（49 页）。

**坑**：`make -j4` 在本 makefile（无模块依赖声明）下会竞态报
`Cannot open module file 'precision_ec.mod'` ⇒ **一律串行 `make`**。

**下一步**：并行改造 P0（跨进程安全网）→ P1（`p` halo）→ … → P6，见 activeContext §0h。

---

## 🗓 2026-10-17 · P0 跨进程安全网（静默不耦合 → 立即 MPI_Abort）

**问题**：接口 12（可压↔低速）、14（固↔固）、17（固↔多孔）、18（多孔↔多孔）、
19（可压↔多孔）只有**同进程**耦合例程。当 `partation.dat` 把两侧分到不同进程时，
`couple_*` 例程直接 `cycle` 跳过 ⇒ 作业“成功结束”但两侧物理上没耦合（最危险失败模式）。

**实现**：`src/sub_multi_region.f90` 末尾新增
- `check_remote_coupling_support(nMesh)`：遍历 `bc_msg2`，只统计
  `B_proc(nb) /= my_id` 且组合无 MPI 版本的条目（12/14/17/18/19），
  已支持的 11/13（`cht_mpi_fluid_solid_face`）、15（U 缓冲）、16（`ls_por_iface_mpi`）
  和 Ts 缓冲交换不计入；`Total_proc<=1` 立即返回。
- `abort_remote_uncoupled(n_local, what)`：`MPI_Allreduce(MPI_SUM)` → 0 号进程打印
  清单 + 修复建议（`-np 1` 或把两侧放同一进程）→ `MPI_Barrier` + `MPI_Abort(2)`。
调用点：`src/multiphysmixsolver.f90` 的 `Init_flow` 之后、`Start ......` 之前
（一次遍历、毫秒级；放在耦合例程里要等第一步算完才报错，故未采用）。

**验证**
- `make` EXIT=0（串行）。
- np=1 `channel_ac_2blk`（2 块 AC，5 步，`Kstep_save=5`）与 P0 前二进制：
  `flow3d.dat` / `field_restart.dat` **逐字节相同** ⇒ np=1 零影响。
- np=2 `cases/fluid_solid` 副本把 `material.in` 改 `0 2`（可压↔低速，界面跨进程）：
  **EXIT=2**，日志
  `ERROR: cross-process coupling is NOT implemented for 12(compressible<->low-speed)`
  `remote interface entries (summed over all ranks) = 2`。
- np=2 原 `cases/fluid_solid`（多孔 3 ↔ 低速 2，界面 16）越过自检、进入
  `AC porous solver` 迭代 ⇒ 无误报。
- 手册 `chap07` 新增《跨进程安全网（启动自检）》小节（含报错样例），
  对照表“强制要求”行补注；`93_changelog.tex` 追加 1.5.1；`./build.sh` 重编译。

**下一步**：P1 —— `p` 跨进程 halo（`p_send/_recv_mpi` 镜像 `Ts_*`）、
`ac_have_pressure_bc` 的 `MPI_LOR`、`ac_zero_mean_p` 的全局均值。

---

## 2026-10-17 — 构建布局重构：`makefile + src/ + build/ + bin/ + cases/`

**需求**（用户）：“将编译环境重新设置为：makefile+src+build+bin+cases，目前中间文件混在代码里面太乱了。”

**改动**
- 新增顶层 `makefile`：`make` / `make -j4` / `make clean` / `make info`；模式规则
  `$(BUILDDIR)/%.o: $(SRCDIR)/%.f90`，`-Jbuild -Ibuild`；链接出 `bin/MultiPhysMixSolver`
  并 `ln -sf` 建 `bin/opencfd-ec1.16a.out`（旧名兼容）。
- `src/makefile` → 转发入口（`cd src && make` 保持可用）；`src/` 只留源码。
- 新增 `util/gen_moddep.sh`：awk 扫 `module`/`use` 生成 **`src/modules.dep`**（36 条文件级依赖，入库）；
  顶层 makefile `include` 它 ⇒ 从干净树 `make -j4` 也安全。
  **不用** gfortran `-MMD`（要求 `-cpp`，会改变编译语义；实测 `-cpp -MMD` 可编译但放弃）。
- `.gitignore`：`/build/`、`/bin/`、`*.d`；`.vscode/{tasks,launch}.json` 改到根目录 cwd + `bin/` 路径。
- 22 个 `cases/**` 的 README/脚本运行路径 → `../../bin/MultiPhysMixSolver`
  （`run_limits.sh` 用 `$R/bin/...`，`run_grid.sh` 用 `../../bin/...`）。
- 清理：`src/` 44 个中间文件（`*.o`/`*.mod`/编译日志/可执行/旧软链）、根目录 2 个 `*.mod`
  与 `snapf.log`、`cases/**` 20+ 个旧 `opencfd-ec1.16a.out` 二进制副本。

**验证**
- `make clean && make -j4` → **EXIT=0 / 89 s**（重构前从干净树的 `-j4` 必失败）。
- `channel_ac_2blk`（2 块 + AC，5 步）`flow3d.dat`、`field_restart.dat` 与重构前**逐字节相同**。
- P0 安全网负例（np=2，码 12 跨进程）仍 **EXIT=2**；np=1 正常 EXIT=0。
- 手册升 **v1.6**（`chap01`/`chap09`/`chap10`/`chap11`/`93_changelog` 1.6，`main.tex` 版本号），
  `./build.sh` **0 error**。

**坑**：`find cases -maxdepth 6 -name opencfd-ec1.16a.out -delete` 会删掉算例目录里的旧二进制副本
（本就是 `*.out` 忽略项，可随时从 `bin/` 重新拷贝）；如需保留请先备份。

**下一步**：P1+P2 合并实施（AC 的 `p` halo + 伪步统一步数，见上一条的下一步）。




## 2026-10-17 — AC 的 `p` halo 统一（方案 A）+ `IF_Debug=3` 缓冲区交换自检

**背景（用户审阅指出）**：`p` 存于 `B%p`（不在 `U(1:5)`），其同进程 ghost 原先由
`Umessage_send_mpi` 内一个“融合块”拷贝：

- 过滤条件窄（`LOWSPEED↔LOWSPEED`），比 `U` 的（`bc<0`，不看类型）和
  `p_send_mpi` 的（`LOWSPEED|POROUS` 两侧）都窄；
- 依附于 `U` 的调用点，而 `p_send_mpi` 只在 `AC_MPI_Sync==1` 分支里调用
  ⇒ 缺省 `AC_MPI_Sync=0` 时 **np=1 有新鲜 `p`、np>1 从不刷新**（`ac_load_state` 把
  `p` 写进 `XW(1)`，直接影响内部重构 → np 相关结果）；
- `AC_MPI_Sync=1` 时同进程 `p` 被拷**两遍**；
- 多孔驱动在 np=1 时**完全不**刷新 halo（`else if(Total_proc==1)` 分支缺失）。

**方案 A 落地**
1. `src/sub_update_buffer_mpi.f90`
   - 删除 `Umessage_send_mpi` 里的融合 `p` 拷贝；
   - `p_send_mpi` 重构：先判 `Send_to_ID`，同进程**直写**（无 `psend` 临时数组、
     无 MPI），跨进程才打包 `MPI_Bsend`，`deallocate` 后置；`associated(B%p)` 守卫。
2. `src/sub_lowspeed_ac.f90` / `src/sub_porous_ac.f90`
   - 伪步循环内统一为 `if(AC_MPI_Sync==1 .or. Total_proc==1)` → 同时
     `update_buffer_onemesh` + `update_p_buffer_onemesh`；多孔补齐 np=1 分支。
3. `src/multiphysmixsolver.f90` + 新例程 `check_buffer_exchange(nMesh)`
   - `IF_Debug=3` 时在 `Init_flow` 之后跑一次自检：把 `U(1:NVAR)`/`p`（含 ghost）
     填成 `((blk*1000+i+16)*1000+j+16)*1000+k+16 + m/8`（整数可精确表示），
     做一次 `U`+`p` 交换，逐 ghost 胞与“映射源胞哨兵值”**严格**比较；
   - 同进程子面直接比；跨进程子面因本 rank 无邻居块，改为把“邻居 ghost 应有的
     内容”按独立 tag（`CK_TAG=900000`、p 加 500000）发给对侧再比；
   - 先保存后恢复 `U`/`p`，**非破坏性**；`MPI_Allreduce` 汇总，末行 `PASS/FAIL`。

**验证**（`/tmp/acbuf/{A,B,C}`，`t_end=20`、`Kstep_save=20`、`AC_MPI_Sync=1`）

| 分布 | 算例 | 子面 同/跨 | ghost 胞 | 结果 |
|---|---|---|---|---|
| np=1，2 块 | `channel_ac_2blk` | 2 / 0 | 320 | `max|dU|=max|dp|=0` PASS |
| np=2，各 1 块 | `channel_ac_2blk` | 0 / 2 | 320 | 0 / 0 PASS |
| np=2，各 2 块 | `channel_ac_4blk` | 4 / 4 | 960 | 0 / 0 PASS |

三例均 EXIT=0 并正常写出 `field_restart.dat`。

**踩坑（写进 activeContext §0j）**
1. 期望值最初用 `nb1`（**接收**块号）→ 必然差 ~1e9（正确应用**发送**块 `B%Block_no`）；
2. 跨进程**不可** `B1 => Mesh%Block(B_n(nb1))`：`B_n` 只对本进程块有效，会指向错误块
   （假报警，看起来像“ghost 没被写”）；
3. `size()` 已含 NVAR 维 → `NVAR*size()` 让 `MPI_Recv` 声明 5 倍缓冲，`MPI` 越界写
   导致 malloc assert / SIGSEGV。

**文档**：`chap04` 新增 §`ac-halo`（调度语义表）与 §`bufchk`（自检说明 + 三分布结果）；
附录 A 补 `AC_MPI_Sync`（`$ac_ec` 15→16）与 `IF_Debug=3`；`93_changelog` 1.6.1；
`docs/程序说明/build.sh` 0 error。

**下一步**：`ac_zero_mean_p` 改全局均值；`AC_MPI_Sync=0 & np>1` 自动升级或明确报错；
`FSH` 注册表跨进程回退（两侧各算一阶 Rusanov）带来的 np 相关精度。

## 2026-10-17 — 手册 §4.7「多块离散一致性」+ 曲面（O 型）网格粘性格式适用性

**触发**：审阅提问 ① 多块 AC 的对流/粘性离散是否需要特殊处理、二阶 Laplacian 为什么
不需要角点数据；② "三个一维二阶差分"在不横平竖直的网格（圆柱 O 型）上会不会失效。

**代码核对（只读）**
- 对流 `ac_face_flux`：法向 5 点 `a..e`（`AC_Recon=1` 只用 `b..e`，1 层 ghost；WENO5
  用 `a`，2 层）；内面循环 `ii=2..nx-1` **排除界面面**（防双计）；界面面由
  `ac_ifc_flux_own`（owner 单值 `FSH`，全局块号定 owner，`ac_ifc_src_cell` 复用缓冲
  交换的 `L1/L2/L3` 映射）或跨进程回退 `ac_face_flux_rus1_own` 处理，**两者都只用一阶
  Rusanov**（界面 ghost 每外步才刷新，二阶重构会激发 checkerboard）。
- 粘性 `ac_viscous_res`（`sub_lowspeed_ac.f90:1830--1893`）：循环 `f=1..nx` **不排除**
  界面面，`if(f-1>=1)/if(f<=nx-1)` 只累加到真实内胞；`H=ν(u_f−u_{f-1})/d`，
  `d=0.5|x(f+1)−x(f-1)|`，斜网格取 `dproj=|Δx·n|`（`>0.2d` 才用）；两侧分子相同
  ⇒ `d_A=d_B` 时**恰好等值反号**，无需 `FSH`。壁面 `d=Δn=2d_wall` ⇒ `H=−νu_1/d_wall`。
- 只含 `∂xx+∂yy+∂zz`、**无混合导数** ⇒ 7 点面邻居壳层闭合，ghost-角点从不被读取；
  即使写成应力张量形式，最坏也只需 ghost-**edge**。

**曲面/非正交网格（回答 ②）**：不失效，前提三条 ——① 速度存**笛卡尔分量**（本代码
如此）⇒ 逐分量标量 Laplacian 与坐标系无关，**不需要**柱坐标曲率项；② ν 常数
（`LS_mu`）；③ 法向导数用"索引线差分 ÷ 法向投影"（过松弛法）：正交网格严格相等，
O 型网格 `i`=径向 / `j`=周向时 skew = `O(Δθ)=O(h)`，仍二阶。FV 用真实面积/体积、
无显式 `1/r` ⇒ **无坐标奇点**。实测（`cases/cylinder_re40_half`，O 型极坐标，
120×80 / 160×120，Re=40，AC）：`Cd` 1.5706→1.5536 vs Rogers&Kwak 1.549（+0.3%），
`Cd_p` −0.2%，壁面积分 ↔ 远场动量平衡互校（1.5706↔1.5708）。会降阶的情形：恒剪切
（非正交度不随 h 减小）→ 一阶；`μ(T)` 未实现；界面非共形/间距跳变 → 界面不守恒。

**文档**：`chap04` 新增 §4.7（表 17 依赖表 + 对流三件事 + 界面粘性等值反号表 +
角点判据 + 曲面网格适用性 + 表 18 实测）；`93_changelog` **1.6.2**；`main.tex` 版本
1.6.2；`./build.sh` **0 error / 0 undefined**（第一遍 `sec:ac-mblock` 前向引用在第二遍
消解，PDF 中无 `??`）。

**下一步**：`ac_zero_mean_p` 改全局均值；`AC_MPI_Sync=0 & np>1` 自动升级或报错；
`FSH` 跨进程回退；若将来上 9 点/紧致/顶点型格式，**必须先扩 ghost 角点交换**。


## 2026-10-17 3 块 AC 方腔：界面通量符号 + `is_high` 判据两处致命 bug 修复

**任务**：`cases/lidcavity_ac_3blk`（Re=100，42×42，三块等分 x，界面 x=1/3、2/3）在
np=1/2/3 下都应给出与同网格单块一致的正确解。

**现象**：界面胞（块 1 的 `i=nx-1`、块 3 的 `i=1`）出现棋盘模态并指数增长
（增长因子≈1.3/步）；压住发散后界面仍整列错位、u 达 1.68/1.83。

**定位**：错误解在 np=1/2/3 **逐位相同** ⇒ 竞争 bug 的可能性排除，改走 np=1 隔离验证
（临时 `DBGS` 打印），直接看到块 1 的 `iG=15>iO=14`（ghost 索引更大）被判成高侧。

**修复（`src/sub_lowspeed_ac.f90::ac_face_flux_rus1_own`）**
1. 两侧原本都 `RAC += FLX*A`，但两侧传入的面法向相同（i 面恒 +x），Rusanov 耗散项
   −λ/2·(q_R−q_L) 在界面上相互抵消 ⇒ 界面棋盘模态零阻尼。改为低侧 `RAC -= FLX*A`、
   高侧 `RAC += FLX*A`。
2. `is_high` 判据写反 ⇒ 中央通量两侧同向反向、界面变抽吸源。改为按 own/ghost 索引关系
   判定（`iG>iO`/`jG>jO`/`kG>kO` ⇒ 低侧），与块号、进程号无关（旋转/反向界面同样成立）。

**验证**
- np=1/2/3（10^4 步）：`max|Δu|=1.16e-6`、`max|Δv|=2.4e-7`、`max|Δp|=9.8e-7`；
  np=2 与 np=3 逐位相同。
- 3×10^4 步收敛：`res_q≈3e-14`，u∈[−0.233,0.926] ⊂ [−1,1]；vs Ghia RMS u 0.02095 /
  v 0.04365。
- 同网格单块 42×42 对照（新增 `gen_single42.py`）：`max|Δu|=9.3e-3`、rms 3.0e-3，
  正是界面**刻意一阶**的预期误差（h=1/42、∂u/∂x≈0.4）；离界面 3 胞以外逐位一致。
  单块 42×42 vs Ghia：u 0.02013 / v 0.04387（与已验证单块 40×40 的 0.02095 同量级）。
- `AC_Max_Iter=1` 是这类多块算例的**要求**（np=1 同伪步刷新与 np≥2 外步末刷新节奏一致）。

**文档/收尾**：`chap04` 新增「两侧各算路径的符号约定」段 + 表 `tab:ac3blk`；
`93_changelog` 1.6.3；`main.tex` 1.6.3；`./build.sh` 0 error；删除临时调试块，
`ac_closed_dom` 与 `ac_dbg` 复位；`make` 干净重建通过。

**下一步**：清理临时目录（`s1`/`s42`/`m42`/`v1`/`v3`/`q*_10k`）；可选：界面一阶误差的
h 减半收敛性表。


## 2026-10-17（收尾）3 块 AC 修复提交入库：4 个 commit + 提交前一致化与冒烟

**提交**
| commit | 内容 |
|---|---|
| `29eae92` | fix(ac)：`ac_face_flux_rus1_own` 符号 + `is_high` 判据；`AC_MPI_Sync`（Ipara(60)）；`update_p_buffer_onemesh`；`init_ac_parallel_flags`（Allreduce + 闭腔 1 号块锚定）；`IF_Debug=3` 缓冲自检；`read_solid_bc` 缺文件时先广播（修 np>1 启动死锁）；`control.ec.template` / `91_appendix_params` / `channel_ac_2blk` 同步 |
| `10ee6e2` | cases(lidcavity_ac_3blk)：`gen_case.py`/`gen_single42.py`/`analyze_3blk.py`/`prep.sh`、两份分区文件、`Mesh3d.x`/`bc3d.inp`/`material.in`/`control.ec`、README、目录内 `.gitignore`（`np*/`、`s42/`、`q*_10k/` 等实测目录不入库） |
| `6a804fa` | docs(manual)：1.6.3（符号约定段 + `tab:ac3blk` + changelog + 两份 PDF） |
| `5094fd8` | chore：删除根目录遗留 `control.ec` |

**提交前一致化（发现并修掉 3 处"文档/配置 ≠ 代码"）**
1. `src/sub_lowspeed_ac.f90::ac_res_norm` 里仍有 `IF_Debug=1` 的临时 `DBG2` 统计块
   （每步多一遍 O(N) 扫描 + 打印）→ 删除；`ac_dbg=1` 的 `IFC DIAG` 保留（默认 0）。
2. `cases/lidcavity_ac_3blk/control.ec` 与实测配置不符（`AC_Max_Iter=50`、`Kstep_save=250`、
   `t_end=1000`）→ 改为 `AC_Max_Iter=1`、`t_end=Kstep_save=10000`，并把"为什么必须 1"写进注释。
3. `control.ec.template` 还说"闭腔 + np>1 启动即报错"，但 `init_ac_parallel_flags` 已改成
   "全局 1 号块锚定压力水平"（np2 日志确认为提示而非 Abort）→ 模板注释同步。

**冒烟（提交前，新二进制）**
- 造 `/tmp/sm1`、`/tmp/sm2`（软链输入，`t_end=2000`、`Kstep_save=2000`、`AC_Max_Iter=1`、
  `AC_MPI_Sync=1`），np1 与 `mpirun -np 2` 各跑 2000 外步：均 EXIT=0、无 NaN、u∈[−0.207,1]；
  `np2 vs np1`：`RMS(u)=4.0e-4`、`max|Δu|=8.7e-4`（瞬态，与收敛后 1.16e-6 同趋势）。
- 注：`cases/lidcavity_ac_3blk/np1/run.log` 是 22:00 的**旧 bug 版**日志（尾部 NaN），
  固定结果在 `v1`/`v3`/`q*_10k`（23:12 之后）；这些目录都不入库。
- `make -j2` EXIT=0；`.tex` 未再改动，PDF（23:36）比所有 `.tex` 新。

**未提交（有意）**：`.vscode/settings.json`（主题）、`cases/ssh/集群信息.md`、`workplan.md`
（另需你确认是否入库 + `workplan.md` 当前是 GBK 编码）。

**下一步**：清临时目录 → 输出含 `p` 的 VTK 后做 h 减半的界面误差表 → 视需要把 3 块算例
接入 `util/` 的自检脚本。

## 2026-10-18 — 集群 `wuzh02` 连接打通（首次成功登录）

- 依据 `cases/ssh/集群信息.md`（GBK，已转 UTF-8）：`wuzh02.hpccube.com:65091`，用户 `lijunyang`。
  本机 DNS 只给 IPv6（`2409:8728:4eff:2001::100`），65091 `TCP_OK`。
- 私钥两份（`/mnt/c/Users/Administrator/.ssh/lijunyang_wuzh02_rsa` 与 Downloads 里的
  `..._RsaKeyExpireTime_2026-11-12_14-15-23.txt`）**md5 相同**；`/mnt/c` 权限 777 ⇒ ssh 拒用，
  故复制到 `~/.ssh/lijunyang_wuzh02_rsa` + `chmod 600`。
- 本机**没有** `~/.ssh/config`（`ssh wuzh02` 解析失败）⇒ 新建 `Host wuzh02`（Port 65091、
  IdentityFile、`IdentitiesOnly yes`、`StrictHostKeyChecking accept-new`）。
- 验证：`ssh -o BatchMode=yes wuzh02 'hostname; whoami; pwd; uname -sr'` →
  `login03` / `lijunyang` / `/work/home/lijunyang` / `Linux 3.10.0-957.nfs.5.x86_64`，EXIT=0。
- 环境：`mpif90` = **gfortran 7.3.1**（hpcx v2.11.0，本地 13.x）；Slurm/GridView
  （`sbatch/squeue/sinfo`，分区 `wzhcnormal`/`wzacnormal`）；`/work` 30 TB 用 93%（剩 2.3 TB）；
  `rsync/scp/tar/md5sum/nproc` 齐全。旧副本 `~/sundong/PorousTest/code/` 仍是单文件老布局。
- 未做：`git push`、源码 `rsync` 同步、提交 Slurm 作业。这些待用户确认后执行。
- 补充：集群没装 `C.UTF-8`，登录脚本 `perl` 每次告警 ⇒ `~/.ssh/config` 加
  `SetEnv LANG=C` / `SetEnv LC_ALL=C`，实测告警消失（`LANG=C`）。
- 补充：`rsync -az -e 'ssh -o BatchMode=yes' <file> wuzh02:./` 上传成功、远端 md5 与本地一致
  （探针文件已删）⇒ 整仓同步通道就绪。
- 作业提交（用户提供 + 集群实读）：帮助页
  <https://www.scnet.cn/help/docs/mainsite/hpc/cmd/slurm/>；模板 `~/sundong/run.slurm`
  （`-p wzacnormal -N 1 -n 1 -J MixSolver -o/-e *_out|err.log` +
  `mpirun -n 1 ../../bin/mixnsolver`）。用法：**拷到算例目录 → `sbatch run.slurm`**；
  跑本仓库需改 `mixnsolver`→`../../bin/MultiPhysMixSolver`、`-n` 改成 np、`export PATH` 可删。
- 核实：本地根目录**没有** `run.slurm`（`find` 全盘无、未被 `.gitignore` 忽略）；
  `~/sundong/MixSolver/MixSolver/`（`mixnsolver`/`mixnsolver_mpi`）是**另一个项目**。
  本仓库仍**没有集群副本**，首次同步要新建目录。



## 2026-10-18 — 集群部署（wuzh02）+ gfortran 7.3.1 编译 + Slurm 复跑 np=1/2/3

- 新规则：`.clinerules` 与 `memory-bank/constraints.md §8` 增加“本地 ↔ 集群同步”，
  工具 `util/sync_cluster.sh {push|pull}`（push 排除 `.git/build/bin`、VTK、重启/网格大文件、`np*/`；
  pull 只取 `memory-bank/`、`docs/`、`cases/**/*.log`、`analyze_*`）。
- 部署：push **1943 文件 / 50 MB** → `/work/home/lijunyang/sundong/MultiPhysMixSolver/`；
  集群 `make clean && make -j4` = **0 error / 2 warning / 68 s**（38 个目标），
  `bin/MultiPhysMixSolver` 3.15 MB；源码 md5 与本地一致（`sub_lowspeed_ac.f90` 0ad488b7…）。
- 复跑：`cases/lidcavity_ac_3blk` 1e4 步，Slurm 46812863/64/65 全部完成（墙钟 27/18/8 s）；
  np2、np3 对 np1 = RMS(u) 5.77e-07 / max|Δu| 1.16e-06（两者数值相同）；
  **集群 np1 与本地 `q1_10k` 逐位相同（RMS=0）**。
- 工具链/脚本坑：`prep.sh` 旧写法 `> np*/control.ec` 截断了软链目标（仓库 `control.ec` 一度 0 字节，
  已 `git checkout` 恢复），现改为 `rm -f` + `.tmp` + `mv`；`run.slurm` 的二进制路径修正为
  `../../../bin/`（运行目录多一层），并加了可执行位检查。
- 集群环境：`python3` 3.6.8 无 numpy，需 `module load apps/anaconda3/2023.09`（py3.11 + numpy1.24）；
  程序日志是 `run_out.log`（无 `run.log`）；`sbatch -n` 必须等于该目录 `partation.dat` 的进程数。
- 文档：`docs/程序说明/chap10_build_run.tex` 新增“集群（Slurm）作业与本地/集群同步”小节，
  `93_changelog.tex` 追加 1.6.4 记录；`./build.sh` 重编 0 error。
- 备注：`cases/lidcavity_ac_3blk/control.ec` 仍是 1000 步冒烟配置（此前记忆里的“已改 10000”与实际不符），
  长跑一律 `./prep.sh 10000` 在 `np*/` 内生成独立 `control.ec`。
- 规则补充（用户追加）：**作业测试放集群（`sbatch`，不占登录节点）、后处理可同步回本地、
  每次集群运行的现场控制参数本地也要留存** ⇒ `pull` 增加 `np*/control.ec`、`output_para.out`，
  并改用 `rsync --update`（本轮踩到一次：`pull` 把未推送的 `constraints.md`/`chap10` 覆盖回旧版后重做）。


## 2026-10-18 — 任务 5（i↔j 旋转界面）+ 任务 8（4/6 块锚定回归）+ 左手块（负体积）保护

- **任务 5 的真凶**：`cases/channel_ac_2blk_rot90` 的 `p` 变体（块 2 `i→+y`、`j→+x`）
  早期让 `k→+z`，于是 `(x_i×x_j)·x_k = -z < 0` —— **左手块**。两块覆盖范围正确、
  `|V|` 数值也正常（只是负号），但体积/面积/`dtime_mesh` 全为负 ⇒ 第一步就发散
  （`u` 到 1e10 量级）。修 `gen_case.py`：`p` 变体 `k→-z`（`L=(2,+1,-3)`，bc 的邻居 k
  区间反写 `-2 -1`），`n` 变体 `i→-y`/`k→+z`（`L=(2,-1,+3)`，本来就右手）。
- **加了硬性保护**：`sub_geometry.f90::check_mesh_quality_onemesh` 逐块检查最小单元体积，
  `<0` 先打印块号/最小体积/`(i,j,k)`，再打 “left-handed block(s) in the input mesh …”
  并 `MPI_Abort(3)`；`negative_vol` 为真时不再打 “Check mesh quality OK”
  （原先先打 OK 再报错，自相矛盾）。用 `/tmp/lefthand.py` 造左手网格实测：exit=3、
  `Min Volume=-1.25e-3 at (4,23,1)`，提示语完整。
- **本地复验（400 步）**：`p_np1` 与 `ref_np1` 逐位相同（RMS(u)=1.2e-16），`n` 差 1e-4；
  界面 max|Δu| 3.1e-4（0.026% of u_max）。顺带修了本机 run 脚本的路径坑
  （`mpirun` 用绝对路径 `/home/sundong/Fortran_Project/MultiPhysMixSolver/bin/...`）。
- **集群（wuzh02，4×10⁴ 步）**：`FRESH=1 ./prep.sh 40000` + 5 个 `sbatch`
  （46827434/36/37/39/40）全部 done。`p_np1` vs 参考 **3.0e-16（逐位）**；
  `n_np1` vs 参考 4.0e-4，差值集中在被反转块内部（x>4 RMS 5.6e-4、x<4 4.4e-5、
  界面无尖峰）；np1 vs np2 两变体各 2.4e-10；出口列和 1.000000。
  墙钟 4m19 / 4m22 / 6m01 / 2m51 / 2m49（`-N 1`，含 IO）。
- **两个集群环境坑（本轮踩到）**：
  (1) `sbatch` 的命令行环境里若 `module load apps/anaconda3/2023.09`（为 numpy），
  作业继承 anaconda 的 `mpirun`（MPICH 内核）⇒ 5 个作业在 `MPI_INIT` 前全部
  `*** on a NULL communicator` 秒退（队列 14:59 提交、14:59 全灭）。**提交用干净 shell**，
  用 `/public/software/apps/anaconda3/2023.09/bin/python3` 或事后单独 load 做分析。
  (2) `make` 也不能带 anaconda（`mpif90` 变成 anaconda 包装脚本 ⇒
  `x86_64-conda-linux-gnu-gfortran: command not found`）；干净 shell 下
  `which mpif90 = /opt/hpc/.../hpcx/v2.11.0/gcc-7.3.1/bin/mpif90`，增量重编通过。
- **任务 8 整理**：`cases/lidcavity_ac_46blk/README.md` 补全回归表（m6 np=1/2/3/6 互差
  3.6e-7、np=6 墙钟 5.81 s/加速比 6.78；只换编号 m6a vs m6b 速度差 2.2e-6、压力只差常数；
  Ghia 对照与界面跳变）。
- **文档**：`chap04_lowspeed.tex` 新增“索引旋转（i↔j）的块与左右手约定”（映射表 + `L`/旋向 +
  报错原文 + 实测）；`chap10_build_run.tex` 启动流程加“单元体积符号检查”；
  `93_changelog.tex` 追加 **1.6.5**；`./build.sh` 后再跑一遍 xelatex，64 页、
  0 error、0 未定义引用，`docs/*.pdf` 同步。
- **同步**：`sync_cluster.sh push`（全仓，14:56）+ `pull`（15:14，只取 memory-bank/docs/日志/
  analyze_*）；提交 `482ee66` 后再次 `push`（15:5x）。
- **提交与两端复验**：`482ee66`（src 2 个 + 新算例 2 个 + docs 3 份 tex 与 2 份 PDF +
  记忆库 3 份）。随后：集群 `make clean && make -j4` = **0 error / 2 warning**、
  产出 `bin/MultiPhysMixSolver`（源文件 md5 两端一致 `84e73434…`/`1f32b650…`）；
  本地 `make clean && make -j4` = 0 error / 29 warning（gfortran 13 的既有风格告警）；
  用新编译的本地二进制重跑左手网格 ⇒ exit=3，日志依次为
  `Check mesh quality ......` → `ERROR: block 2 has a NEGATIVE cell volume` →
  `ERROR: left-handed block(s) in the input mesh`（**没有** “OK”），
  `Min Volume=-1.25e-3 at (4,23,1)`；`pdftotext` 校验 PDF 含“索引旋转”与“1.6.5”。
- 新增 `.gitignore`：`cases/channel_ac_2blk_rot90/.gitignore`、`cases/lidcavity_ac_46blk/.gitignore`
  （隐藏 `ref_np1/ p_np1/ …`、`prep_cluster.log`、`analyze_*.txt`，与 3blk 算例同约定）。



