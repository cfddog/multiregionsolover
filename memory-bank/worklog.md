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



