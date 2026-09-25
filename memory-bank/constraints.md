# Constraints — 关键约束（工作流 + 代码约定 + 环境）

## 0. 工作流约束（来自 `.clinerules`，必须遵守）
- **程序说明手册（LaTeX）必须随功能同步更新**（细则见 `.clinerules` 与 `docs/程序说明/README.md`）：
  ① 更新对应章节；② 新增参数同步**附录 A** 与 `control.ec.template`；③ **附录 C 追加更新记录**；
  ④ `cd docs/程序说明 && ./build.sh` 重编译（xelatex 两遍；须 0 error、0 undefined reference）。
  手册源为 ASCII 文件名（`main.tex` + `chap01..11` + `91/92/93`），中文目录名 `docs/程序说明/`；
  PDF 交付物 `docs/OpenCFD-EC-1.16a-程序说明.pdf`。
  LaTeX 注意：`\code{...}` 内必须转义下划线（写 `\_`），数学符号要放进 `$...$`；新增章节记得在
  `main.tex` 里 `\input`。
- 每次新任务开始，**先按序读取**：`memory-bank/projectbrief.md` → `activeContext.md` →
  `progress.md` → `worklog.md`（最近 5 条）→ `constraints.md`；
  然后**先用不超过 10 行**总结（目标 / 已完成 / 当前任务 / 关键约束 / 建议下一步），
  **并问用户是否继续，不要直接改代码**。
- 未读完上述文件前**不要重新扫描整个项目**；**不要重新通读整个代码库**（除非用户明确要求）。
- **不要凭猜测补信息**：缺失就**问用户**。
- 每次任务结束或关键决策后，更新 `activeContext.md` / `progress.md` / `worklog.md`（追加）。

## 1. 代码约定
- 全部 Fortran 90 自由格式、`implicit none`；实数为 `real(PRE_EC)`（=8 字节，`sub_modules.f90`）。
- 文件内**缩进风格不统一**（有的块 5 空格、有的 6 空格）——修改时**先读实际内容**再写补丁，
  不要按记忆猜缩进；跨行锚点匹配容易失败，优先“按行号/唯一子串”定位。
- **面编号**：`1=i-,2=j-,3=k-,4=i+,5=j+,6=k+`（全局统一）。
- **接口码**：11..19 表“块类型对”（不含方向）；`is_interface_bc(bc) = (bc<0) .or. (11<=bc<=19)`；
  耦合分派**按块类型**（`Bn%Block_type`），不由接口码决定方向。
- **连接描述符 `L1,L2,L3`**：本块第 k 维 ↔ 邻块第 `|Lk|` 维，符号=正/反向
  （`Convert_bc`；`sub_convert_inp.f90`）。新写耦合例程应沿用
  `couple_compressible_fluid_solid_face` 的 L 表用法，而非硬编码下标。
- **新增耦合能力时必须向后兼容**：旧分支（如 `(2,5)` 拓扑、码 13/16）代码**逐字保留**，
  只在外层加 `if/else`；改完要能重跑既有算例（porous_fluid_phaseB、high_low_fluid_800K 等）。
- **求解器内部约定**：可压缩块 `U=(ρ, ρu, ρv, ρw, E)`（无量纲，`T*=T/T_inf`，`p*=p/(ρ∞U∞²)`）；
  低速块 `U(1)=ρ=LS_rho`（常密度）、`U(2..4)=速度[m/s]`、`U(5)=T[K]`，压力在 `B%p`；
  多孔块同样 `U(2..4)=速度`、`U(5)=T[K]`，骨架温度在 `B%Ts`。**注意低速块存速度不存动量**。
- 命名沿既有风格（`run_staggered_*`、`couple_*_face`、`*_one_block`、`*_one_mesh`）；
  注释中英混用，新增注释用**中文**、并说明“为什么”。
- 单次编辑上限：`editor` 工具约 6000 字符/次；大改动用“写片段文件 + python 补丁脚本”分块拼接。

## 2. 单位与无量纲化（易错）
- 网格坐标保持**网格单位**（本算例 mm，`Lscale=0.001`）；`x*Lscale` 得 SI。
- 可压缩参考量：`a_ref=sqrt(γR T_inf)`、`U_ref=Ma*a_ref`、
  `RHO_REF = Re·μ∞(T_inf) / (Ma·a_ref·Lscale)`、`p_ref=RHO_REF·U_ref²`。
  ⇒ **`Re` 以 `Lscale` 为参考长度**（`Lscale=0.001` 时 Re=5000 ⇒ 板长 Re_L≈8.3e5）。
- 固体/多孔骨架物性用 SI（K、W/(m·K)、kg/m³、J/(kg·K)）。
- `LS_*`/`Porous_*` 均为 SI（m/s、K、Pa、kg/m³）。
- **`LS_Inlet_Type=2`（质量入口）现按“i 面 + 网格单位面积”实现** ⇒ 用于 j 面入口会
  方向/量纲错；本仓库算例改用等效速度入口 `LS_V_in=G/LS_rho`（见 activeContext 决策 3）。

## 3. 环境与运行约束
- 编译：`cd src && make`（`mpif90`；`-freal-4-real-8` 等来自 `src/makefile`）；
  产物 `src/opencfd-ec1.16a.out`。改 `src/*.f90` 后必须重新 `make` 再验证。
- **交错耦合（`Iflag_Couple_Scheme=1`）必须单进程**（`mpirun -np 1`，目录里不要放 `partation.dat`）：
  跨进程耦合目前只对码 11 实现（`cht_mpi_fluid_solid_face`），**12/19 未实现跨进程**。
- **一个文件同时只能有一个连接（2026-09-25 教训，务必遵守）**：不得对**已打开的文件**用另一个
  unit 再 `open`（违反 F2018 §12.5.6）。本地 gfortran 13 放宽了该检查（同一文件开两个 unit
  `iostat=0`，**本地测不出来**），集群旧 libgfortran 会**致命报错**
  `File already opened in another unit`；若该 `open` 未带 `iostat` 则直接 `Error termination`，
  连捕获都做不到。⇒ 复用调用者已打开的 unit（`rewind` 后读），或先 `close` 再开。
  `control.ec` 只允许 1 处 open（现为 `sub_read_parameter.f90:293` unit 99）。
  **已知未修的同类点**：`util/readflow3d-ver2.5.f90`（99@707 + 96@791）、
  `util/readflow3d-ver2.4a.f90`（99@681 + 96@765）。
- **本地验证通过 ≠ 集群能用**：集群是另一份源码副本（`/work/home/lijunyang/sundong/PorousTest/code/`），
  改完须**同步文件 + 重新 `make`**；工具链版本差（libgfortran）会暴露本地测不到的标准违规。
- 交错模式由 `run_staggered_multiregion` 按块类型对自动分派（11/12/13/19，单配对）；
  `Kstep_Couple_Comp` 为每轮气体段步数，前 `Niter_Couple_Warm` 轮满步、之后**每轮减半**至
  `Kstep_Couple_Min`；`Niter_Couple_Outer` 为外层层数。
- 热流计算**必须有粘性**（`If_viscous=1`）；`Twall<0` / `LS_T_wall<0` 表示绝热壁。
- 交错模式不使用 `t_end`（用 `Niter_Couple_Outer` 控制），`t_end` 给个大值即可。
- 运行产物（`*.vtk`、`flow3d.dat`、`Residual.dat`、`*.inc`、`force*.log`、`mesh-quality*.dat` 等）
  是运行副产物，仓库里大量存在（git 未跟踪）；清理时保留 `run_trial.log` / `iface_*.dat` 作证据。
- 终端偶发“卡住/输出捕获失败”（尤其 heredoc、未闭合引号、后台 make/mpirun）：
  写 python 脚本文件代替 `python3 -c "…"`；长命令加 `timeout`；优先**读日志文件**而非终端回显。

## 4. 禁止 / 谨慎
- 不要为“让算例跑起来”而**改动既有已验证算例的物理行为**（如需改，先记录并回归）。
- 不要凭猜测填 `memory-bank/` 或算例参数；不确定就问用户（尤其：Re 参考长度、
  冷却剂密度/压力水平、hv/dp、是否需要跨进程）。
- 不要提交未验证的改动；提交前至少：`make` 通过 + 相关算例短跑无 NaN/EXIT=0。
- 大文件（`Mesh3d.x` 1.3MB、`flow3d.vtk` ~10MB）不要塞进对话；用脚本读取统计量。

## 5. `cases/fluid_solid` 实测硬约束（2026-09-18 现场测试）
- **低速块（码 12）必须走 AC**：`LS_Algorithm=1`(SIMPLE) 在本算例**第 1 步即发散**
  （`res_p`→6e31、`u`→-1e6）；`AC_CFL=20` 会让 AC 数轮后 **NaN** ⇒
  用 `LS_Algorithm=3, AC_CFL=2, AC_beta=10, AC_Max_Iter=5000`。
- **多孔块预算 = `AC_Max_Iter × Porous_Chunk_Iter`**：实测 `5000×3` ≈170 s/轮（`res_q`≈6e-5）；
  `20000×10` ≈65 min/轮（无必要）。**改动这两个参数前先按 170 s/轮 估时间。**
- **单轮耗时参考（单进程）**：气体段 1000 步 ≈40 s；case1 ≈100 s/全 12 轮（含固体 GS）；
  case2 ≈105 s/轮（LS-AC 主导）；case3 ≈170 s/轮（多孔 AC 主导）。
- **`Iflag_vtk_onefile=1` 每次 `Kstep_save` 写合并 VTK ≈11 MB**（外加 `flow3d.dat` 3.3 MB）；
  集群上建议 `Kstep_save` 取大值或关掉单文件输出。
- **外层判据要与物理时间尺度匹配**：界面温度仍在漂移时（如 case1 固体未达稳态、`T_w` 每轮
  缓降 ~0.1 K），`Tol_Couple_Tw=1e-2` 永远不达 ⇒ 要么放宽判据、要么增外层轮数；
  不要误判为“设置错误”。
- **12/19 无跨进程耦合** ⇒ 必须 `mpirun -np 1`（目录里不要放 `partation.dat`）。

## 6. control.ec 约定（2026-09-23 起：7 组 namelist）
- `control.ec` = **7 组**：`freestream_ec / flow_ec / lowspeed_ec / ac_ec / solid_ec / porous_ec / couple_ec`；
  每个变量**必须在且只在 1 组**；旧组 `$control_ec`（128 变量 + `Solid_*`）保留，
  **新旧并存时新组优先**；组存在但解析失败 ⇒ `stop 1`（不静默回落默认值）。
- **新增参数时必须 6 步**：①`sub_modules.f90` 声明 → ②`set_default_parameter` 给默认值
  （**必须**，这是“算例可省略”成立的前提）→ ③加入**对应分组** namelist → ④**同时**加入旧组
  `$control_ec` 列表（旧文件仍可读）→ ⑤`bcast_para` 广播（`Ipara/rpara` 空闲槽位从 **55** 起）→
  ⑥`output_para.out` 回显。
- 读取顺序：扫描（只有首个非空字符为 `$`/`&` 的行才算组头；`!` 之后为注释）→ 旧组 → 7 组。
- `control.ec.template` 与 `docs/control.ec-说明.md` 必须与代码同步。
- `util/readflow3d-ver2.4a/2.5` 也解析 control.ec：**改分组/加变量后必须同步这两个工具**。
- `cases/**/control.ec` **尚未迁移**（第1档决定）：靠双读兼容，勿擅自批量改写；
  如要精简，只能删“与默认值完全相同”的行，并用 `output_para.out` 逐项 diff 验证。

## 7. 重启 / 节点输出的约定（2026-09-23 起）
- 开关都在 `$flow_ec`：`Iflag_restart`(0=存在即自动续算/1=强制/-1=关)、`Kstep_restart`(<=0→`Kstep_save`)、
  `Iflag_flow_node`(1=写 `flow3d_node.dat`)。**加了开关必须同时进 `bcast_para`**：
  否则各 rank 取值不一致 → 在 I/O 处 MPI 死锁。当前占 `Ipara(55/56/57/58)`
  （58 = `Iflag_Couple_Restart`，2026-09-25 加）。
- `field_restart.dat`：**必须含完整 LAP=4 鬼点缓冲**（`1-LAP:nx+LAP-1`），否则高精度格式续算不稳；
  头记录校验 `nblock/NVAR/LAP/块1维数`，不一致即回退原初始化（`Iflag_restart=1` 时 `stop 1`）。
- 重启/节点输出的 MPI 模式必须与 `output_flow`/`read_flow_data` 一致：
  **只有 0 号进程写/读文件并遍历全局块；其它 rank 只遍历自己拥有的块**（`B_n(m)` 只在本进程有效，
  非 0 号进程访问全局 `B_n` 会越界挂死）。
- 交错模式：三个驱动入口的 `tt/Kstep` 归零已改为“有重启则恢复”；**改写这些驱动时注意
  `call restart_step_output(nMesh)` 必须在 `contains` 之前、且由所有 rank 调用**。
- 输出文件较大：case1 重启 ≈17.5 MB/次、节点文件 ≈2.2 MB/次 ⇒ 集群上注意 `Kstep_save/Kstep_restart`。
- `flow3d_node.dat` 每条记录 = **6 个节点场 `d,u,v,w,T,Ts`**（逐块同长、块序同 `Mesh3d.x`）：
  `T` = 流体温度（固体块 = 固体温度，兼容 `flow3d.vtk` 约定）；
  `Ts` = **固体/骨架温度**（固体/多孔为真实值，LTNE 温差可见；流体/低速块无固相 → 镜像 `T`）。
  **新增/修改块类型分支时必须同时填这两个温度槽位**；校验工具 `util/check_flow3d_node.py`。
- **`field_restart.dat` 版本**：`iver=1`（2026-09-23 起）无耦合信息；`iver=2`（2026-09-25 起）
  在块数据后追加 **coupling-state trailer**（`mode/conv/pair/iter/nf,nk` + `max|dT_w|/tol`
  + 界面量 `T_w,q_w,u,p_w`）。**读侧必须保持对 `iver=1` 与残缺 trailer 的容错**
  （状态视为"未知" ⇒ 旧行为），新增块类型/驱动时同步在 `set_couple_state` 里登记。
- **重启后的耦合模式由 `Iflag_Couple_Restart`（`$couple_ec`）决定**：0 自动（上次已满足
  `Tol_Couple_Tw` ⇒ 切逐步强耦合）/ 1 强制强耦合 / 2 强制交错 / -1 忽略（=旧行为，回归用）。
  交错续算**必须**从 trailer 恢复界面量（`fp_Tw/fp_qw/fp_u/fp_pw`）而不是用 `Twall_Couple_Init`
  重播，否则每次重启都有界面跳变。
- **切到逐步强耦合（`Iflag_Couple_Scheme=0`）必须同时降固体预算**：实测默认
  `Solid_Max_Iter=20000/Solid_Tol=1e-9` 时**每步**固体 GS ≈1400 次扫掠、≈15–20 s/步（case1）；
  现策略为 `Solid_Max_Iter=min(现值,20)`、`Solid_Tol=max(现值,1e-6)`（**只动扫掠上限与容差**，
  `Solid_GS_Omega` 与物理不变），并在日志/`output_para.out` 回显实际值。
  改这段前请先按"每步是否可承受"评估。
- 逐步强耦合续算路径退出前会**补写一次重启文件**（`tight_continuation` 分支），
  否则 `t_end` 落在两次 `Kstep_save` 之间就没有文件可续。
- 详细说明与验证记录见 `docs/重启与节点流场输出说明.md`。
