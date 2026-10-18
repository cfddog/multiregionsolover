# Active Context — 当前任务 / 最近决策 / 下一步

> 更新时间：2026-10-18（**最新 §0r：P3 并行固体 GS 推进（`Solid_Sync_Chunk`）已收口**；
> 上一版 §0q：P1 收口 + 并行验证框架 `cases/parallel` 实测全 PASS**；
> 上一版 §0o：任务 5 旋转界面 + 任务 8 锚定回归 + 左手块保护）。
> **HEAD = `d96c61d`**（`main` 领先
> `origin/main` 15 个提交，尚未 push）
>（本轮 `29eae92`→`10ee6e2`→`6a804fa`→`5094fd8` 四个提交；随后一个记忆库提交）。
> 本次全部改动（control.ec 7 组拆分 + 重启/节点输出 + 三算例迁移 + 文档/记忆库）**已提交并推送**；
> 仅 `.vscode/settings.json`（主题设置）有意未提交。
> 上一版（2026-09-18，fluid_solid 三算例）已归档为 `memory-bank/archive_activeContext_2026-09-18.md`。
> **2026-10-17：程序已更名 `MultiPhysMixSolver`**（旧名 `opencfd-ec1.16a.out` 由 `make` 建兼容软链），
> 详见 §0h；并行改造 P0（跨进程安全网）已完成并提交（§0h 末），P1 起未开始。
> **2026-10-17（最新）：手册 §4.7「多块离散一致性」（§0k）** —— 对流界面单值通量/刻意
> 一阶、粘性两层差分自然等值反号（**不需要**角点数据）、O 型极坐标网格**仍二阶**
> （`cylinder_re40_half` 实测 Cd +0.3%）；手册升 **1.6.2**。
> **2026-10-17：AC 的 `p` halo 统一 + `IF_Debug=3` 缓冲区自检（§0j）** ——
> 三种进程/块分布实测 `max|dU|=max|dp|=0`；手册 1.6.1。
> **2026-10-17：构建布局重构（§0i）** —— 顶层 `makefile` + `src/`（仅源码）+
> `build/`（`*.o`/`*.mod`）+ `bin/`（可执行）+ `cases/`；`make -j4` 已可用。
> 提交：`c8e423d`（更名）→ `c7537f1`（P0）→ `d92f720`（PDF 名引用）→ 本轮（构建布局，v1.6）。
> **本轮（第四轮）改动尚未提交**：`src/sub_IO.f90`、`util/check_flow3d_node.py`、两份 docs、手册与记忆库。
> **第五轮（2026-09-26）改动同样尚未提交**：`Iflag_Couple_Stop` 停机判据（6 个源文件 +
> 2 个 util 工具 + 模板 + 4 份文档 + 记忆库）。
> **第六轮（2026-09-26）改动尚未提交**：区域残差监视输出（`src/sub_debug.f90`、
> `sub_multi_region.f90`、`sub_lowspeed(_ac).f90`、`sub_porous(_ac).f90` + 手册 4 处 + 总结文档 + 记忆库）。

## 0q. 2026-10-18：P1 收口 —— 闭合域 AC 压力基准统一（低速 $+$ 多孔）✅

## 0r. 2026-10-18：P3 并行固体 GS 推进（`Solid_Sync_Chunk`）✅

**目标**：np>1 时各进程持有的固体块扫掠次数互不相同（`|dTs|<Solid_Tol` 各自早退），
只要有跨进程 $T_s$ 交换就会错位挂死 ⇒ 改成「按全局判据同步推进」。

**代码**（`src/sub_multi_region.f90` 为主）：
1. `solid_gs_sweep(mBlock, sweeps, resmax)`：一次完整 GS 扫掠（所有本地固体块，含界面耦合调用）
   ＋ 局部 `|dTs|` 判据，纯局部无 MPI。
2. `solid_solver_one_block`（np=1 或 `Solid_Sync_Chunk<=0`，旧路径）与新增
   `solid_solve_all_blocks`（np>1 且 `Chunk>0`）**共用** `solid_gs_sweep`。
3. 并行路径每 `Chunk` 次扫掠：`MPI_Allreduce(MPI_MAX)`（全局 `|dTs|`）＋
   `MPI_Allreduce(MPI_LOR)`（是否所有进程都收敛），
   终止 `(gmax<=Tol .and. all_conv) .or. sweeps>=Solid_Max_Iter`，
   每外步打印一行 `solid sync (parallel): sweeps= ... (identical on every MPI rank)`。
4. `sub_read_parameter.f90`/`sub_modules.f90`/`control.ec.template`：新增 `Solid_Sync_Chunk`（默认 100）。

**顺带修掉的 4 个 np>1 既有缺陷**：
① `read_solid_bc` 拿**全局**块号索引**局部** `Block(:)` ⇒ 越界（np>1 静默写坏网格，Ts 对不上）；
② `check_remote_coupling_support` 把 Fortran 空语句 `continue` 当 `cycle` 用 ⇒ np>1 任何界面都被判「不支持」；
③ `output_vtk` 的 MPI tag 收发两侧用了不同的块号（局部 vs 全局）＋ 阻塞发送顺序不一致
   ⇒ np=2 在 `MPI_Recv` 处 glibc 堆损坏断言（改为两侧都用全局块号、master 两趟发送顺序一致）；
④ `output_Ts` 用**局部**块号命名 `Ts_block_N.dat` ⇒ 两进程互相覆盖同名文件
   （**P3 最初看到的“np 不一致”其实全是这个输出竞争**，不是求解器差异）；
   现用 `my_blocks(mBlock)`（np=1 时二者相同），并加 `(rank<N>)` 提示 + 注释。

**新回归算例** `cases/solid_4blk/`：4 个**互不耦合**固体块（每块 1-D 导热：i⁻=400 K、i⁺=350 K，
其余对称 ⇒ 精确线性解），np=2 时 blk1/2→rank0、blk3/4→rank1；四块 x 分辨率不同
（41/21/61/15 节点，k=16.3/200/16.3/1）⇒ 各块所需扫掠数 18100/5300/18700/5200 差别极大。
`gen_case.py`（生成 `bc3d_interface.inp` 也必须存在，否则 `read_parameter` 打不开）、
`prep.sh`（np1/np2 运行目录）、`analyze_solid_4blk.py`（对精确解 + 两目录互比）、`README.md`。

**验证**（本地，2 步，`Solid_Tol=1e-9`）：
- np=2 + `Chunk=100`：每步一行 `solid sync (parallel)`，四个块**全部 18700** 扫掠，无死锁；
  `flow3d.dat`/`Residual.dat`/`Step_mess.dat` 与 np=1 **逐字节相同**；
  `Ts` max|Δ|=3.9e−7 K（blk3 逐字节相同；差异只来自“18700 vs 各自扫掠数”的退出时机，
  = `Solid_Tol` 量级）；两者对精确线性解 ≤9.0e−7 K。
- 反证（旧路径）：同一二进制 `Solid_Sync_Chunk=0` 跑 np=2 ⇒ rank0=18100、rank1=18700
  各自早退（相差 600 次扫掠），`Ts` 与 np=1 **逐字节相同** ⇒ 旧路径行为未被改动，
  也说明“若界面跨进程交换 $T_s$，两进程必然错位”这一死锁机理是真实的。
- 编译：本地 `make -j4` OK；集群 `make -j4`（gfortran 7.3.1）EXIT=0。
- 集群复验（`prep.sh 2` + `sbatch -n 1/2`，46850363/46850364，两个作业 30 s 内完成）：
  `analyze_solid_4blk.py np2 np1` 的四个 max|dT| 与本地**一字不差**（3.9e−7/1.5e−7/0/1.7e−7 K），
  `flow3d.dat` md5 `2d0e4d22e55e5dbc617c0b23640e3e03`、`Ts_block_1.dat` md5
  `77d54efb3c867383e942951e81a12bff` **与本地逐位相同**（本地 gfortran 13 vs 集群 7.3.1）
  ⇒ 结果是“算法决定”而非“编译器决定”。现场 `np*/control.ec`、`output_para.out`、
  `Residual.dat` 已 pull 回本地留存。

**遗留（与本任务无关，本轮新发现）**：码 14 固--固**局部**耦合界面温度不连续
（两段 k=16.3 / k=200 串联板，np=1 收敛后界面两侧胞心温度差 ~38 K）⇒ 待查
`couple_solid_solid_interfaces` 的镜像量/面距；P3 验证算例因此改为互不耦合的 4 块。

**文档**：手册 `93_changelog` 1.6.7；`chap06_solid.tex` + `91_appendix_params.tex` +
`control.ec.template` 增 `Solid_Sync_Chunk`。

**背景**：P1（`p` 跨进程 halo ＋ AC 全局判据）的三项代码在 2026-10-17 已落地，但①从未单独收口
（`progress.md` 还写着“P1/P2 待做”）；②复核时发现多孔侧与低速侧规则**互相矛盾**。

**问题（真 bug）**：`sub_porous_ac.f90` 闭合域时是
`if(.not. have_pbc) call ac_zero_mean_p(...)` ⇒ **每个多孔块各减自己的均值**；
低速侧是“只让全局 1 号块锚定，其余块保留水平”。两者在多孔界面等价于加一层常数压差膜
（低速侧 2026-10-17 实测 31→131→324→704 后发散的那条机理）。现有 66 个 `material.in`
无一含 2 个以上多孔块、且多孔算例都有出口（`have_pbc=.true.`）⇒ 一直没暴露，
但 P5（闭合域低速$+$多孔）一做就踩。

**改法**（4 个源文件）：
1. `sub_modules.f90` 新增 `AC_Anchor_Block`（默认 0）；
2. `init_ac_parallel_flags`：在原有 `AC_PBC_Global`（MPI_MAX）之外，统计
   `LS_Algorithm=3` 低速/多孔块中**全局块号最小**者，本进程无 AC 块时用哨兵
   `ANCH_NONE` 参与一次 `MPI_MIN`（所以各进程块数不等也一致）；启动行加打印
   `pressure anchor block=`；
3. 新增 `ac_anchor_pressure(nMesh, mBlock, have_pbc)`（`sub_lowspeed_ac.f90`，紧跟
   `ac_zero_mean_p`）：有压力边界 ⇒ 直接返回；闭合域 ⇒ **只有锚定块**去均值；
   若 `AC_Anchor_Block` 未设置（驱动没调用过 init）则退回旧的“块 1 锚定”。**无集合通信**；
4. 低速（`lowspeed_ac_solver_one_block`）与多孔（`sub_porous_ac.f90`）都改调它。

**向后兼容**：block 1 是 AC 块时行为逐位不变（现有全部算例均属此类，或有出口）。

**验证（本地，改前二进制由 `git archive HEAD` 在 `/tmp/p1old` 构建）**：
| 检查 | 结果 |
|---|---|
| `lidcavity_ac_3blk`（闭合 3 块）np=1，5 步 | `flow3d.dat`/`field_restart.dat`/`Residual.dat`/`Step_mess.dat` **5 个文件 md5 相同** |
| `beavers_joseph/run_ac`（低速$+$多孔，有进出口）np=1，5 步 | 同上 $+$ `Ts_block_2.dat`，**6 个文件 md5 相同** |
| 闭合 3 多孔块对照（lid 网格 $+$ `material 3 3 3` $+$ `porous.inp`，200 步） | 旧版三块均值恒 $0$；新版仅 1 号块为 $0$（块 2 均值 $-1.68\times10^{-3}$）⇒ 规则只在“闭合 $+$ 多孔”上生效 |
| 闭合 3 块低速 `np=1 vs np=2`，**1e4 步**（`analyze_para.py`） | **PASS**：$\max|\Delta u|=1.15\times10^{-6}$、RMS $7.3\times10^{-7}$、压力去常数 $3.4\times10^{-7}$；两侧打印 `anchor block= 1` |
| 集群 `wuzh02` `make -j4`（gfortran 7.3.1） | EXIT=0（4 个改动文件重编 $+$ 链接） |

**新发现（不属 P1，记给 P5）**：闭合域 3 个多孔块、np=2（1$+$2 分布）跑到约 22 个外步后**停滞**；
新/旧二进制**完全一样**（都在 67 行 `Region res` 处停）⇒ 与本次改动无关，是“多孔 AC 在
各进程块数不等 $+$ `AC_MPI_Sync=1`”的既有问题（低速侧 1$+$2 分布已验证可用）；np=1 正常。
另外踩到：`porous.inp` **第一行是条目数**，漏写会让块号/参数整体错位（块退回默认
`eps=0.9/dp=1e-3`）；在算例目录里直接跑会**重写 `bc3d.inp`**（由 `bc3d_interface.inp` 重新生成），
本次已 `git checkout` 复原。

**文档**：手册 `chap04` §sec:ac-halo 增「闭合域的压力基准：全局锚定块（2026-10-18）」，
`93_changelog` 追加 **1.6.6**，`cases/parallel/README.md` 增「给 P3–P6 加测试」第 5 条，
`./build.sh` 重新编译无错误。

**下一步**：P2（`AC_Max_Iter` 跨进程统一步数 / 伪步语义）→ P3（固体块温度 halo）。


## 0p. 2026-10-18：并行验证框架 `cases/parallel`（P3–P6 的回归底座）✅

**做了什么**：新增 `cases/parallel/`，把「np=1 与 np>1 是否同一个解」变成一条命令。

- `run_para_verify.sh`：对 `tests.conf` 的每个测试，为 np=1..k 各建一个运行目录（网格/边界/
  物性软链，**不动算例目录**），用 `gen_control.py` 生成 `control.ec`、`util/mkpart.py` 生成
  `partation.dat`，跑完调 `analyze_para.py` 断言；`MODE=slurm` 时每个目录 `sbatch --wait`。
- `gen_control.py`：base 里已有的键**就地替换**，缺的键按前缀追加到
  `$flow_ec/$ac_ec/$lowspeed_ec/$solid_ec/$porous_ec/$couple_ec`；base 已有同名组时**报错**
  （Fortran 只读第一个同名 namelist 组，静默忽略第二个）。
- `analyze_para.py`：`flow3d_node.dat` 逐点（d,u,v,w,T,Ts）+ VTK 单元中心压力（去常数）。
- 默认阈值 1e-5（速度 max）/5e-6（RMS）/5e-6（压力 rel）/1e-6（T rel）。

**基线结果**（wuzh02，gfortran 7.3.1，10000 步，`AC_MPI_Sync=1` + `AC_Max_Iter=1`）：

| 测试 | np | max\|Δu\| | RMS | 压力 rel | 结论 |
|---|---|---|---|---|---|
| `t01_two_ac`（channel_ac_2blk） | 2 | 4.03e-6 | 2.67e-6 | 2.31e-6 | PASS |
| `t02_four_ac`（46blk m4a） | 2/3/4 | 1.37e-6 | 7.47e-7 | 1.22e-6 | PASS |

**关键认识**（以后别再当 bug 追）：
1. np=1 在伪时间步内刷鬼点（块 Gauss-Seidel），np>1 只能每个外步刷一次（Jacobi），
   **两者只收敛到同一个解**；瞬态差 2000 步时 4.1e-5/8.5e-4，1e4 步降到 1e-6 量级。
2. np>1 之间（np=2/3/4）**逐位相同**（max|Δ|=0）⇒ 跨进程交换本身没有不确定性。
3. `channel_ac_2blk` 自带的 `control.ec` 是 `AC_Max_Iter=20` 的**单进程**配方，在 np=2 下差
   6.6e-3；框架必须在 `tests.conf` 里覆盖成 1。这也是「框架能抓到真问题」的证据。

**踩过的坑**：
- 集群 `LANG=C` ⇒ python3 stdout 是 ascii，脚本里的中文提示抛 `UnicodeEncodeError`，
  4 个比对全被判 FAIL。三个脚本都已加 `sys.stdout.reconfigure(encoding="utf-8")`，
  驱动再 `export PYTHONIOENCODING=utf-8`。
- 驱动里 `while read` 从 `tests.conf` 读配置，子进程（`mpirun`）会继承该 stdin 把剩余配置
  行吃掉（t01 跑完 t02 直接消失）。改成 `mapfile` 先读进内存 + 子进程 `< /dev/null`。
- wuzh02 有**多个登录节点**（login03/login05/…），`/tmp` 不共享；日志必须写在
  `cases/parallel/`（共享盘）里才能跨 ssh 轮询。
- `sync_cluster.sh push` 会用本地旧 `runs/` 覆盖集群结果（`rsync -a` 保留 mtime）⇒ push
  增加 `--exclude 'runs/'`。

**下一步**：P3（固体块）落地后，在 `tests.conf` 加 `t03_*`；`flow3d_node.dat` 增加固体位移/
应力字段时，同步扩展 `analyze_para.py` 的 `VARNAMES` 与断言（脚本按 NVAR=6 校验长度，
字段数变了会报错而不是静默错位）。


## 0o. 2026-10-18：任务 5（i↔j 旋转界面）+ 任务 8（4/6 块锚定回归）+ 左手块保护 ✅

**任务 5：`cases/channel_ac_2blk_rot90`（块 2 转置，界面 = blk-1 `i+` ↔ blk-2 `j-`）**
- 关键发现：**`gen_case.py` 早期把 `p` 变体写成左手块**（`i→+y`、`j→+x` ⇒
  `(x_i×x_j)·x_k = -z`，却让 `k→+z`，`L=(2,+1,+3)`）。两块覆盖范围正确、`|V|` 也正常，
  但体积/面积/`dtime_mesh` 全为负 ⇒ 迭代发散（`u` 到 1e10）。修复：`p` 变体 `k→-z`
  （`L=(2,+1,-3)`，bc 里把邻居 k 区间反序写：`-2 -1` 而不是 `-1 -2`）；`n` 变体
  `i→-y`、`k→+z`（`L=(2,-1,+3)`，本来就右手）。
- 结论（wuzh02，4×10⁴ 步）：`p_np1` 与未旋转参考 `ref_np1`（直接软链
  `channel_ac_2blk` 的网格/bc）**逐位相同**（RMS(u)=3.0e-16、max|Δp|=1.3e-15）；
  `n_np1` 与参考差 RMS(u)=4.0e-4，**差值集中在被反转块内部**（x>4 RMS 5.6e-4，
  x<4 仅 4.4e-5，界面无尖峰）⇒ 同格式在不同索引走向下的舍入/残差地板差，不是映射
  错误；np=1 vs np=2 两变体各差 2.4e-10；出口列和 = 1.000000。
- 墙钟：ref 4m19s、p_np1 4m22s、n_np1 6m01s、p_np2 2m51s、n_np2 2m49s。

**新增强制检查（防复发）**：`sub_geometry.f90::check_mesh_quality_onemesh` 逐块求最小
单元体积，`<0` 先打印块号/最小体积/其 `(i,j,k)`，再汇总
“left-handed block(s) in the input mesh …” 并 `MPI_Abort(3)`（exit=3），不再进入 `Init`；
`negative_vol` 时也不再打印 “Check mesh quality OK”（原先会先打 OK 再报错，自相矛盾）。
实测：本地故意造左手 `p` 网格 ⇒ 启动即 abort（exit=3），报错定位
`(i,j,k)=(4,23,1)`、`Min Volume=-1.25e-3`。

**任务 8：`cases/lidcavity_ac_46blk` 回归表（README 已补全）**
- 进程数：m6 np=1/2/3/6 互差 RMS(u)=3.6e-7（np≥2 之间**逐位相同**），m4 np=1/2/4 = 3.1e-7；
  墙钟 m6 39.4/14.7/10.3/5.81 s（np=6 加速比 6.78），m4 34.8/12.7/9.85 s。
- 压力锚点（只换块编号）：m2a/m2b 1.28e-6、m4a/m4b 1.92e-6、m6a/m6b 2.23e-6；
  `max|Δp| == |mean(p−p′)|` ⇒ 只差一个常数；m6a vs m6c 4.7e-9、m6a vs m6d 逐位相同。
- Ghia(1982)：m2a 0.02058/0.04410、m4a 0.02082/0.04363、m6a 0.02121/0.04311
  （3 块参考 0.02095/0.04365）。界面跳变 5.9e-3–1.9e-2 = “界面刻意一阶”的预期。

**文档**：`chap04_lowspeed.tex` 新增“索引旋转（i↔j）的块与左右手约定”（含三种映射表、
`L`/旋向、保护报错原文、实测数据）；`chap10_build_run.tex` 启动流程加“单元体积符号检查”；
`93_changelog.tex` 追加 **1.6.5**；三个 PDF 重新生成（64 页、0 error、0 未定义引用）。

**集群坑（务必记住）**：`sbatch` 前若 `module load apps/anaconda3/2023.09`（为了 numpy），
作业继承 anaconda 的 `mpirun`（MPICH 内核），与 hpcx 链接的可执行文件不匹配 ⇒ 5 个作业
全部在 `MPI_INIT` 之前报 `*** on a NULL communicator` 立刻死掉。**提交作业用干净 shell**
（`which mpirun` 应为 `/opt/hpc/software/mpi/hpcx/v2.11.0/gcc-7.3.1/bin/mpirun`），
分析时再单独 `module load`/用 `/public/software/apps/anaconda3/2023.09/bin/python3`。
另：`make` 时也不能带 anaconda（`mpif90` 会被解析成 anaconda 的包装脚本 ⇒
`x86_64-conda-linux-gnu-gfortran: command not found`）。

## 0n. 2026-10-18：集群 `wuzh02` 连接打通 + 环境核对 ✅

**做了什么**
- 按 `cases/ssh/集群信息.md`：DNS 只解析出 IPv6（`2409:8728:4eff:2001::100`），
  65091 端口 TCP 可达；两份私钥（`/mnt/c/Users/Administrator/.ssh/lijunyang_wuzh02_rsa`
  与 Downloads 里的 `..._RsaKeyExpireTime_2026-11-12_14-15-23.txt`，**md5 相同**、PEM RSA、
  LF 行尾）⇒ 复制到 `~/.ssh/lijunyang_wuzh02_rsa` 并 `chmod 600`
  （`/mnt/c` 下是 777，ssh 会以 `UNPROTECTED PRIVATE KEY FILE` 直接拒绝）。
- **本机原先没有 `~/.ssh/config`**，`ssh wuzh02` 无法解析 ⇒ 新建 `Host wuzh02`
  （HostName `wuzh02.hpccube.com`、Port 65091、User `lijunyang`、IdentityFile、
  `IdentitiesOnly yes`、`StrictHostKeyChecking accept-new`）。
- 实测：`ssh -o BatchMode=yes wuzh02` → 节点 **`login03`**（第二次落到 `login05`，登录节点不固定）、
  家目录 `/work/home/lijunyang`、`EXIT=0`（免密可用）；首次连接把 ED25519 主机键写入 `known_hosts`。

**集群环境（影响后续同步/编译/排队）**
- `mpif90` = `/opt/hpc/software/mpi/hpcx/v2.11.0/gcc-7.3.1/bin/mpif90`，
  即 **gfortran 7.3.1**（本地 13.x）⇒ 同步后必须重编；老编译器对标准违规更严格
  （当年 `File already opened in another unit` 就是集群先暴露的）。
- 排队 = **GridView + Slurm**（`sbatch/squeue/sinfo`，另有 `bsub/qsub` 兼容）；
  分区 `wzhcnormal`（约 4005 空闲）/`wzacnormal`；`/work` 30 TB 用 93%（剩 2.3 TB）。
- 旧副本 `/work/home/lijunyang/sundong/PorousTest/code/` 仍是**单文件老布局**
  （`opencfd_ec3d_v1.16a.f90` + `*.mod`），没有新的 `makefile+src/+build/+bin/+cases/`
  ⇒ 同步应**新建目录**（如 `~/sundong/MultiPhysMixSolver/`），不要覆盖旧副本。
- `cases/ssh/集群信息.md` 已由 **GBK 转 UTF-8** 并补上上述实测结论。

**作业提交（2026-10-18 记录）**
- 帮助文档：<https://www.scnet.cn/help/docs/mainsite/hpc/cmd/slurm/>。
- 模板在**集群** `~/sundong/run.slurm`（`/work/home/lijunyang/sundong/run.slurm`）：
  `#SBATCH -p wzacnormal / -N 1 / -n 1 / -J MixSolver / -o|e *_out|err.log`，
  末尾 `export PATH=.../MixSolver/bin` + `mpirun -n 1 ../../bin/mixnsolver`。
  **用法：拷到哪个算例目录就在哪 `sbatch run.slurm`**；跑本仓库要把 `mixnsolver` 换成
  `../../bin/MultiPhysMixSolver`、`-n` 改成目标 np，`export PATH` 可删。
- ⚠️ **本仓库本地根目录没有 `run.slurm`**（`find` 全盘无、也不在 `.gitignore` 里），
  模板只在集群侧；`~/sundong/MixSolver/MixSolver/`（`mixnsolver*`）是**另一个项目**，
  `~/sundong/PorousTest/code/` 是单文件老布局 ⇒ 本仓库**尚无集群副本**。

**未做（下一步候选）**：`git push`（15 个提交）、用 `rsync -az` 首次同步整仓 + 集群端
`make -j4`、按 Slurm 提交 3 块 AC 算例的 np=1/2/3 复验作业。

## 0m. 2026-10-17（收尾提交）：3 块 AC 修复入库 + 冒烟复验 ✅

**提交（`main`，本轮 4 个）**
- `29eae92` fix(ac)：界面通量符号 + `is_high` 判据、`AC_MPI_Sync`、`p` halo、
  `IF_Debug=3` 缓冲自检、`init_ac_parallel_flags`、`read_solid_bc` 启动死锁。
- `10ee6e2` cases(lidcavity_ac_3blk)：输入 + `gen_single42.py` + `analyze_3blk.py`
  + `prep.sh` + README + 目录内 `.gitignore`（实测目录不入库）。
- `6a804fa` docs(manual)：手册 **1.6.3**（符号约定段 + `tab:ac3blk` + changelog）。
- `5094fd8` chore：删除根目录遗留 `control.ec`（输入按算例分发）。

**提交前修正（都是为了"文档=代码"一致）**
- 删除 `ac_res_norm` 中 `IF_Debug=1` 的临时 `DBG2` 统计块；`ac_dbg=1` 的界面诊断保留
  （默认 0）。
- `cases/lidcavity_ac_3blk/control.ec`：`AC_Max_Iter` 50→**1**、`t_end`/`Kstep_save`
  1000/250→**10000/10000**，与实测配置一致。
- `control.ec.template`：闭腔 + np>1 **不再报错**，改为"压力水平由全局 1 号块锚定"
  （代码 `init_ac_parallel_flags` 已改，模板注释是旧的）。
- 冒烟（新二进制，2000 外步）：np1/np2 均无 NaN、u∈[−0.207,1]，
  `np2 vs np1 max|Δu|=8.7e-4`（瞬态量级，收敛后 1e-6）；`make -j2` EXIT=0。

**未提交（有意）**：`.vscode/settings.json`（主题）、`cases/ssh/`（集群信息）、`workplan.md`。

## 0l. 最新（2026-10-17）：3 块 AC 方腔「界面通量符号 + is_high 判据」两处致命 bug 修复 ✅

**现象**：`cases/lidcavity_ac_3blk`（Re=100，42×42，三块等分 x，界面 x=1/3、2/3）在任何
进程数下都不收敛：块 1 的 `i=nx-1` 与块 3 的 `i=1` 出现棋盘模态并指数增长（增长因子
≈1.3/步）；即便压住发散，界面仍整列错位（u 达 1.68/1.83）。**np=1/2/3 结果逐位相同**
⇒ 先判为离散 bug（而非竞争 bug），并靠 np=1 隔离验证定位。

**根因（`ac_face_flux_rus1_own`，`src/sub_lowspeed_ac.f90`）**
1. 两侧都做 `RAC += FLX*A`，而两侧传入的面法向**相同**（i 面恒 +x），Rusanov 耗散项
   −λ/2·(q_R−q_L) 在界面上相互抵消 ⇒ 界面棋盘模态零阻尼。正解：低侧 `RAC -= FLX*A`、
   高侧 `RAC += FLX*A`。
2. `is_high` 判据写反（把"ghost 索引更大"当成高侧），中央通量两侧同时反向 ⇒ 界面变抽吸源。
   正解：按 own/ghost 索引关系判（`iG>iO`、`jG>jO`、`kG>kO` ⇒ 低侧），与块号、进程号无关。

**验证**
- np=1/2/3（10^4 步）：`max|Δu|=1.16e-6`、`max|Δv|=2.4e-7`、`max|Δp|=9.8e-7`；
  np=2 与 np=3 **逐位相同**（进程数无关）。
- 3×10^4 步收敛（`res_q≈3e-14`），u∈[−0.233, 0.926] ⊂ [−1,1]；vs Ghia RMS：u 0.02095 /
  v 0.04365（与已验证的单块 40×40 答案 0.02095 一致）。
- 同网格单块 42×42 对照（新增 `gen_single42.py`，目录 `s42`）：`max|Δu|=9.3e-3`、
  rms 3.0e-3 —— 与界面**刻意一阶**的预期误差同量级（h=1/42、∂u/∂x≈0.4 ⇒ ~5e-3）；
  离界面 3 胞以外两者逐位一致。单块 42×42 vs Ghia：u 0.02013 / v 0.04387。
- 关键约定：**`AC_Max_Iter=1` 是这类多块算例的要求**（np=1 同伪步刷新与 np≥2 外步末刷新
  的耦合节奏才一致；`>1` 时界面滞后于块内）。

**文档**：`chap04` 新增「两侧各算路径的符号约定」段（含两条 bug 说明）+ 表 `tab:ac3blk`；
`93_changelog` **1.6.3**；`main.tex` 1.6.3；`./build.sh` **0 error**。

**下一步**：其余待办不变；可选：把「3 块 vs 单块」的界面误差做成 h 减半的收敛性表。

## 0k. 最新（2026-10-17）：手册 §4.7「多块离散一致性」+ 曲面网格适用性 ✅

**审阅提问**：① 多块 AC 的对流/粘性离散是否需要特殊处理？二阶 Laplacian 为什么不需要
角点数据？② 这种"三个一维二阶差分"在不横平竖直的网格（圆柱 O 型）上会不会失效？

**核对结论（源码只读）**
- 对流 `ac_face_flux`：法向 5 点 `a..e`（`AC_Recon=1` 只用 `b..e`/1 层 ghost；WENO5 用
  `a`/2 层）。内面循环 `ii=2..nx-1` **排除界面面**（防双计）；界面面用 `ac_ifc_flux_own`
  （`FSH` owner 单值，全局块号定 owner，`ac_ifc_src_cell` 复用缓冲交换的 `L1/L2/L3` 映射）
  或跨进程回退 `ac_face_flux_rus1_own`，**两者都只用一阶 Rusanov**（界面 ghost 每外步才刷
  新，二阶重构会激发 checkerboard）。
- 粘性 `ac_viscous_res`（`sub_lowspeed_ac.f90:1830--1893`）：`f=1..nx` **不排除**界面面，
  `if(f-1>=1)/if(f<=nx-1)` 只累加到真实内胞；`H=ν(u_f−u_{f-1})/d`，
  `d=0.5|x(f+1)−x(f-1)|`，斜网格取 `dproj=|Δx·n|`（>0.2d 才用）。两侧分子相同 ⇒
  `d_A=d_B` 时**恰好等值反号**，界面动量守恒，**不需要** `FSH`。壁面 `d=Δn=2d_wall` ⇒
  `H=−νu_1/d_wall`（半个胞心距口径，二阶）。
- 只含 `∂xx+∂yy+∂zz`、**无混合导数** ⇒ 7 点面邻居壳层闭合，ghost-角点从不被读取；
  写成应力张量形式最坏也只需 ghost-**edge**。

**曲面/非正交（不会失效，三条前提）**：① 速度存**笛卡尔分量**（本代码如此）⇒ 逐分量标量
Laplacian 与坐标系无关，不需要柱坐标曲率项；② ν 常数（`LS_mu`）；③ 法向导数用"索引线差分
÷ 法向投影"（过松弛法）：正交网格严格相等，O 型网格 `i`=径向 / `j`=周向时 skew =
`O(Δθ)=O(h)` ⇒ 仍**二阶**；FV 用真实面积/体积、无显式 `1/r` ⇒ **无坐标奇点**。
实测 `cases/cylinder_re40_half`（O 型极坐标，120×80 / 160×120，Re=40，AC）：`Cd`
1.5706→1.5536 vs Rogers&Kwak 1.549（**+0.3%**），`Cd_p` −0.2%，壁面积分↔远场动量平衡
互校（1.5706↔1.5708）。会降阶：恒剪切网格（非正交度不随 h 减小）→ 一阶；`μ(T)` 未实现；
界面非共形/间距跳变 → 界面不守恒。

**文档**：`chap04` §4.7（表 17 依赖表、对流三件事、界面粘性等值反号表、角点判据、曲面网格
适用性、表 18 实测）；`93_changelog` 1.6.2；`main.tex` 1.6.2；`./build.sh` 0 error /
0 undefined（第二遍消解前向引用；PDF 无 `??`，新增内容在 pp.18–20）。

**下一步**：见 §0j 末（`ac_zero_mean_p` 全局均值、`AC_MPI_Sync=0 & np>1` 升级/报错、`FSH`
跨进程回退）；另：若将来上 9 点/紧致/顶点型格式，**必须先扩 ghost 角点交换**（当前只填面
slab）。

## 0j. 最新（2026-10-17）：AC 的 `p` halo 统一 + 缓冲区交换自检（`IF_Debug=3`）✅

**问题（用户指出）**：`p` 不在 `U(1:5)` 内。同进程 `p` 的 ghost 原先由
`Umessage_send_mpi` 里一个“融合块”拷贝，条件只有 `LOWSPEED↔LOWSPEED`，且依附于
`U` 的调用点；`p_send_mpi` 又只在 `AC_MPI_Sync==1` 分支里调用 ⇒ 缺省
`AC_MPI_Sync=0` 时 np=1 有新鲜 `p`、np>1 永远不刷新；`AC_MPI_Sync=1` 时同进程
`p` 被拷两遍；多孔驱动在 np=1 时完全不刷新 halo。

**方案 A（已实施）**
1. 删除融合块：`p` 只走 `p_send_mpi`/`p_recv_mpi`。同进程直写（无临时数组、无 MPI），
   跨进程才打包 `MPI_Bsend`；两侧过滤条件一致（`LOWSPEED|POROUS` × `LOWSPEED|POROUS`）。
2. 驱动条件对齐：`AC_MPI_Sync==1 .or. Total_proc==1` 时**同时**刷新 `U` 与 `p`
   （`sub_lowspeed_ac.f90`、`sub_porous_ac.f90`，后者新增 np=1 分支）。
3. 新增 `check_buffer_exchange` + `IF_Debug=3`：哨兵值逐 ghost 胞**严格**比对，
   非破坏性（保存/恢复）；同进程直接比，跨进程把“邻居 ghost 应有内容”发给对侧再比。

**实测**（`/tmp/acbuf`，`t_end=20`，`AC_MPI_Sync=1`，三例 EXIT=0）

| 分布 | 算例 | 子面 同/跨 | ghost 胞 | max\|dU\| / max\|dp\| |
|---|---|---|---|---|
| 1 进程 2 块 | `channel_ac_2blk` | 2 / 0 | 320 | 0 / 0 PASS |
| 2 进程 各 1 块 | `channel_ac_2blk` | 0 / 2 | 320 | 0 / 0 PASS |
| 2 进程 各 2 块 | `channel_ac_4blk` | 4 / 4 | 960 | 0 / 0 PASS |

**踩坑**：① 自检第一版用 `nb1`（接收块号）算期望值 → 必然差 ~1e9；
② 跨进程**不能** `B1 => Mesh%Block(B_n(nb1))`（`B_n` 只对本进程块有效，会指到错块）；
③ `size()` 已含 NVAR 维，`NVAR*size()` 使 `MPI_Recv` 声明 5 倍缓冲 → 堆损坏（malloc assert）。

**未做 / 下一步**：`ac_zero_mean_p` 仍是**块局部**均值（无压力 BC 的多块 AC 会出现块间
`p` 偏置，建议改为全局均值 `MPI_Allreduce`）；`AC_MPI_Sync=0 & np>1` 伪步内仍不刷新（仅告警）；
P2 的“所有块统一步数”只在 `AC_MPI_Sync=1` 下实现。
文档：`chap04` 新增 §`ac-halo`/§`bufchk`、附录 A 补 `AC_MPI_Sync` 与 `IF_Debug=3`、
`93_changelog` 1.6.1，PDF 已重编译（0 error）。


## 0i. 最新（2026-10-17）：构建布局重构 `makefile + src/ + build/ + bin/ + cases/` ✅

**动机**：`src/` 里混着 44 个中间文件（`*.o`、`*.mod`、可执行、20 份编译日志、旧二进制备份），
"源码目录"名不副实；编译产物/可执行文件也没有统一位置。

**新布局**

| 路径 | 内容 |
|---|---|
| `makefile`（顶层） | 唯一构建规则：`make` / `make -j4` / `make clean` / `make info` |
| `src/` | **只放源码**：37 个 `*.f90` + `makefile`（转发入口，`cd src && make` 仍可用） |
| `build/` | `*.o` / `*.mod`（`-Jbuild -Ibuild`），可随时 `rm -rf` |
| `bin/` | `MultiPhysMixSolver` + 兼容软链 `opencfd-ec1.16a.out` |
| `src/modules.dep` | 文件级模块依赖（36 条），由 `util/gen_moddep.sh` 扫描 `module`/`use` 生成 ⇒ **`make -j4` 首次构建即安全** |
| `cases/` | 脚本/README 运行路径统一 `../../bin/MultiPhysMixSolver`；`run_limits.sh` 用 `$R/bin/...`；`.vscode` 的 build/run/debug 同步 |

**关键实现**
- 顶层 makefile 用 `$(BUILDDIR)/%.o: $(SRCDIR)/%.f90` 模式规则 + `-J$(BUILDDIR) -I$(BUILDDIR)`；
  链接产物写 `bin/`，随后 `ln -sf` 建旧名软链。
- 依赖**不用** gfortran 的 `-MMD`：其输出要求 `-cpp`（会改变编译语义），改用
  `util/gen_moddep.sh`（awk 扫 `module`/`use`，排除自有模块/`iso_*`）生成 `src/modules.dep`，
  该文件入库随源码一起维护。**改 `use` 关系或增删源文件后必须重跑该脚本。**

**验证（本地）**
- `make clean && make -j4` → **EXIT=0，89 s**（此前从干净树 `make -j4` 必失败）。
- 数值惰性：`channel_ac_2blk`（2 块 + AC，5 步）`flow3d.dat`、`field_restart.dat`
  与重构前二进制 **逐字节相同**。
- P0 安全网负例（np=2 + 码 12 跨进程）仍 **EXIT=2**；np=1 正常 EXIT=0。
- 手册 v1.6：`chap01`（目录树/上手流程）、`chap09`、`chap10`（编译/运行/集群注意）、
  `chap11`、`93_changelog` 1.6；`./build.sh` 0 error。

**顺带清理**：`src/` 下 44 个中间文件、根目录 2 个 `*.mod`、`cases/**` 里 20+ 个旧
`opencfd-ec1.16a.out` 二进制副本（`*.out` 已在 `.gitignore`）。`.gitignore` 新增
`/build/`、`/bin/`、`*.d`。


## 0h. 最新（2026-10-17）：程序更名 **MultiPhysMixSolver** ✅（后续并行改造从 P0 开始）

**动机**：多物理（可压 / 低速 AC / 固体 / 多孔）混合求解，旧名 `OpenCFD-EC-1.16a` 已不达意。
**原则**：更名**数值完全惰性**；旧名保留为兼容软链，旧算例目录零改动。

**改动清单**
| 类别 | 内容 |
|---|---|
| 源码 | `git mv src/opencfd_ec3d_v1.16a.f90 src/multiphysmixsolver.f90`；文件头注释 + 启动横幅（5 行，含 “based on OpenCFD-EC 1.16a”） |
| 构建 | `src/makefile`：`srs` 换名、`EXE = MultiPhysMixSolver`、链接后 `ln -sf $(EXE) opencfd-ec1.16a.out`、`clean` 同时删两者 |
| 输出 | `sub_read_parameter.f90:383` 的 `output_para.out` 版本行改为 `MultiPhysMixSolver (based on OpenCFD-EC Ver 1.1t)` |
| **未改（有意）** | `flow3d.dat` / `Ts_block_*.dat` / `bc3d.inc` 等**数据文件头字符串仍为 `OpenCFD-EC`**（后处理兼容）；`field_restart.dat` 是纯数值记录，与文件名无关 |
| 算例 | 22 个 `cases/**` 的 README / `gen_case.py` / `run_*.sh` 中可执行名改为 `MultiPhysMixSolver`（`run_limits.sh`、`run_grid.sh` 等按脚本位置推导 `src/` 的也一并改） |
| 编辑器 | `.vscode/tasks.json`（run MPI）、`.vscode/launch.json`（gdb program/name） |
| 文档 | 手册 `main.tex`（标题/页眉/摘要）、`chap01`、`chap09`、`chap10`、`93_changelog.tex` 1.5；`docs/程序能力与算例考核总结.md`、`docs/重启与节点流场输出说明.md`、`cases/beavers_joseph`、`cases/porous_fluid_phaseB/staggered_tw` 的源码引用；`build.sh` 产物名 |
| 产物 | `docs/MultiPhysMixSolver-程序说明.pdf`（新）+ `docs/OpenCFD-EC-1.16a-程序说明.pdf`（同名旧文件同步更新，保留兼容）；49 页、0 error |

**验证（本地）**
- `cd src && make clean && make` → EXIT=0，产出 `MultiPhysMixSolver` + 软链 `opencfd-ec1.16a.out -> MultiPhysMixSolver`。
- **逐位一致回归**：用 `git show HEAD:...` 复原更名前源码另建 `/tmp/oldsrc` 编译旧二进制，
  与新版在 `cases/channel_ac_2blk`（**2 块 + AC，5 步，`Kstep_save=5`**）同输入对跑：
  `flow3d.dat`、`field_restart.dat`、`Residual.dat`、`force.log`、两块 `*.vtk`、`Step_mess.dat`
  **全部逐字节相同**；唯一差异是 `output_para.out` 第 1 行版本串（预期）。
- `docs/程序说明/build.sh` 两遍通过，0 error / 0 undefined reference。
- ⚠️（**2026-10-17 已解决，见 §0i**）旧 `src/makefile` 无模块依赖声明，`make -j4` 会竞态失败
  （`Cannot open module file`）；现由顶层 makefile + `src/modules.dep`（`util/gen_moddep.sh`）
  提供文件级依赖，`make -j4` 从干净树也能通过。

**下一步（并行改造，按 P0→P6 逐级验证）**
P0 跨进程安全网（12/19/14/18 远端界面计数 + `MPI_Allreduce` + `MPI_Abort(2)`）→
P1 `p` 的跨进程 halo（`p_send/p_recv_mpi`）+ AC 全局判据 → P2 驱动层逐步 halo →
P3 固体扫描数统一 → P4 码 12 → P5 码 19 → P6 码 14/18 + 手册《并行能力矩阵》；
每级都在 `cases/parallel/` 下用 np=1 vs np=2 的 `max|d|=0` 验收。详见 `workplan.md` 与
`cases/parallel/`（待建）。

### P0（2026-10-17，已完成）：跨进程安全网
- 新增 `check_remote_coupling_support(nMesh)`（`src/sub_multi_region.f90`，文件末尾）：
  遍历全部界面，统计**邻居在别的进程且该组合没有 MPI 版本**的条目
  （12 可压↔低速、14 固↔固、17 固↔多孔、18 多孔↔多孔、19 可压↔多孔），
  `MPI_Allreduce(MPI_SUM)` 后 `MPI_Abort(2)` 并打印清单 + 修复建议（改 `partation.dat` 或 `-np 1`）。
  `Total_proc<=1` 直接返回 ⇒ np=1 零影响。
- 调用点：`src/multiphysmixsolver.f90` 在 `Init_flow` 之后、时间推进之前
  （顺序：`Init` → `set_control_para` → `check_mesh_quality` → `Init_flow` → **自检** → `Start`）。
- 有意**不**在耦合例程里累加计数（那样要等第一步算完才报错），启动自检失败更快、
  且不会在交错驱动里重复计数。
- 验证：
  * np=1 `channel_ac_2blk`（2 块 AC，5 步，`Kstep_save=5`）与 P0 前二进制
    `flow3d.dat` / `field_restart.dat` **逐字节相同**；
  * np=2 + `cases/fluid_solid` deck 把 `material.in` 改成 `0 2`（可压↔低速跨进程）
    → **EXIT=2**，日志打印 `ERROR: cross-process coupling is NOT implemented for 12(...)`
    `remote interface entries = 2`；
  * np=2 + 原 `fluid_solid`（多孔 3 ↔ 低速 2，界面 16，`ls_por_iface_mpi`）
    → 越过自检并进入 `AC porous solver` 迭代（无假阳性）。
  * 手册 `chap07` 新增《跨进程安全网》小节 + 1.5.1 更新记录。

**下一步：P1**（`p` 跨进程 halo：镜像 `Ts_send/_recv_mpi` 写 `p_send/p_recv_mpi`；
`ac_have_pressure_bc` 用 `MPI_LOR` 全局化；`ac_zero_mean_p` 用全局均值），
验收 `cases/parallel/t02_ac_2blk` np=1 vs np=2 `max|d|=0`。

## 0g. 最新（2026-10-17）：第四～七轮改动分组提交 + ε 极限“负结果”澄清 ✅

**提交（HEAD 由 `ff7066a` → 三个提交）**：
| 提交 | 内容 |
|---|---|
| `c974e81` feat(output) | `flow3d_node.dat` 3 条 PLOT3D 头 + `Iflag_restart<0` 冷启动确定性 + `check/compare_flow3d_node.py` |
| `ebe7053` feat(solver) | `region_stat_mod` 区域残差监视 + `Iflag_Couple_Stop` 停机判据 + 多孔 ε 极限 + 模板/读参/util 同步 |
| `a29a6ee` cases(limits) | `cases/fluid_solid/limits/` 5 个短调度算例 + `run_limits.sh` + README（含实测表与负结果分析） |

**头记录的真实历史（重要，勿再误记）**：HEAD 里的 `output_flow_node` **原本没有头**；
工作树一度写成**过渡的 2 条**（块数 + 每块 `ni,nj,nk,nvar`），2026-09-26 的 limits 文件
都是这个布局；现恢复为**标准 3 条**（块数 / 各块维数 / 变量数=6）。⇒ 任何读端必须同时
兼容 2 条与 3 条；`util/check_flow3d_node.py` 已按“记录数 + 长度”自动判别并回显
`format=... 2-record variant` 或 `3-record PLOT3D-Q header`。

**ε 极限实测（2026-09-26 一轮，2026-10-17 复核并写入手册 §11.4 / README）**：
- ✅ 5 例 EXIT=0、NaN=0、`STRUCTURE OK`；ε=0 的块1 `u,v,w` 在**全部节点精确为 0**；
  ε=1 的阻力自零（ε=1 有 O(10) m/s 流动，ε=1e-6 时 |u|≲1e-9）。
- ❌ ε=1 ↔ 低速参照不逐点一致：主分量 **u（±17.6）** vs 参照 **v（17.7）**；
  `--blocks 1`：max|Δu|=17.64、max|Δv|=28.02、max|ΔT|=96.34 K，max 位置在
  (59,5,1)/(12,2,1)/(65,2,1)——**排除界面首层后指标完全不变** ⇒ 差值不只来自界面一排。
- ❌ ε=0 ↔ 固体参照温度不重合：max|ΔT|=498.6 K（Ts 303 vs 800 K），两 deck 仅接口码
  19 vs 11 不同 ⇒ 这一步实际在检验“界面壁模型”，不是块内导热。
- ❌ 多孔块能量几乎不推进：三例流体 T 均停在 300 K、Ts 最多 +3 K（外层 2 轮 × 块内 2 次）。
- 结论口径：limits 只是**回归快照**；要拿到定量证据需 ① 同接口模型（12↔12 或 19↔19）
  ② 外层/多孔预算提高一个量级（或跑到稳态）。

**未提交（有意）**：`.vscode/settings.json`（主题）、`workplan.md`（用户笔记，**GBK 编码**，
未转 UTF-8 前不建议入库）；6 个新 `Mesh3d.x` 单独提交（`cases: 补齐网格输入`）。

**下一步**：① 集群同步（`/work/home/lijunyang/sundong/PorousTest/code/`）后 `make`，
重点看旧 libgfortran 的双 unit 检查与 `-freal-4-real-8` 行为；② 若要真做 ε 定量对照，
新增“ε=1 + 低速入口（码 5/lowspeed_inlet_velocity）”或“同接口 12↔12”的 deck；
③ 三例 `cases/fluid_solid/*/README.md` 仍是旧待办。

## 0f.（2026-09-26 第六轮）区域残差监视输出（固体/低速/多孔输出精简）✅

**需求（用户提出）**：「固体、低速和多孔介质减少屏幕输出，输出区域最大 rms、平均 rms 以及最大 rms 的位置」。
用户选定**方案甲**：不加 namelist 参数，硬编码「每 100 次迭代 1 行 + 结束汇总」。

**改动（6 个源文件）**：
- `src/sub_debug.f90`：新增模块 `region_stat_mod`（常量 `RES_PRINT_STEP=100` + 例程
  `region_res_stat(nMesh, mBlock, what, kk, res_max, imax, jmax, kmax, res_rms, a1, a2)`），
  只由 0 号进程打印一行：`Region res [名] block n (块类型) step=/final= k : max= … rms= … at (i,j,k)= … xyz=[…] m`；
  `kk>0` 过程行 / `kk<0` 结束汇总行；`what` 选残差名（1 固体 GS、2/3 低速/多孔 AC、4/5 SIMPLE、6 AC 能量、7/8 LTNE）；
  行尾保留旧字段 `res_q/res_m/res_p/dT(window)`。坐标用格心 `%xc/%yc/%zc`×`Lscale`。
- `src/sub_multi_region.f90`（固体 GS）：循环内同时累加 `Σr²` 与 argmax；每 100 扫掠 1 行；
  结束 1 行；删除 `GS iter`（前 5 次）与 `DBG: Ts(...)`；`Couple_Res_Solid` 登记保持。
- `src/sub_lowspeed.f90` / `src/sub_porous.f90`：`*_pressure_correction` 追加返回
  `imax,jmax,kmax,res_p_rms`（在既有 `|div F|` 循环内累加）；驱动每 100 迭代 1 行 + 结束行；
  删除 `DBG iter/maxU2/div max`、`DEBUG T iter`、`SIMPLE iter`、末尾 `SIMPLE iterations` 行。
- `src/sub_lowspeed_ac.f90`：新增 `ac_res_stat`（口径与 `ac_res_norm` 一致，返回合并幅值的
  max/位置/rms + `res_q/res_m`）；AC 结束打 1 行汇总并删掉重复的 `iter==1` 打印；
  `ac_lowspeed_energy` 在 50 步窗口内累加统计、删掉每 1000 步进度行、结束打 1 行。
- `src/sub_porous_ac.f90`：AC 结束打 1 行汇总（删 `iter==1` 行）；`ac_porous_ltne` 在 100 步窗口内
  统计流体/骨架两组量、删掉每 2000 步进度行、结束各打 1 行。

**口径（与可压块 `Mesh%Res_max/Res_rms` 一致）**：`max`=逐胞残差幅值区域最大值；
`rms`=`sqrt(Σr²/N)`（N = (nx-1)(ny-1)(nz-1)）；`at (i,j,k)`=max 所在胞 + 胞心坐标 [m]。
**SIMPLE 路径**（`LS_Algorithm≠3`）用 `|div F|`（`res_p`）；**AC 路径**过程行仍按 `AC_Print`（默认 500）。

**验证（`mpirun -np 1`，全部 EXIT=0）**：
| 用例 | 覆盖 | 结果 |
|---|---|---|
| `cases/solid_1d`（t_end=0.3） | 固体 GS | 每 100 扫掠 1 行 + `final=`；`max` 与旧 `final res` 一致（3.9e-7@20000） |
| `cases/channel_ac_4blk` | AC 低速 4 块 | 4 条 `low-speed AC flow` 行，行尾 `res_q/res_m` 与 `AC iter` 自洽 |
| `cases/porous_ltne_1d`（SIMPLE，iters↓） | 多孔 SIMPLE | `porous SIMPLE div(F)` 过程行 + `final=` 单调下降 |
| `cases/porous_ltne_1d/run_ac` | AC 多孔 + LTNE | `porous AC flow` + `LTNE fluid/frame` 三行 |
| `cases/porous_plug_ac/eps06` | AC 多孔（hv=0） | `porous AC flow` 汇总行 |
| `fluid_solid/run_case1_solid`（码 11，2 轮） | 耦合内固体 | 轮内固体汇总行正常，退出信息不变 |
| `fluid_solid/run_case2_lowspeed`（码 12，1 轮） | 耦合内 AC + AC 能量 | `low-speed AC flow` / `AC energy` 行 + 接口度量行 |
| `fluid_solid/run_case3_porous`（码 19，1 轮） | 耦合内多孔 AC + LTNE | `porous AC flow` + `LTNE` 两行，`porous chunk done` 正常 |
临时目录 `/tmp/rt1/*`；`cd src && make` EXIT=0。

**文档**：手册 `chap06`（固体日志小节改写为统一格式）、`chap04`/`chap05` 各新增
「残差监视输出」小节、`92_appendix_codes.tex` 日志关键字表（新增 `Region res ...` 行，
替换已删除的 `Solid solver Block ... GS iterations:` 行）、`93_changelog.tex` 新增 **1.3** 行；
`./build.sh` 重编译 **45 页、0 error**。`docs/程序能力与算例考核总结.md` 新增 §2.9。

**记住（手手相传）**：① 打印间隔是硬编码常量 `RES_PRINT_STEP=100`（`region_stat_mod`），
要改节奏只改这一处；② 调用 `region_res_stat` 的例程必须 `use region_stat_mod`（否则链接期
`undefined reference`，本轮踩过一次）；③ 长稳态算例（如 `solid_1d` 20000 扫掠）会出现
~200 行趋势表——若嫌多，把 `RES_PRINT_STEP` 调到 1000 或改成“decade 式”即可；
④ 收敛与监视节奏重合时用局部 `conv_ok` 跳过 `step=` 行，避免与 `final=` 汇总行重复。

## 0e. 最新（2026-09-26 第五轮）：交错外层停机判据 `Iflag_Couple_Stop` ✅

**需求（用户提出）**：「以 `t_end` 为结束标志，或者以全局所有的残差都小于残差标准」——
即交错耦合（`Iflag_Couple_Scheme=1`）原先既不看 `t_end`、也不看求解器残差，只能靠界面量判据
（`Tol_Couple_Tw/p/u`）或跑满 `Niter_Couple_Outer` 收口。

**改动（用户选定方案 A：新开关控制，默认 0 保留旧行为）**：
- `src/sub_modules.f90`：新 `Iflag_Couple_Stop`（默认 0：旧行为；1：每轮末测试
  `tt>=t_end` **或** 全局残差全达标）+ 3 个每轮残差登记量 `Couple_Res_AC_q/AC_m/Solid`。
- `src/sub_multi_region.f90`：新 `couple_stop_test`（可压块 RMS 现场统计 + MPI 归约；按“块类型是否存在”
  逐类比对 `Res_Inner_Limit` / `AC_Tol` / `Solid_Tol`；`-1.` 表示该类块不存在）；三个交错驱动
  （码 11/13、12、19）在轮末调用它，并新增：轮首清零残差登记量、气体段按剩余时间**截断**
  （`nstep=min(nstep, ⌊(t_end-tt)/dt_global⌋+1)`）、退出提示区分“到 `t_end`”与“残差达标”。
- 残差登记：`sub_lowspeed_ac.f90` / `sub_porous_ac.f90`（AC `rq/rm` 取本轮 `max`）、
  `sub_multi_region.f90::solid_solver_one_block`（GS 末次 `res`）。
- `src/sub_read_parameter.f90`：6 步流程（默认值 / `$couple_ec` / 旧 `$control_ec` / `Ipara(59)` /
  `output_para.out` 回显）；`control.ec.template` 同步；两个 util 工具（`readflow3d-ver2.4a/2.5`）
  的旧组 namelist 补 `Iflag_Couple_Stop` + `Iflag_Couple_Restart`（后者属历史遗漏）。
- **不做**：可压块残差标准不新增参数（沿用其自身标准 `Res_Inner_Limit`）；SIMPLE 路径不导出残差
  ⇒ 该路径“残差达标”恒假（只可能由 `t_end`/轮数上限结束），避免误判收敛。

**验证（`mpirun -np 1`，全部 EXIT=0；临时目录 `/tmp/stg*`）**：

| 用例 | 配置要点 | 结果 |
|---|---|---|
| 码 19 `t_end` 停机 | `t_end=25, Kstep_Couple_Comp=10, Warm=99` | 第 3 轮气体段自动截断 10→6，`tt=26` ⇒ `stopped: physical time reached t_end` |
| 码 19 残差停机 | `t_end=1e6`、`AC_Tol=Res_Inner_Limit=1e30` | 第 1 轮 `all_res_ok=T` ⇒ `stopped: all global residuals below standard` |
| 码 11 case1 | `Solid_Max_Iter=2000` | `solid=7.09e-10` 实测登记；`Res_Inner_Limit=1e-10` 未达 ⇒ 继续跑，跑满轮数按旧提示退出 |
| 旧行为回归 | `Iflag_Couple_Stop=0` | 日志无 `stop test` 行，与旧版一致 |
| 代码 | `cd src && make` | EXIT=0 |

- 验证用 smokes 的日志副本：`cases/fluid_solid/run_case3_porous/run_stoptest_*.log`（t_end/残差/旧行为）、
  `cases/fluid_solid/run_case1_solid/run_stoptest_case1.log`。
- 手册：`chap07` 新增 §7.2「交错外层的停机判据」（含残差口径表 + 示例）、`91_appendix_params.tex`
  （组计数 12→13 + 新行）、`93_changelog.tex`（2026-09-26/1.2 一行）；`./build.sh` 重编译通过
  （`main.toc` 含 `subsection.7.2`，0 error）。`docs/control.ec-说明.md §9`、
  `docs/程序能力与算例考核总结.md` 同步。


## 0d. 最新（2026-09-25 第四轮）：`flow3d_node.dat` 补 PLOT3D 解(Q)/函数型文件头 ✅

**问题（用户提出）**：`flow3d_node.dat` 原本**只有裸数据记录**（`write(99) ((((U(m1,i,j,k)...`），
没有块维数与变量数 ⇒ 按 PLOT3D 解(Q)/**函数型文件**读法（Tecplot/ParaView）必然报格式错误，
读端只能用外部 `Mesh3d.x` 猜维数。

**改动（方案 A，2 个文件）**：
- `src/sub_IO.f90::output_flow_node`：数据记录前新增 **3 条头记录**
  ①`Total_block` ②`(bNi(m),bNj(m),bNk(m), m=1,Total_block)` ③`NVAR_NODE`（=6，新 `integer,parameter`）；
  仅 0 号进程写；注释块同步说明布局与「旧文件无头」。
  数据段写法**未改**（变量序 `d,u,v,w,T,Ts`、i 最快 → 与 PLOT3D Q 顺序一致），
  故**旧二进制输出 = 新文件去掉 3 条头记录**（逐字节一致）。
- `util/check_flow3d_node.py`：自动识别「有头 / 无头」两种布局，校验头记录
  （`nblock`、维数与 `Mesh3d.x` 一致、`nvars=6`）与每条记录长度 `6·ni·nj·nk`。

**验证（`mpirun -np 1`，EXIT=0）**：
| 算例 | 结果 |
|---|---|
| `cases/solid_1d`（1 块，51×9×9） | 头 = nblock 1、dims (51,9,9)、nvars 6（xxd 逐字节确认）；`STRUCTURE OK`；T=Ts=294–353 K |
| `cases/channel_ac_4blk`（4 块） | 7 条记录 = 3 头 + 4 数据，维数记录 48 B（4×3 整数）；`STRUCTURE OK` |
| 旧无头文件（删掉 3 条头） | checker 判为 `legacy (no header)`，`STRUCTURE OK`（兼容） |
| 等价性（HEAD 版旧二进制 A/B） | 两算例 `cmp` 与「新文件去头」**逐字节一致** ⇒ 数据段零变化 |

**文档**：`docs/重启与节点流场输出说明.md §3`（新增布局 + 修掉旧 `5*ni*nj*nk` 示例）、
`docs/程序能力与算例考核总结.md §2.8`、手册 `chap08`（表 + 节点小节加 lstlisting）、
`93_changelog.tex` 追加记录；`./build.sh` 重编译 **40 页、0 error/0 undefined**。

**注意**：集群副本需同步 `src/sub_IO.f90` 后 `make`；旧的 `flow3d_node.dat` 不用删（checker 兼容），
但要用 PLOT3D 读法直接加载则必须重跑一次生成带头文件。

## 0. 最新（2026-09-25 第三轮）：LaTeX 程序说明手册 v1.0 建立 ✅

**交付物**：`docs/程序说明/`（LaTeX 源码，`main.tex` + 12 个章节文件 + `build.sh` + `README.md`）
→ PDF `docs/OpenCFD-EC-1.16a-程序说明.pdf`（**40 页，xelatex 两遍编译，0 undefined**）。

**内容**：
- 第 1–2 章：块类型（0/1/2/3）、接口码 11–19、物理边界码（含 21/22）、面编号、**单位与无量纲化**、
  目录结构、上手流程；`control.ec` 7 组 namelist 规则、网格/bc3d/bc3d_interface/material/
  porous/solid_bc 各文件格式与示例。
- 第 3–6 章：**可压高速**（时间推进/空间格式/湍流/来流参考/输出节奏 + 选参建议）、
  **低速不可压**（变量约定、SIMPLE/SIMPLEC/AC-FV、`LS_*` 与 `AC_*` 参数表、**AC 标定经验表**、
  质量入口限制与等效速度换算）、**多孔**（DBF+Ergun/LTNE 公式、porous.inp/material.in、
  预算 `AC_Max_Iter × Porous_Chunk_Iter`、参数选取）、
  **固体导热**（FVM 非正交修正、SOR 参数、`Tw/Qw` 边界、CHT 界面温度公式与两种 CHT 模式、成本）。
- 第 7–8 章：跨区域耦合（同步耦合 vs 分段交错、调度表公式、判据、性能对比、
  接口配对/L 表机制、发汗壁、FAQ）+ 输出与重启（文件总览、`field_restart.dat` 布局 `iver=2`、
  **重启耦合状态与 `Iflag_Couple_Restart`**、`flow3d_node.dat`、MPI 约定）。
- 第 9–11 章：**三个 fluid_solid 算例详解**（拓扑/L 描述符、共同设置、case1/2/3 逐项参数与
  预期日志、时间预算、续算实测 2995→5990、调参清单）、构建运行与集群注意事项、验证/回归方法。
- 附录 A：**`control.ec` 全部参数与默认值总表**；附录 B：面编号/接口码/物理码/文件格式/单位/日志速查；
  **附录 C：更新记录（后续新功能追加于此）**。

**维护规则（已写入 `.clinerules` 与 `constraints.md §0`）**：新功能必须
① 更新对应章节 ② 同步附录 A（新参数） ③ 在附录 C 追加记录 ④ 重跑 `./build.sh` 确认无错误。

**状态**：**已提交并 push 到 `origin/main`**（两个提交：代码/功能 + 手册；详见 `worklog.md` 同日「收尾」条目）。

## 0a. 最新（2026-09-25 第二轮）：重启后按耦合状态决定「继续交错 / 切逐步强耦合」✅

**需求（用户提出）**：重启后不该把交错调度（暖机 + 每步减半）从头重播，而应基于上次状态判断；
**若已满足强耦合，就该用强耦合（逐步同时耦合）继续算**。

**改动（5 个源文件 + 3 个文档；`make` EXIT=0）**：
- `src/sub_modules.f90`：新开关 `Iflag_Couple_Restart`（默认 0）+ 耦合状态全局量
  （`Couple_State_Mode/Conv/Pair/Iter/nf/nk/Twmax/Tol/F found` + `Couple_Tw/qw/u/pw_save`）。
- `src/sub_restart.f90`：重启文件 `iver` 1→**2**，在所有块后追加 **coupling-state trailer**
  （mode/conv/pair/iter/nf,nk + `max|dT_w|/tol` + 界面量 `T_w,q_w,u,p_w`）；
  读侧 `iostat` 容错 + 0 号进程广播；新增 `set_couple_state` / `set_couple_tight_state`。
  **旧 `iver=1` 文件仍可读**（状态视为未知 ⇒ 旧行为）。
- `src/sub_multi_region.f90`：三个交错驱动（11/13、12、19）① 若 trailer 带同 pair、同形状的界面量
  ⇒ **直接恢复 `fp_Tw/fp_qw/fp_u/fp_pw`**（跨重启**零跳变**，不再从 `control.ec` 重播）；
  ② 每轮判据算完后 `set_couple_state(...)` 登记状态。`Iflag_Couple_Restart=-1` 时全部忽略（旧行为）。
- `src/multiphysmixsolver.f90`（更名前为 `opencfd_ec3d_v1.16a.f90`）：分派处按状态决定是否走交错驱动——
  `0`=自动（上次已满足 ⇒ 切强耦合）/`1`=强制强耦合/`2`=强制交错/`-1`=忽略；
  切强耦合时把 `Iflag_Couple_Scheme` 内部置 0、**自动限制固体每步预算**
  （`Solid_Max_Iter=min(现值,20)`、`Solid_Tol=max(现值,1e-6)`；只动扫掠上限与容差，物理/松弛因子不变）
  并回显实际值；主循环退出前**补写一次重启文件**（便于再次续算）。
- `src/sub_read_parameter.f90`：`Iflag_Couple_Restart` 走完 6 步（默认值 → `$couple_ec` → 旧组
  `$control_ec` → `bcast_para` 用 `Ipara(58)` → `output_para.out` 回显）。
- 文档：`control.ec.template`、`docs/control.ec-说明.md §8`、`docs/重启与节点流场输出说明.md §8`。

**验证（`mpirun -np 1`，全部 EXIT=0 / NaN=0）**：
| 测试 | 结果 |
|---|---|
| ct1 run1（case1 原设置等价性） | 12 轮、`max|dT_w|` 序列与改动前一致（0.0995 K 平台）、写重启 Kstep=2995 |
| ct1 run2（未收敛续算） | `mode=1 satisfied=0 pair=11 66x1` → **界面量已恢复（零跳变）** → 继续交错 12 轮 |
| ct2 run1（`Tol=1000` 促收敛） | it=2 判为已满足、写重启 Kstep=150（conv=1） |
| ct2 run2（自动切换） | `satisfied=1` → **切逐步强耦合**（`Solid_Max_Iter=20/Tol=1e-6`）→ 150→170 逐步推进 → 补写重启 |
| ct3a / ct3d | 状态已是强耦合 → 保持强耦合续跑（190→210）；`Iflag_Couple_Restart=1` 强制切强耦合 |
| ct3b / ct3c | `=2` / `=-1` 均**走交错驱动**（与旧行为一致） |
| 等价性 | `output_para.out` 仅多 1 行新回显；其余逐行一致 |

**注意 / 遗留**：
- case1 的 `max|dT_w|` 平台 ~0.0995 K > `Tol_Couple_Tw=1e-2` ⇒ **自动判定永不触发**；
  要用新行为需放宽 `Tol_Couple_Tw` 或直接 `Iflag_Couple_Restart=1`。
- 逐步强耦合即使 capped 后仍明显贵于交错（交错每外层轮只解一次固体）；长时间续算建议按需调
  `$solid_ec`。
- 本次改动**尚未提交**；集群需同步 `src/` 重新 `make`。

## 0a. 最新（2026-09-25）：修复集群 run1 启动崩溃 —— `control.ec` 同时连两个 unit ✅

**现象（集群 `run_case1_solid`，`./OPENCFDEC`）**：
`At line 976 of file sub_read_parameter.f90 (unit = 98)` →
`Fortran runtime error: File already opened in another unit`；
backtrace：`scan_control_ec_groups_ ← read_parameter_ec_ ← read_parameter_ ← MAIN__`。

**根因**：`src/sub_read_parameter.f90` 第 293 行 `open(99,file="control.ec")`（为读 7 组 namelist），
第 976 行扫描组头时又 `open(98,file="control.ec",status='old')` ⇒ **同一文件同时连到两个 unit**，
违反 Fortran 标准（F2018 12.5.6：一个文件同时只能有一个连接）。新版 libgfortran
（本地 gfortran 13.3.0，最小程序实测 `iostat=0`、不报错）放宽了该检查 ⇒
2026-09-23/24 的本地验证全部通过；集群的旧 libgfortran 按标准**致命报错**，
且该 `open` 未带 `iostat` ⇒ 运行时直接 `Error termination`（不是可捕获的 ios）。
**与算例、control.ec 内容、并行方式无关，纯代码可移植性缺陷**（2026-09-23 引入组扫描时带进）。

**改动（仅 `src/sub_read_parameter.f90`，1 file / +20 / −5）**：
`scan_control_ec_groups` 增加 `unit` 形参，内部 `rewind(unit)` + `read(unit,'(A)',iostat=ios)`，
**删除 `open(98,...)` 与 `close(98)`**（不 close：调用者还要用同一 unit 逐组
`rewind(99); read(99,nml=...)`）；调用点改为 `call scan_control_ec_groups(99, ...)`。
⇒ 全 `src/*.f90` 中 `control.ec` 只剩 1 处 open（第 293 行，unit 99），逻辑完全等价。

**验证**：`make` EXIT=0（改动文件零告警）；`grep` 静态审计确认唯一 open、unit 98 不再使用；
`cases/fluid_solid/run_case1_solid` 输入 + 新二进制在 `/tmp/val_case1` **整段跑完 12 轮**
（`max|dT_w|` 499.9 → 0.0995 K，与历史一致），写出 `field_restart.dat`(Kstep=2995) 与
`flow3d_node.dat`(6 变量)，**NaN=0**；`output_para.out` 与改动前二进制产物**逐字节一致**
（组回显仍为 `legacy=F freestream=T flow=T lowspeed=T ac=T solid=F porous=T couple=T`）。

**遗留**：① 两个 util 工具**有同类缺陷**（`util/readflow3d-ver2.5.f90` 99@707 + 96@791；
`util/readflow3d-ver2.4a.f90` 99@681 + 96@765），用户选择本轮不修 ⇒ 下次动它们时按同法修
（传 unit，不自己 open）；② **集群侧必须同步**修好的 `src/sub_read_parameter.f90` 到
`/work/home/lijunyang/sundong/PorousTest/code/` 并重新 `make`，否则 run1 仍报同样错误；
③ **已提交并 push**（与手册同一批；详见 `worklog.md` 同日「收尾」条目）。

**FAQ（同日实测）：case1 为什么只算到 Kstep≈2995？** —— 这是**设计终点**，不是卡死/发散：
`Iflag_Couple_Scheme=1` 时主程序调 `run_staggered_multiregion` 后**直接 `stop`，不进 `do while(tt<t_end)`**
⇒ `t_end=1e6` 无效；气体步数由调度表决定（前 `Niter_Couple_Warm` 轮满 `Kstep_Couple_Comp`，
之后每轮减半至 `Kstep_Couple_Min`），case1 `Warm=2/Outer=12`、`Comp/Min=1000/1(默认)`
⇒ 1000,1000,500,250,125,62,31,15,7,3,1,1 = **2995**；`max|dT_w|` 平台 ~0.0995 K > `Tol_Couple_Tw=1e-2`
⇒ 判据不触发，跑满 12 轮后正常退出（EXIT=0）。**续算**：直接再跑一次即可（`Iflag_restart=0`，
有 `field_restart.dat` 就恢复）—— 实测 run2：`Kstep 2995 → 3995 → … → 5990`（+2995 步，NaN=0）。


## 0. 最新状态（2026-09-23 第二轮）：重启文件 + 节点流场输出 ✅

**已完成并验证**（详见 `docs/重启与节点流场输出说明.md`）：
- **重启 `field_restart.dat`**（新文件 `src/sub_restart.f90`，开关 `$flow_ec:Iflag_restart/Kstep_restart`）：
  写 **U（含 LAP=4 鬼点缓冲区）+ Ts/Tsn + p + Kstep/tt**（头记录含 `Ma/Re/Lscale` 校验）；
  默认与 `Kstep_save` 同节奏；**启动时存在即自动续算**，不可用则回退原初始化，`Iflag_restart=1` 时 `stop 1`；
  读后刷新块间缓冲。交错驱动入口 `tt/Kstep` 归零已改为“有重启则恢复”，并在保存点与结束处写出。
- **节点流场 `flow3d_node.dat`**（`src/sub_IO.f90::output_flow_node`，开关 `Iflag_flow_node`）：
  8 格心平均→**网格节点**、**SI 单位**、**unformatted 逐块一条记录 `d,u,v,w,T`**、
  块序与 `Mesh3d.x` 一致（固体块 `d=u=v=w=0, T=Ts[K]`）⇒ 可与 `Mesh3d.x` 组合显示。
- 挂接：`sub_init.f90::Init_flow`（重启优先）、主循环、3 个交错驱动（补丁脚本 3+3+3 处）、
  `makefile`（+`sub_restart.o`）、`.gitignore`（忽略两个大文件）。
- 关键 bug 修复：`output_flow_node`/`read_restart` 原先在**非 0 号进程**上也遍历全局块并访问
  只在本进程有效的 `B_n(m)`（`-np>1` 会挂死/越界）→ 改为「主机走全局块、从机只走自有块」。
- 新增开关已加入 `bcast_para`（`Ipara(55/56/57)`）——**这是多进程必须的**（否则各 rank 取值不一致会死锁）。
- 验证：单块续算往返（Kstep 11→20，无 NaN）；节点文件维度与 `Mesh3d.x` 完全一致；
  **多进程 2 块** `-np 2` 写出并续算成功（SI 数值合理：Ma=3 → ρ≈0.1、u≈1700 m/s、T≈800–1174 K）；
  交错 case1 在 Kstep=50/100/150 写出两文件、无 NaN。

## 0c. 最新：三个 fluid_solid 耦合算例已迁移到新写法并验证重启 ✅（2026-09-24）

- `cases/fluid_solid/{run_case1_solid,run_case2_lowspeed,run_case3_porous}/control.ec`
  **已改写为 7 组写法**（只留非默认项；原文件备份为同目录 `control.ec.legacy`），
  并在 `$flow_ec` 里显式启用 `Iflag_restart=0, Kstep_restart=0, Iflag_flow_node=1`。
- **等价性**：三例「旧 `control.ec.legacy` vs 新 `control.ec`」短跑比对 `output_para.out`
  → 除 5 行回显外 **42 行参数逐项完全一致** ✅
- **重启可用**：三例 run1 均在交错推进中写出 `field_restart.dat`（case1 Kstep=2500/2995、
  case2 2000/2500、case3 500/1000）与 `flow3d_node.dat`；run2 自动续算并打印
  `restart from … Kstep=… ` + `continue staggered run at Kstep=…`，**NaN=0** ✅
- **节点文件**：三例均 2 记录、`(67,101,2)+(101,201,2)`、点数 13534/40602 与 `Mesh3d.x` 一致；
  固体块 `d=u=v=w=0`、T∈[300,799.9] K；可压块 SI 合理（case1 ρ=0.041–0.115 kg/m³、
  u=38.9–1706 m/s、T=793–2106 K ≈ Ma=3 滞止温 2240 K）✅
- 详见 `docs/重启与节点流场输出说明.md` §6。
- **补（同日）**：`flow3d_node.dat` 由 5 变量扩为 **6 变量 `d,u,v,w,T,Ts`** ——
  原先只填一个温度槽位，**多孔块骨架温度 `B%Ts` 没写进文件**（LTNE 温差不可见）。
  现在 `T`=流体温度（固体块=固体温度）、`Ts`=固体/骨架温度（有固相块为真实值，
  流体/低速块镜像 `T`）。实测 `porous_ltne_1d`：`Tf=300.7–1039 K`、`Ts=407.7–1174 K`、
  **`max|T−Ts|=157.7 K`**；校验工具 `util/check_flow3d_node.py`。

## 0b. 上一轮（2026-09-23）control.ec 拆分为 7 组 namelist ✅

**已完成并验证**（详见 `docs/control.ec-说明.md`、`control.ec.template`）：
- `control.ec` 由单个 `$control_ec`（128 变量）拆成 **7 组**：
  `$freestream_ec`(24) / `$flow_ec`(49) / `$lowspeed_ec`(21) / `$ac_ec`(15) /
  `$solid_ec`(4 新增) / `$porous_ec`(8) / `$couple_ec`(11)。
- **旧组 `$control_ec` 完全保留**（变量表 + `Solid_*`），**旧算例零改动**；
  新旧并存时新组覆盖旧组（支持逐组渐进迁移）。
- 读取逻辑：文本扫描组名（只有首个非空字符为 `$`/`&` 的行算组头）→ 读旧组 → 逐组可选读；
  **组存在但解析失败 ⇒ 报错 + `stop 1`**（不静默用默认值）；`output_para.out` 末尾回显各组来源。
- `$solid_ec` 的 4 个量原为 `solid_solver_one_block` 内硬编码（1.7/20000/5/1e-9），
  默认值与原值全等 ⇒ 行为不变。
- 同步修复：`util/readflow3d-ver2.4a/2.5` 现在也能读新格式，且补齐了原先缺失的
  `Lscale` / `LS_*` / `AC_*` / `Porous_*` / `couple` 变量（**历史遗留：这两个工具以前读不了现代算例**）。

**验证记录**：`make` EXIT=0；`cases/solid_1d` 与 `cases/channel_ac` 旧↔新二进制 A/B
`output_para.out` 逐字节一致（仅多 4 行回显）；新格式取值/覆盖/严格报错/模板 7 组全测通过；
两个 util 工具新旧格式都能跑完（EXIT=0）。

## 1. 本会话的关键决策（务必遵守）
1. **只改结构、不改物理**：`set_default_parameter` 与 `bcast_para` 不动，只改“从哪个组读”。
2. **兼容优先（第1档）**：`cases/**/control.ec` 一律不改动；不精简、不做每组最小模板、不做
   `$advanced_ec`（这些是后续可选档位，纯增量）。
3. **固体组必须非空**（Fortran 不允许空 namelist）⇒ 把 4 个硬编码量暴露为 `Solid_*`。
4. **严格失败**：组存在却解析失败必须 `stop 1`，不能静默回落默认值。
5. 新增参数走 6 步流程（见 `constraints.md §6`）。

## 2. 下一步（按优先级）
1. 三例 `cases/fluid_solid/*/README.md`（旧待办，仍未写）。
2. 可选：在三个 `cases/fluid_solid/*/control.ec` 里打开 `Iflag_flow_node=1`（本次**未改算例文件**，
   保持“既有算例零改动”；需要节点流场时在 `$flow_ec` 里加一行即可）。
3. 清理算例目录运行产物；提交（建议消息：`feat(restart): field_restart.dat + node-centred SI flow3d_node.dat`，
   若与拆分一起提交则 `refactor(control.ec)+feat(restart)`）。
4. 可选：`util/split_control_ec.py`（旧单组→7 组一键重排，值不变）。
5. ~~可选（改进）：把交错耦合的**界面量**（`fp_Tw/fp_u/fp_pw/fp_qw`）也写入重启文件，使
   `Iflag_Couple_Scheme=1` 的续算完全零跳变。~~ → **已于 2026-09-25 第二轮完成**
   （见 §0：ivre=2 trailer + `Iflag_Couple_Restart` 自动/强制强弱耦合切换）。
   剩余可做：`util/split_control_ec.py`；`cases/**/README.md`。

## 3. 历史上下文（仍生效的硬约束，详见 constraints.md）
- 交错耦合（`Iflag_Couple_Scheme=1`）**必须单进程**（码 12/19 无跨进程）。
- `cases/fluid_solid` 低速块必须 AC（`LS_Algorithm=3, AC_CFL=2`）；多孔预算
  `AC_Max_Iter×Porous_Chunk_Iter`；单轮耗时参考见 `constraints.md §5`。
- 面编号 1..6；接口码 11..19；`cases/fluid_solid` 拓扑 = `gas i+ (4) ↔ 对侧 j- (2)`。
- 未提交改动：`src/{sub_read_parameter,sub_modules,sub_multi_region,makefile,opencfd_ec3d_v1.16a}.f90`、
  `util/readflow3d-ver2.4a/2.5.f90`、`control.ec.template`、`docs/control.ec-说明.md`、`memory-bank/`。

## 0n. 2026-10-18 集群部署与 Slurm 复跑（已完成，后续按此执行）

- **集群副本已建立并验证**：`/work/home/lijunyang/sundong/MultiPhysMixSolver/`（`ssh wuzh02`）。
  同步一律 `bash util/sync_cluster.sh push|pull`；规则见 `.clinerules` 与 `constraints.md §8`。
- 集群 clean 重编（gfortran 7.3.1 + hpcx 2.11）：**0 error / 2 warning / 68 s / 38 个目标**；
  `bin/MultiPhysMixSolver` 3.15 MB，源码 md5 与本地一致。
- `cases/lidcavity_ac_3blk` 1e4 步 Slurm 复跑：np=1/2/3 墙钟 27/18/8 s；
  np2、np3 对 np1 = RMS(u) 5.77e-07、max|Δu| 1.16e-06；**集群 np1 与本地 `q1_10k` 逐位相同**。
- 用法：`cd cases/lidcavity_ac_3blk && ./prep.sh 10000`，然后 `cd np1 && sbatch -n 1 ../run.slurm`
  （np2/np3 同理，`-n` 必须等于 `partation.dat` 进程数）；分析先
  `module load apps/anaconda3/2023.09`（登录节点自带 python3 无 numpy），或拉回 VTK 本地分析。
- 脚本约定：`prep.sh` 的“步数”参数在运行目录内生成**独立** `control.ec`；
  **禁止** `sed ... > np*/control.ec`（会把软链目标截断成空文件）；`np*/` 到仓库根是三级
  （`../../../bin/`）。
- 待办不变：任务 5（L 形 90° 两块的 i↔j 界面）、任务 8（4/6 块 anchoring 灵敏性）、
  任务 7（强扩展，现可直接在集群做）；另 `git push`（本地领先 origin/main 15+ 提交）。

