# Progress — 已完成 / 待办

> 更新：2026-10-19（**码 14（固--固）耦合修好，见 `activeContext.md` §0s**；
> 上一轮：2026-10-18 P3 收口；P0/P1/P3 已完成，P2/P4–P6 待做）。
> 更早的开发历史见 `docs/工作日志.md`、`docs/程序能力与算例考核总结.md`、`memory-bank/worklog.md`。

## ✅ 已完成（2026-10-19）：码 14（固--固）耦合 —— 对接描述符符号约定

- [x] 定位根因：`bc3d*.inp` 对接面两条切向范围**都取正** ⇒ `Convert_bc` 只能猜对应关系，
      `solid_2blk_chain`（两切向 21 / 3 节点）猜成 `k↔k`，界面量整错位，界面两侧差 ~38 K。
- [x] 约定固化：`bc3d.inp` **只给面类型**（同类对写 14，无配对行、无正负号；只有 `-1` 才带配对行）；
      `bc3d_interface.inp` 是块间连接**唯一**来源，对接面三对范围必须
      「一对待值=法向 / 一对全负=对接方向 / 一对全正=另一条切向边」。
- [x] `src/sub_convert_inp.f90`：`Convert_bc` 重写（`classify_pair` + `Convert_bc_badpair`），
      按类别配对得 `L1..L3` 并按 `abs(ke)-abs(kb)` 与同侧关系决定是否取反；不合约定直接中止。
- [x] `util/fix_iface_sign.py`（新，`--check`/`--fix`，列宽保持）；复核 `cases/*/bc3d*.inp`：
      本算例 changed=0/bad=0，另 3 例报「type-only 面（码 15）带符号」（无功能影响，待清理）。
- [x] 生成器：`solid_2blk_chain`（`with_link` 开关）+ `beavers_joseph`、`porous_fluid_phaseB`、
      `cylinder_re40_200`、`bl_cht`（对接方向取负）；`cases/solid_2blk_chain/run.slurm`（只能 np1）。
- [x] 集群 wuzh02 验证：`make -j4` EXIT=0；`prep.sh 2000` + `sbatch -n 1`（作业 46903753，1 m 03 s），
      `bc3d.inc` 两对接面为 `14` 且配对栏全 0，`bc3d_interface.inc` 得 `L = 1 2 3`；
      `analyze_solid_2blk_chain.py`：**PASS**（max 误差 3.7e−9 K，界面跳变 0.625000 K = 解析），
      结果存 `cases/solid_2blk_chain/analyze_np1.txt`。
- [x] 文档：`chap02_inputs.tex` 两小节重写（正负号规则 + 实例 + 工具用法）、`93_changelog.tex` 1.6.8，
      `./build.sh` 0 error。
- [ ] 待办：`high_low_fluid`/`high_low_fluid_800K`/`interface19_phaseA` 的码 15 仅类型面清理（待确认）；
- [ ] 待办：重生成并回归 `solid_4blk`、`fluid_solid`、`beavers_joseph`、`lidcavity_ac_3blk`、
      `channel_ac_2blk`、`bl_cht`（确认 `bc3d_interface.inc` 的 `L` 未变、结果不变）；
- [ ] 待办：码 14 的跨进程耦合（当前 np>1 直接中止）。


## ✅ 已完成（2026-10-18）：P3 并行固体 GS 推进（`Solid_Sync_Chunk`）

- [x] `sub_multi_region.f90`：扫掠体抽成 `solid_gs_sweep(mBlock, sweeps, resmax)`（纯局部，
      一次完整 GS 扫掠所有固体块 + 局部 `|dTs|` 判据）；`solid_solver_one_block`
      （np=1 / `Solid_Sync_Chunk<=0` 的旧路径）与新增 `solid_solve_all_blocks`
      （np>1 且 `Chunk>0`）共用它 ⇒ 两条路径物理完全一致，只有判据时机不同。
- [x] 并行路径：每个进程**同序**推进所有本地块，每 `Solid_Sync_Chunk` 次扫掠做一次
      `MPI_Allreduce(MPI_MAX)` + `MPI_Allreduce(MPI_LOR)`，**全局扫掠数/界面交换次数完全一致**
      （不再出现“某进程先满足 `|dTs|<Solid_Tol` 早退 ⇒ 阻塞消息错位挂死”）。
      终止条件 `(全局 max<=Tol .and. 所有进程都收敛) .or. sweeps>=Solid_Max_Iter`，
      并在每步打印一行 `solid sync (parallel): sweeps= ... (identical on every MPI rank)`。
- [x] 修 4 个 np>1 既有缺陷：
      ①`read_solid_bc` 用全局块号索引局部 `Block(:)`（越界写坏网格）；
      ②`check_remote_coupling_support` 用 Fortran 空语句 `continue` 当 `cycle`（np>1 一律拒界面）；
      ③`output_vtk` 的 MPI tag 用局部块号但收发两侧不一致（nc=2 时 `MPI_Recv` 处 glibc 堆损坏）
      ＋ 阻塞发送顺序不一致（交叉死锁）；
      ④`output_Ts` 用局部块号命名文件（两进程写同名 `Ts_block_1.dat`，`replace` 覆盖 ⇒
      误判“结果与进程数有关”，**这才是 P3 最初看到的 np 不一致的根因**）。
- [x] 新回归算例 `cases/solid_4blk/`（互不耦合 4 块，精确线性解；`gen_case.py`+`prep.sh`
      +`analyze_solid_4blk.py`+`README.md`）。
- [x] 验证（2 步，`Solid_Tol=1e-9`）：np=1 各块扫掠数 18100/5300/18700/5200，
      np=2 (`Chunk=100`) 全部 18700，无死锁；`flow3d.dat`/`Residual.dat`/`Step_mess.dat`
      与 np=1 **逐字节相同**；`Ts` 差 ≤3.9e-7 K（退出时机差异，= `Solid_Tol` 量级），
      blk3 逐字节相同；对精确解 ≤9.0e-7 K。反证：`Chunk=0` 跑 np=2 时各进程
      rank0=18100 / rank1=18700 各自早退，`Ts` 仍与 np=1 逐字节相同（旧路径未变）。
- [x] 集群 wuzh02 复验（`make -j4` EXIT=0；`prep.sh 2` + `sbatch -n 1/2`，作业
      46850363/46850364）：`analyze_solid_4blk.py` 数字与本地**完全一致**
      （3.9e-7 / 1.5e-7 / 0 / 1.7e-7 K，blk3 逐字节相同），且
      `flow3d.dat` md5 `2d0e4d22…`、`Ts_block_1.dat` md5 `77d54efb…`
      **与本地 gfortran 13 的结果逐位相同**（换机器/换编译器都不变）。
      现场参数 `np*/control.ec`（`Solid_Sync_Chunk=100`）+ `output_para.out` 已 pull 留存。
- [x] **文档一致性核查**（`Solid_Sync_Chunk` 七处落点全部对齐：默认值 / `$solid_ec` /
      旧组 `$control_ec` / `Ipara(61)` / `output_para.out` / 模板 / 附录 A + `chap06`），
      补掉两处漏网：`docs/control.ec-说明.md`（4→5 变量，顺带补漏记的 `AC_MPI_Sync`）、
      `util/readflow3d-ver2.4a/2.5.f90`（namelist 补 `Solid_Sync_Chunk`，否则旧组文件
      含该变量时工具报 `Cannot match namelist object name`；`gfortran -c` EXIT=0）。
      手册 `./build.sh` 重编译 EXIT=0。P1/P3 已 `git push` 到 GitHub `origin/main`。
- [ ] 遗留（与 P3 无关，2026-10-18 新发现）：**码 14 固--固局部耦合界面温度不连续**
      （np=1 收敛后两侧差 ~38 K）⇒ 待查 `couple_solid_solid_interfaces`。


## ✅ 已完成（2026-10-18）：P1 收口 —— 闭合域 AC 压力基准统一（低速 + 多孔）

- [x] `sub_modules.f90` 新增 `AC_Anchor_Block`；`init_ac_parallel_flags` 用一次 `MPI_MIN`
      取 `LS_Algorithm=3` 低速/多孔块中**全局块号最小**者（无 AC 块用哨兵 `ANCH_NONE`），
      启动打印 `pressure anchor block=`。
- [x] 新增 `ac_anchor_pressure`（`sub_lowspeed_ac.f90`）：有压力边界 ⇒ 不去均值；
      闭合域 ⇒ **只有锚定块**去均值。低速与多孔 AC 都改调它 ⇒ 修掉多孔侧
      “每个多孔块各减自己均值”造成的界面常数压差（闭合域 + 多孔块才会暴露）。
- [x] 回归（编译单元级，改前二进制在 `/tmp/p1old`）：`lidcavity_ac_3blk`（闭合）与
      `beavers_joseph/run_ac`（低速+多孔，有出口）np=1 5 步，5～6 个输出文件 **md5 相同**。
- [x] 闭合 3 块低速 `np=1 vs np=2` **1e4 步 PASS**（max|Δu|=1.15e-6，压力去常数 3.4e-7）。
- [x] 集群 `make -j4` EXIT=0；手册 `chap04` 新增“闭合域的压力基准”段 + `93_changelog` **1.6.6**。
- [ ] 待办：P2（`AC_Max_Iter` 跨进程统一步数语义）；P5 前处理新发现的问题——
      **闭合 3 个多孔块、np=2（1+2 分布）跑约 22 外步后停滞**（新旧二进制一致 ⇒ 既有问题，
      与 P1 无关；低速 1+2 分布正常）。


## ✅ 已完成（2026-10-18）：并行验证框架 `cases/parallel`（P3–P6 的回归底座）

- [x] 新增 `cases/parallel/`：`tests.conf`（测试清单）+ `run_para_verify.sh`（驱动）
  + `gen_control.py`（由算例 `control.ec` 生成本次控制参数，就地替换 + 缺键按前缀追加，
  同名 namelist 组冲突时**报错**）+ `analyze_para.py`（断言）+ `run.slurm` + `README.md`。
- [x] 比较面：`flow3d_node.dat`（d,u,v,w,T,Ts 全局块序，与进程数无关）+ VTK 单元中心压力
  （**去常数**：封闭域压力水平是 gauge，由全局 1 号块锚定）。
- [x] 两个基线：`t01_two_ac`（`channel_ac_2blk`，2 块，进口/出口）/ `t02_four_ac`
  （`lidcavity_ac_46blk` m4a，4 块，全封闭），都强制 `AC_Max_Iter=1`。
- [x] 集群实测（wuzh02，gfortran 7.3.1，10000 步，**全部 PASS**）：
  t01 np=2 max|Δu|=4.03e-6 / RMS 2.67e-6 / 压力 2.31e-6；
  t02 np=2,3,4 max|Δu|=1.37e-6 / RMS 7.47e-7 / 压力 1.22e-6，
  且 **np=2/3/4 之间逐位相同（max|Δ|=0）**。
- [x] 关键认识：np=1 是「逐伪步刷鬼点」（块 Gauss-Seidel），np>1 是「每外步刷一次」（Jacobi），
  **只收敛到同一个解**，瞬态必然不同 ⇒ 判据用收敛后的容差，不用逐位相同；
  2000 步的差距（t01 4.1e-5、t02 8.5e-4）只配做冒烟。
- [x] `util/sync_cluster.sh`：push 增加 `--exclude 'runs/'`（本地旧结果不得覆盖集群新结果，
  当天踩过一次：`rsync -a` 保留 mtime，把集群刚写的 `analyze_np*.txt` 盖回旧版）。



## ✅ 已完成（2026-10-18）：集群 `wuzh02` 连接打通
- [x] 私钥（`/mnt/c/Users/Administrator/.ssh/lijunyang_wuzh02_rsa`，与 Downloads 里
      `..._2026-11-12_14-15-23.txt` md5 相同）复制到 `~/.ssh/lijunyang_wuzh02_rsa`（600）。
- [x] 新建 `~/.ssh/config`（本机原先**没有**该文件）：`Host wuzh02` + HostName/Port 65091/
      User `lijunyang`/IdentityFile/`IdentitiesOnly yes`/`accept-new`。
- [x] 实测 `ssh -o BatchMode=yes wuzh02` 免密登录 `login03`，EXIT=0；`hostname/whoami/pwd` 均正确。
- [x] 环境核对：`mpif90` = **gfortran 7.3.1**（hpcx v2.11.0）；Slurm（`sbatch/squeue/sinfo`），
      分区 `wzhcnormal`/`wzacnormal`；`/work` 剩 2.3 TB；`rsync/scp/tar/md5sum` 齐全。
- [x] 确认旧副本 `~/sundong/PorousTest/code/` 仍是**单文件老布局**（≠ 现在的 `src/` 布局）。
- [x] 消掉集群 `perl: Setting locale failed` 噪声：`~/.ssh/config` 加 `SetEnv LANG=C/LC_ALL=C`。
- [x] `rsync -az` 上传探针成功、远端 md5 一致（探针已删）⇒ 同步通道就绪。
- [x] `cases/ssh/集群信息.md` GBK → UTF-8，并补进上述实测结论。
- [x] 记录作业提交方式（2026-10-18 用户提供）：帮助页
      <https://www.scnet.cn/help/docs/mainsite/hpc/cmd/slurm/>；模板在**集群**
      `~/sundong/run.slurm`，拷到目标算例目录后 `sbatch run.slurm` 提交
      （跑本仓库需把 `mixnsolver`→`../../bin/MultiPhysMixSolver`、`-n` 改成 np）。
- [x] 核实：本地根目录**没有** `run.slurm`（模板只在集群侧）；`~/sundong/MixSolver/MixSolver/`
      是另一个项目（`mixnsolver*`）；本仓库**尚无集群副本**。
- [ ] 待办：`git push`（`main` 领先 `origin/main` 15 个提交）。
- [ ] 待办：`rsync -az` 首次同步整仓到 `~/sundong/MultiPhysMixSolver/`，集群端 `make -j4` 验证。
- [ ] 待办：写 Slurm 脚本，集群上复跑 `lidcavity_ac_3blk` 的 np=1/2/3（老编译器 + 真网络）。

## ✅ 已完成（2026-10-17，收尾）：提交 4 个 commit + 提交前冒烟复验
- [x] `29eae92` fix(ac)：符号/is_high + `AC_MPI_Sync` + `p` halo + `IF_Debug=3` 自检
      + `init_ac_parallel_flags` + `read_solid_bc` 死锁（10 文件）。
- [x] `10ee6e2` cases(lidcavity_ac_3blk)：11 个输入/脚本/README + 目录内 `.gitignore`。
- [x] `6a804fa` docs(manual)：1.6.3（`chap04` 符号约定段 + `tab:ac3blk` + changelog + 2 份 PDF）。
- [x] `5094fd8` chore：删除根目录遗留 `control.ec`。
- [x] 提交前一致化：`ac_res_norm` 临时 `DBG2` 块删除；3 块算例 `AC_Max_Iter=1`、
      `t_end=Kstep_save=10000`；`control.ec.template` 闭腔→1 号块锚定（不再报错）。
- [x] 冒烟（新二进制，2000 外步）：np1/np2 无 NaN、u∈[−0.207,1]、`np2 vs np1 max|Δu|=8.7e-4`；
      `make -j2` EXIT=0。
- [ ] 待办：清临时目录（`np*/`、`s1`、`s42`、`m42`、`v1`、`v3`、`q*_10k`）；
- [ ] 待办：h 减半的界面一阶误差收敛性表（需 `p` 也写出 VTK 才能定量 `max|Δp|`）；
- [ ] 待办：`cases/ssh/`、`workplan.md` 是否入库待定（本轮有意未提交）。

## ✅ 已完成（2026-10-17，最新）：3 块 AC 方腔跨进程界面通量两处致命 bug 修复
- [x] 根因 1：`ac_face_flux_rus1_own` 两侧都做 `RAC += FLX*A`，而两侧面法向相同 ⇒
      Rusanov 耗散项在界面上相互抵消、界面 checkerboard 零阻尼、指数发散。改为
      低侧 `RAC -= FLX*A`、高侧 `RAC += FLX*A`。
- [x] 根因 2：`is_high` 判据写反（"ghost 索引更大"被当成高侧）⇒ 中央通量两侧同向反、
      界面变抽吸源（u≈1.68/1.83、界面整列错位 0.15）。改为按 own/ghost 索引关系判
      （`iG>iO`/`jG>jO`/`kG>kO` ⇒ 低侧），与块号、进程号无关。
- [x] 定位方法：错误解 np=1/2/3 **逐位相同** ⇒ 先排除竞争，走 np=1 隔离验证（含临时
      `DBGS` 打印），确认块 1 的 `iG=15>iO=14` 被判成高侧。
- [x] 进程数无关性：np=1/2/3（10^4 步）`max|Δu|=1.16e-6`、`max|Δv|=2.4e-7`、
      `max|Δp|=9.8e-7`，np=2 与 np=3 逐位相同。
- [x] 正确性：Re=100 方腔 `cases/lidcavity_ac_3blk` 3×10^4 步收敛（`res_q≈3e-14`，
      u∈[−0.233,0.926]）；vs Ghia RMS u 0.02095 / v 0.04365（与单块 40×40 已验证答案一致）。
- [x] 同网格单块对照：新增 `cases/lidcavity_ac_3blk/gen_single42.py`（42×42 单块），
      对照 `max|Δu|=9.3e-3`、rms 3.0e-3 —— 界面**刻意一阶**的预期误差，离界面 3 胞外逐位一致；
      单块 42×42 vs Ghia：u 0.02013 / v 0.04387。
- [x] 约定：`AC_Max_Iter=1` 为多块 AC 算例的**要求**（统一 np=1 与 np≥2 的 ghost 交换节奏）。
- [x] 文档：`chap04`「两侧各算路径的符号约定」段 + 表 `tab:ac3blk`；`93_changelog` **1.6.3**；
      `main.tex` 1.6.3；`./build.sh` **0 error / 0 undefined**。
- [x] 收尾：删除临时调试块、恢复 `ac_closed_dom = (.not. have_pbc) .and. Total_proc .gt. 1`、
      `ac_dbg` 回 0；`make` 干净重建通过。
- [ ] 待办：清理临时目录（`s1`/`s42`/`m42`/`v1`/`v3`/`q*_10k`）；可选：把界面一阶误差做成
      h 减半的收敛性表。

- [x] 顶层 `makefile`：`make` / `make -j4` / `make clean` / `make info`；模式规则
      `$(BUILDDIR)/%.o: $(SRCDIR)/%.f90`，`-Jbuild -Ibuild`，链接到 `bin/MultiPhysMixSolver`
      并 `ln -sf` 建兼容软链 `bin/opencfd-ec1.16a.out`。
- [x] `src/makefile` 改为转发入口（`cd src && make` 仍可用），`src/` 只剩 37 个 `*.f90` + makefile。
- [x] 新增 `util/gen_moddep.sh` + 入库的 `src/modules.dep`（36 条文件级依赖）⇒ `make -j4`
      **首次构建**即安全（此前必竞态失败）；不引入 `-cpp/-MMD`（会改变编译语义）。
- [x] `.gitignore` 加 `/build/`、`/bin/`、`*.d`；`.vscode/{tasks,launch}.json` 指向新路径；
      22 个 `cases/**` 的 README/脚本运行路径改为 `../../bin/MultiPhysMixSolver`。
- [x] 清理：`src/` 下 44 个中间文件（`*.o`/`*.mod`/`*.log`/可执行/旧软链）、根目录 2 个
      `*.mod`、`cases/**` 下 20+ 个旧 `opencfd-ec1.16a.out` 副本。
- [x] 验证：`make clean && make -j4` → EXIT=0（89 s）；`channel_ac_2blk` 5 步
      `flow3d.dat`/`field_restart.dat` 与重构前**逐字节相同**；P0 负例（np=2 码 12）仍 EXIT=2。
- [x] 手册 v1.6：`chap01`/`chap09`/`chap10`/`chap11`/`93_changelog`（1.6）+ `main.tex` 版本号；
      `./build.sh` 0 error；`docs/程序能力与算例考核总结.md` §7 同步。
- [ ] 待办：P2（AC 伪步统一步数语义；P1 的 `p` 跨进程 halo + AC 全局判据已于 2026-10-18 收口，
      见本文件顶部与 `activeContext.md` §0q）；集群副本同步（现在还要带上顶层 `makefile` 与
      `src/modules.dep`）。


## ✅ 已完成（2026-10-17）：分组提交 + ε 极限澄清 + 文档补全
- [x] 提交 `c974e81` feat(output)：`flow3d_node.dat` 恢复 **3 条 PLOT3D 头**
      （块数 / 各块维数 / 变量数=6）；`check/compare_flow3d_node.py` 兼容过渡的 **2 条**写法；
      `Iflag_restart<0` 时 `restarted` 未初始化 → 先置 `.false.`（冷启动确定性）。
- [x] 提交 `ebe7053` feat(solver)：`region_stat_mod` 区域残差监视、`Iflag_Couple_Stop` 停机判据、
      多孔 ε 极限（`[0,1]` 接受 / 阻力 `max(eps,1e-8)` + 速度冻结 / `k_eff` 界面导热），
      模板与 `util/readflow3d-2.4a/2.5` 同步。
- [x] 提交 `a29a6ee` cases(limits)：5 个短调度算例 + `run_limits.sh`（按脚本位置推导仓库根目录）+
      README（实测表、负结果分析、后续方案）。
- [x] 手册：`chap05`（ε 极限小节改为“回归快照 + 定量未通过”）、`chap09`（极限算例加注意段）、
      `chap11` 新增 §11.4「多孔 ε 极限的实测」表 + §11.5「2026-10-17 验证记录」、
      `93_changelog.tex` 追加 **1.4 / 1.4.1** 两行；`./build.sh` → **0 error / 0 undefined**。
- [x] `docs/程序能力与算例考核总结.md`：§2.4 极限处理改为实测口径、§5 新增 ε 极限“部分通过”条目、
      §6 新增第 15 条（负结果与下一步）。
- [x] 验证：`make -j4` EXIT=0（无 warning）；`solid_1d`/`channel_ac_4blk` 节点文件 `STRUCTURE OK`；
      `Iflag_restart=-1` 冷启动无 NaN；limits 5 例复查（`--blocks 1` 数字：17.64 / 28.02 / 96.34 / 498.6）。
- [ ] 待办：集群副本同步 + `make`（双 unit / `-freal-4-real-8`）；三例 `cases/fluid_solid/*/README.md`；
      ε 定量对照需同接口模型的新 deck（12↔12 或 19↔19）+ 更高外层/多孔预算。

## ✅ 已完成（2026-09-26 第六轮）：区域残差监视输出（固体/低速/多孔输出精简）
- [x] 需求：「固体、低速和多孔介质减少屏幕输出，输出区域最大 rms、平均 rms 以及最大 rms 的位置」；
      用户选定**方案甲**（不加参数，硬编码「每 100 次迭代 1 行 + 结束汇总」）。
- [x] `src/sub_debug.f90`：新模块 `region_stat_mod`（`RES_PRINT_STEP=100` + `region_res_stat`），
      统一打印一行 `Region res [名] block n (块类型) step=/final= k : max= … rms= … at (i,j,k)= … xyz=[…] m`。
- [x] 固体 GS（`sub_multi_region.f90`）、SIMPLE 低速/多孔（`sub_lowspeed.f90`/`sub_porous.f90`，
      `*_pressure_correction` 追加返回 `imax/jmax/kmax/res_p_rms`）、AC 低速/多孔
      （`sub_lowspeed_ac.f90::ac_res_stat`、`sub_porous_ac.f90`）、AC 能量与 LTNE（窗口统计）全部接入；
      删除所有逐迭代 `DBG/DEBUG/GS iter/SIMPLE iter/AC-* iter` 调试行。
- [x] 口径与可压块 `Res_max/Res_rms` 一致：`rms = sqrt(Σr²/N)`；`at (i,j,k)` + 胞心坐标 [m]。
- [x] 验证（`mpirun -np 1`，EXIT=0）：`solid_1d`、`channel_ac_4blk`、`porous_ltne_1d`（SIMPLE 与 `run_ac`）、
      `porous_plug_ac/eps06`、`fluid_solid` 码 11/12/19 短跑；`cd src && make` EXIT=0。
- [x] 文档：手册 `chap04/05/06` 新增/改写“残差监视输出/日志”小节、`92_appendix_codes.tex` 日志关键字表、
      `93_changelog.tex` 新增 **1.3**；`./build.sh` 重编译 **45 页 0 error**；
      `docs/程序能力与算例考核总结.md` 新增 §2.9。
- [ ] 待办：提交第四～六轮改动；同步集群副本（`src/sub_debug.f90` 等 5 个源文件）并 `make`。

## ✅ 已完成（2026-09-26 第五轮）：交错外层停机判据 `Iflag_Couple_Stop`
- [x] 需求：交错耦合需要「以 `t_end` 为结束标志，或以全局所有残差 < 残差标准为结束标志」。
- [x] `src/sub_modules.f90`：新 `Iflag_Couple_Stop`（默认 0 旧行为 / 1 新判据）+
      `Couple_Res_AC_q/AC_m/Solid`（每轮残差登记量）。
- [x] `src/sub_multi_region.f90`：新 `couple_stop_test`（可压块 `B%Res` RMS + `MPI_ALLREDUCE`；
      按存在块型比对 `Res_Inner_Limit`/`AC_Tol`/`Solid_Tol`）；三个驱动（11/13、12、19）轮末调用、
      轮首清零、气体段按剩余时间截断、退出提示区分 `t_end`/残差；固体 GS 登记末次 `res`。
- [x] `src/sub_lowspeed_ac.f90`、`src/sub_porous_ac.f90`：AC `res_q/res_m` 登记（本轮 `max`）。
- [x] `src/sub_read_parameter.f90` 6 步流程（`Ipara(59)` + `$couple_ec` + 旧组 + 回显）；
      `control.ec.template`、`util/readflow3d-ver2.4a/2.5.f90`（旧组 namelist 补 2 个变量）同步。
- [x] 验证（`mpirun -np 1`，EXIT=0）：码 19 `t_end=25` ⇒ 第 3 轮截断 10→6、`tt=26` 停机；
      放宽容差 ⇒ 第 1 轮残差判据停机；码 11 固体残差 `7.1e-10` 正确登记；`Iflag_Couple_Stop=0`
      输出与旧版一致。
- [x] 文档：手册 `chap07` §7.2 + 附录 A（12→13）+ `93_changelog.tex`；`./build.sh` 通过；
      `docs/control.ec-说明.md §9`、`docs/程序能力与算例考核总结.md`。
- [ ] 待办：提交本轮（含第四轮未提交内容）；同步集群副本并 `make`。

## ✅ 已完成（2026-09-25 第四轮）：`flow3d_node.dat` 增加 PLOT3D 解(Q)/函数型文件头
- [x] 起因：用户报「`flow3d_node.dat` 函数型格式错误，应补充维度和变量数」——
      文件原先**无任何头信息**，只有逐块 `6·ni·nj·nk` 裸数据，PLOT3D Q/函数型读法无法识别。
- [x] `src/sub_IO.f90::output_flow_node`：新增 3 条头记录
      `Total_block` / `(bNi,bNj,bNk, m=1,Total_block)` / `NVAR_NODE(=6)`（新 `integer,parameter`），
      仅 0 号进程写；数据段写法不变（变量序与 PLOT3D Q 一致）。
- [x] `util/check_flow3d_node.py`：识别「有头/无头」两种布局，校验头记录 + 记录长度。
- [x] 验证：`solid_1d`（1 块）与 `channel_ac_4blk`（4 块）跑通、`STRUCTURE OK`；
      旧无头文件仍 `STRUCTURE OK`；与 HEAD 版旧二进制 A/B **数据段逐字节一致**。
- [x] 文档：`docs/重启与节点流场输出说明.md §3`、`docs/程序能力与算例考核总结.md §2.8`、
      手册 `chap08` + `93_changelog.tex`；`./build.sh` 重编译 40 页 0 error/0 undefined。
- [ ] 待办：把改动同步到集群副本（`src/sub_IO.f90` + `make`）；`cases/fluid_solid/*/README.md` 仍未写。

## ✅ 已完成（2026-09-25 第三轮）：LaTeX 程序说明手册 v1.0
- [x] 新建 `docs/程序说明/`（`main.tex` + 12 章节 + `build.sh` + `README.md`），
      xelatex（ctex/Fandol 字体）两遍编译通过，**0 undefined**，PDF 40 页：
      `docs/OpenCFD-EC-1.16a-程序说明.pdf`。
- [x] 覆盖：块类型/接口码/物理边界码/面编号/单位与无量纲化；`control.ec` 7 组规则；
      网格/bc3d/bc3d_interface/material/porous/solid_bc 文件格式；
      可压高速、低速（SIMPLE/SIMPLEC/AC-FV，含 AC 标定经验表）、多孔（DBF+LTNE+Ergun）、
      固体导热（含 CHT 界面公式）四类求解器与**参数选取建议**；
      跨区域耦合（同步 vs 交错、调度表、判据、成本对比、FAQ）；输出与重启（`iver=2` + 耦合状态）；
      **三个 fluid_solid 算例详解**；构建运行与集群注意事项；验证/回归方法与已知问题。
- [x] 附录 A `control.ec` 全部参数与默认值总表；附录 B 编号/格式/单位/日志速查；
      **附录 C 更新记录（后续新功能追加于此）**。
- [x] `.clinerules`、`constraints.md §0`、`projectbrief.md`、`README.md` 写入
      **“新功能必须同步更新手册 + 附录 A/C + 重编译”** 的维护规则；
      `.gitignore` 忽略 LaTeX 中间文件。

## ✅ 已完成（2026-09-25 第二轮）：重启后按耦合状态续算（继续交错 / 切逐步强耦合）
- [x] `field_restart.dat` `iver` 1→2：追加 coupling-state trailer（mode/conv/pair/iter/nf,nk +
      `max|dT_w|/tol` + 界面量 `T_w,q_w,u,p_w`）；旧 `iver=1` 文件仍可读；
      读侧 iostat 容错 + 0 号进程广播；新增 `set_couple_state` / `set_couple_tight_state`。
- [x] 三个交错驱动：界面量可从 restart 恢复（**跨重启零跳变**）；每轮登记状态。
- [x] 新开关 `Iflag_Couple_Restart`（`$couple_ec`，默认 0 自动 / 1 强制强耦合 / 2 强制交错 / -1 忽略）；
      6 步流程齐全（`Ipara(58)`）。
- [x] 主程序分派：自动判定已满足 ⇒ 内部 `Iflag_Couple_Scheme=0` 走**逐步强耦合**（`t_end` 生效），
      并自动 cap 固体每步预算（`Solid_Max_Iter≤20`、`Solid_Tol≥1e-6`）、退出前补写重启文件。
- [x] 验证（`mpirun -np 1`）：ct1（12 轮等价 + 续算零跳变）、ct2（收敛→自动切强耦合→逐步推进）、
      ct3a–d（四种开关取值）全部 **EXIT=0 / NaN=0**；`output_para.out` 仅多 1 行新回显。
- [x] 文档：`control.ec.template`、`docs/control.ec-说明.md §8`、`docs/重启与节点流场输出说明.md §8`。

## ✅ 已完成（2026-09-25）：修复集群 run1 启动崩溃（`control.ec` 同一文件双 unit）
- [x] 定位：`src/sub_read_parameter.f90:976` `open(98,file="control.ec",status='old')` 与第 293 行
      `open(99,file="control.ec")` 冲突 —— 同一文件同时连两个 unit，违反 F2018 §12.5.6。
- [x] 解释“本地通过 / 集群崩”的差异：本地 gfortran 13.3.0 实测同一文件开两个 unit `iostat=0`
      （放宽了检查）；旧 libgfortran 按标准**致命报错**，且该 open 无 `iostat` ⇒ 直接终止。
- [x] 修复：`scan_control_ec_groups(unit, ...)` 传 unit 99，删 `open(98)/close(98)`
      （**1 文件 / +20 / −5**，逻辑等价；`control.ec` 全仓库只剩 1 处 open）。
- [x] 验证：`make` EXIT=0、改动文件零告警；`/tmp/val_case1`（case1 输入副本）**12 轮跑完、NaN=0**、
      写出 `field_restart.dat`/`flow3d_node.dat`；`output_para.out` 与旧二进制参考**逐字节一致**。

## ✅ 已完成（2026-09-23 会话）：control.ec 拆分为 7 组 namelist
- [x] 盘点：`$control_ec` 共 **128 个变量**、**全部有代码默认值**；实测算例中 70~75% 的赋值行与
      默认值相同（case1 121行/112条→75条冗余；channel_ac 93/91→65；lidcavity 89/87→66）。
- [x] `src/sub_modules.f90`：新增 `Solid_GS_Omega/Solid_Max_Iter/Solid_Min_Iter/Solid_Tol`。
- [x] `src/sub_read_parameter.f90`：1 组 → **7 组 namelist**（freestream 24 / flow 49 /
      lowspeed 21 / ac 15 / solid 4 / porous 8 / couple 11，合计 128+4，变量一个不漏）；
      旧组 `$control_ec` 完整保留（+`Solid_*`）。
- [x] 读取逻辑：**扫描组头**（首个非空字符为 `$`/`&`；`!` 后为注释）→ 旧组 → 7 组可选读，
      **新组覆盖旧组**；组存在但解析失败 ⇒ 报错 + `stop 1`；`output_para.out` 末尾回显各组来源。
- [x] 新增 `scan_control_ec_groups` / `str_to_lower` / `check_nml_ios` 三个辅助例程。
- [x] `src/sub_multi_region.f90`：`solid_solver_one_block` 的 4 个硬编码常量改为读 `$solid_ec`
      （默认值与原值全等 ⇒ 物理行为不变）。
- [x] 顺手修掉一个**既有潜在 bug**：`read_parameter_ec` 里 `a0/d0/mu0/mu1` 在非叶轮机械模式下
      未初始化就被写进 `output_para.out`（garbage，逐次运行不同）→ 现已显式清零。
- [x] `util/readflow3d-ver2.5.f90`、`util/readflow3d-ver2.4a.f90`：同样双读；并补齐
      `Lscale` + `LS_*`/`AC_*`/`Porous_*`/`couple`/`Solid_*` 声明
      （**历史遗留：这两个工具原先连 `Lscale` 都没有，读现代算例必报 iostat=5010**）。
- [x] 新增 `control.ec.template`（7 组带中文注释、标注默认值与“可整组省略”）与
      `docs/control.ec-说明.md`（分组表/兼容规则/示例/迁移建议/验证记录）。
- [x] 记忆库同步：`projectbrief.md §4`、`constraints.md §6`（新增参数 6 步流程）、
      `activeContext.md`（旧版归档为 `archive_activeContext_2026-09-18.md`）。
- [x] **验证**：`make` EXIT=0；`cases/solid_1d` 与 `cases/channel_ac` 旧↔新二进制 A/B
      `output_para.out` **逐字节一致**（仅多 4 行回显）；新格式取值/覆盖/严格报错/模板 7 组
      全测通过；`Solid_Max_Iter=30` 生效（日志 GS iterations 31=30+1）；两个 util 工具
      新旧格式均 EXIT=0。

## ✅ 已完成（2026-09-23 第二轮）：重启文件 + 节点流场输出
- [x] 新增 `src/sub_restart.f90`：`write_restart()` / `read_restart(ok)` / `restart_step_output(nMesh)`。
      重启文件含 **U(1:NVAR,含 LAP=4 鬼点缓冲) + Ts + Tsn + p + Kstep/tt**（头记录 + 每块头），
      校验 `nblock/NVAR/LAP/块1维数/Ma/Lscale`；不匹配 → 回退原初始化（`Iflag_restart=1` 时 `stop 1`）。
- [x] `src/sub_IO.f90`：新增 `output_flow_node()` + `flow_node_one_block()`：
      8 格心平均→**网格节点**、**SI 单位**、**unformatted 逐块一条记录 `d,u,v,w,T`**、
      块序与 `Mesh3d.x` 一致（固体块 `d=u=v=w=0, T=Ts[K]`）。
- [x] 控制量：`Iflag_restart`(0 自动/1 强制/-1 关)、`Kstep_restart`(<=0 用 `Kstep_save`)、
      `Iflag_flow_node`；进 `$flow_ec` 与旧组 `$control_ec`；`bcast_para` 用 `Ipara(55/56/57)`；
      `output_para.out` 回显；`control.ec.template` 同步。
- [x] 挂接：`sub_init.f90::Init_flow`（重启优先于 `Iflag_init` 路径）、主循环每步调用
      `restart_step_output(1)`；三个交错驱动「`tt/Kstep` 归零→有重启则恢复」+ 保存点调用 +
      运行结束写出（补丁脚本 3+3+3 处，含 `contains` 位置修正）。
- [x] `src/makefile` 增加 `sub_restart.f90`；`.gitignore` 忽略 `cases/**/field_restart.dat`、
      `cases/**/flow3d_node.dat`。
- [x] **修 bug**：`output_flow_node`/`read_restart` 原先在非 0 号进程上也遍历全局块、访问只在本进程
      有效的 `B_n(m)` ⇒ `-np>1` 越界/挂死；改为「主机走全局块、从机只走自有块」（与 `output_flow` 同构）。
- [x] **验证**：单块续算往返（Kstep 11→20，无 NaN）；节点文件维度与 `Mesh3d.x` 一致
      （4131/13534/40602 点）；**`-np 2` 2 块**（固体+可压）写出并续算成功，SI 数值合理
      （Ma=3：ρ≈0.10–0.13 kg/m³、u≈788–1701 m/s、T≈800–1174 K）；`-np 2` 空闲 rank 不死锁；
      交错 case1 在 Kstep=50/100/150 写出两文件、无 NaN。
- [x] 文档：新增 `docs/重启与节点流场输出说明.md`；`docs/程序能力与算例考核总结.md` 增 §2.8；
      `docs/control.ec-说明.md` 更新 `$flow_ec` 变量表。

## ✅ 已完成（2026-09-24）：fluid_solid 三算例迁移 + 重启实测
- [x] `cases/fluid_solid/{run_case1_solid,run_case2_lowspeed,run_case3_porous}/control.ec`
      改写为 7 组新写法（只留非默认项；原文件备份 `control.ec.legacy`）；
      `$flow_ec` 显式启用 `Iflag_restart=0, Kstep_restart=0, Iflag_flow_node=1`。
- [x] 等价性 A/B（旧 vs 新，`output_para.out` 过滤回显后逐行比对）：三例 **42 行参数完全一致**。
- [x] 重启实测（单进程）：三例 run1 写出 `field_restart.dat` + `flow3d_node.dat`；
      run2 自动续算（`restart from … Kstep=…` + `continue staggered run at Kstep=…`），**NaN=0**。
- [x] `flow3d_node.dat` 结构校验：2 记录、点数 13534/40602 与 `Mesh3d.x` 一致；
      固体块 T∈[300,799.9] K；可压块 SI 合理（Ma=3 峰温 2106 K ≈ 滞止温 2240 K）。

## ✅ 已完成（2026-09-24 补）：节点文件补上固体/骨架温度 `Ts`（5→6 变量）
- [x] 起因：`flow_node_one_block` 只往一个温度槽位填值（固体填 `Ts`、多孔填 `Tf`），
      **多孔块的骨架温度 `B%Ts` 根本没写进 `flow3d_node.dat`**（LTNE 温差不可见）。
- [x] 修复：`flow3d_node.dat` 每条记录由 `d,u,v,w,T` 扩为 **`d,u,v,w,T,Ts`**：
      `T`=流体温度（固体块为固体温度，沿用 `flow3d.vtk` 约定）、
      `Ts`=固体/骨架温度（固体与多孔为真实值；流体/低速块无固相 → `Ts` 镜像 `T`）。
      改动：`src/sub_IO.f90::flow_node_one_block`（新增 `Tsc`，4 个分支分别赋值）与
      `output_flow_node`（`Num_data/allocate/write` 用 6）。
- [x] 新增校验工具 `util/check_flow3d_node.py <算例目录>`（比对 `Mesh3d.x` 块数/维数、
      记录长度 `6·ni·nj·nk`，打印各块 `d/|u|/T/Ts` 范围与 `max|T-Ts|`）。
- [x] 验证：`case1`（固体+可压）固体块 `T=Ts=300–800 K`、可压块 `T=Ts` 镜像；
      **`porous_ltne_1d`（`q=2e6 W/m²`, `hv=2e6 W/(m³K)`）→ `Tf=300.7–1039 K`、
      `Ts=407.7–1174 K`、`max|T−Ts|=157.7 K`**（修复前该温差在文件里不存在）；
      `case3`（多孔+可压）多孔块 `T`/`Ts` 分列；三例均 `STRUCTURE OK`、`NaN=0`。

## ⏳ 待办（本次未做，按优先级）
- [ ] **集群同步（立即）**：把修好的 `src/sub_read_parameter.f90` 拷到
      `/work/home/lijunyang/sundong/PorousTest/code/` 并重新 `make` —— 否则集群 run1 仍会报
      `File already opened in another unit`；建议整份 `src/` 同步，避免版本漂移。
- [ ] 同类隐患（2026-09-25 按用户要求**未修**）：`util/readflow3d-ver2.5.f90`（99@707 + 96@791）、
      `util/readflow3d-ver2.4a.f90`（99@681 + 96@765）——同一 `control.ec` 连两个 unit，
      旧编译器下同样崩；修法与 `scan_control_ec_groups` 相同（传 unit、不自己 open）。
- [x] ~~提交本次修复（建议消息：`fix(io): 不再二次打开 control.ec`）~~ → **已提交并 push**
      （与 LaTeX 手册同批，两个提交；见 `worklog.md` 2026-09-25「收尾」条目）。
- [ ] **（上一任务）重启文件 `field_restart.dat`**：写流场 `U`(含 LAP ghost 缓冲) + `Ts/Tsn` + `p`
      + `Kstep/tt`；`Kstep_save` 节奏定期保存；启动时存在即自动续算（`Iflag_restart`）。
      新开关放入 `$flow_ec`。交错驱动入口的 `tt/Kstep=0` 需改为“续算时恢复”。
- [ ] **（上一任务）节点/SI 流场输出**：`flow3d.dat` 8 格心平均 → 网格节点、SI 单位、unformatted、
      与 `Mesh3d.x` 组合可显示（`output_flow_node`，开关 `Iflag_flow_node` 放 `$flow_ec`）。
- [ ] 三例 `cases/fluid_solid/*/README.md`（几何/边界表、G→V、限制、复现命令、单进程要求）。
- [ ] 文档同步：`docs/程序能力与算例考核总结.md` 补“control.ec 7 组 + 双读兼容”。
- [ ] 清理算例目录运行产物；提交（消息建议：`refactor(control.ec): split into 7 namelists with legacy dual-read`）。
- [ ] 可选：`util/split_control_ec.py`（旧单组→7 组一键重排，值不变）；
      精简现有算例文件（只能删“与默认值相同”的行 + `output_para.out` 逐项 diff）。
- [ ] 已知未实现（历史遗留）：接口 17/18；12/19 跨进程耦合；固体/骨架 Ts 与 LS 温度的重启读回。

## 📌 长期待办（来自 docs/程序能力与算例考核总结.md §9，摘录）
- 多孔速度入口“流量锁定”（~28% 通量缺口）；低速块网格病态（high_low_fluid/flatplate）；
  SIMPLE 压力求解改用 BiCGStab/ILU；耦合外迭代加 Aitken/Anderson；合并 VTK 改分区输出等。

## 📁 历史记录：已完成（2026-09-18 会话，fluid_solid 三算例）
- [x] 摸清 `cases/fluid_solid` 几何：Block1 67×101×2 板（j=1 为对接面、j=101 冷端/入口）、
      Block2 101×201×2 通道（i=101 为对接面、i=1 天花板、j=1 入口 / j=201 出口、k 对称），
      界面共形 67 节点，网格单位 mm（`Lscale=0.001`）。
- [x] 用运行时诊断确认界面连接描述符：`blk2 L=(2,-1,3)`（gas i↔对侧 j、切向 j↔i 反向）。
- [x] 新增 `dbg_dump_interfaces`（`IF_Debug=1`）+ 主程序调用（`src/opencfd_ec3d_v1.16a.f90`）。
- [x] 扩展码 11：`run_staggered_fluid_solid` 面检查放宽为 `solid j- vs comp i-/i+`。
- [x] 扩展码 12：`couple_highlow_fluid_face` 新增 `(i+,j-)` 分支；`run_staggered_highlow` 支持并改界面量索引。
- [x] 扩展码 19：`couple_compressible_porous_blowing_face` 新增 `(i+,j-)` 分支；
      `run_staggered_fluid_porous` 重构为 `pair42` 双分支（壁 ghost / q_w / 多孔 j- ghost / Ts 热流 / 输出）。
- [x] 新增发汗壁 `wall_bound_blowing`（pair42 下气体界面壁等温 + 法向喷射）。
- [x] 编译通过（`src/opencfd-ec1.16a.out`，2026-09-18 13:44，EXIT=0）。
- [x] 新建三个算例目录 `run_case1_solid/ run_case2_lowspeed/ run_case3_porous/`
      （Mesh3d.x、bc3d.inp、bc3d_interface.inp、material.in、control.ec；case1 加 solid_bc.inp，
      case3 加 porous.inp，二进制软链）。
- [x] 交付版 `control.ec` 参数已用短时限运行验证（三例 `flow/gas chunk steps=1000` 正确打印）。
- [x] 三例短调度试跑：EXIT=0、无 NaN（case1 CHT `max|dT_w|` 499.9→0.094 K；case3 `max|q_w|` 5.3e6→3.9e6 W/m²）。
- [x] 写 `docs/进度_fluid_solid三算例.md`（中断点记录）并在 `docs/工作日志.md` 末尾加指针条目。
- [x] 建 `.clinerules` 与 `memory-bank/`（本记忆库）。
- [x] **第二轮现场测试（交付版设置整段跑）**：三例均能启动/推进/无 NaN ——
      case1 **12 轮全跑完 99 s**（调度 1000→1000→500→250→…，`max|dT_w|` 499.9→0.0996 K 平台）；
      case2 界面度量逐轮下降（`max|dT|` 321→9.2 K）但 420 s 限时退出；
      case3 首轮气体段完成、多孔 AC `res_q`=5.6e-5，但原多孔预算在 437 s 内跑不完首轮。
- [x] **case3 多孔段预算压缩并复验**：交付 `control.ec` 改为
      `AC_Max_Iter=5000, AC_Print=1000, AC_CFL=20, AC_beta=100, Porous_Chunk_Iter=3`；
      复验（`run_check_new.log`）：`pair42=T`、1000 步气体段、多孔 AC 3×5001 迭代
      `res_q` 9.7e-5→6.4e-5、无 NaN。（原 `20000×10` 太慢，每轮 ≈65 min。）
- [x] **case2 的 LS-AC 标定结论**（保留原稳定设置 `LS_Algorithm=3, AC_CFL=2, AC_beta=10`）：
      `AC_CFL=20, AC_beta=100` 使 `res_q` 降到 3.3e-3 但**第 3 轮 NaN**；
      `LS_Algorithm=1`(SIMPLE) **第 1 步即发散** ⇒ 必须用 AC + `AC_CFL=2`。

## 📁 历史记录：当时待办（2026-09-18；多数仍未做，见上方 09-23 待办）
- [ ] 三个算例目录的 `README.md`（几何/边界表、G→V 换算、限制、复现命令、单进程要求）。
- [ ] ~~case2 的 AC 标定~~ → 已得结论（见上，保留 `AC_CFL=2`；`res_q` 平台 ~1.4e-2 属已知限制）；
      如要更紧收敛需改进 LS-AC 求解器本身（非设置问题）。
- [ ] ~~case3 外层收敛加速~~ → 已完成预算压缩；如仍偏慢可再降 `Porous_Chunk_Iter` 或 `AC_Max_Iter`。
- [ ] case1 若希望外层早停：`Tol_Couple_Tw` 由 1e-2 放宽到 ~0.2 K（当前平台 ~0.0996 K，属固体未达稳态）。
- [ ] case3 关注 `coolant exit |v|_max`（冷启动 10.3→15.4 m/s，入口 1.699 m/s）是否随外层收敛而回落。
- [ ] 清理算例目录运行产物（保留 `run_trial.log`、`iface_*.dat`）；注意 case1 曾有含空格怪文件 `'force-debug-  0'`。
- [ ] 文档同步：`docs/程序能力与算例考核总结.md` §1 接口表补“11/12/19 支持 gas i+ ↔ 对侧 j- 拓扑”。
- [ ] 提交（消息建议：`feat(staggered): support gas i+ <-> partner j- topology for codes 11/12/19 in cases/fluid_solid`）。
- [ ] 可选修复：`sub_lowspeed.f90::lowspeed_inlet_velocity` / `sub_porous.f90::porous_inlet_velocity`
      按**面法向**赋速度分量并明确面积量纲（使 `LS_Inlet_Type=2` 可直接按 kg/(m²·s) 输入）。
- [ ] 已知未实现（历史遗留）：接口 17/18；12/19 的**跨进程**耦合；固体/骨架 Ts 与 LS 温度的重启读回。


## 2026-10-17 · 程序更名 MultiPhysMixSolver

- [x] `git mv src/opencfd_ec3d_v1.16a.f90 src/multiphysmixsolver.f90`；文件头注释 + 启动横幅。
- [x] `src/makefile`：`EXE = MultiPhysMixSolver`，链接后建兼容软链 `opencfd-ec1.16a.out`；
      `clean` 同时删除两者。**串行 `make`**（`-j4` 会因缺模块依赖竞态失败）。
- [x] `sub_read_parameter.f90` 的 `output_para.out` 版本行；`.gitignore` 加 `/src/MultiPhysMixSolver`。
- [x] 22 个 `cases/**` 的 README / `gen_case.py` / `run_*.sh`；`.vscode/{tasks,launch}.json`。
- [x] 手册 v1.5（`main/chap01/chap09/chap10/93_changelog` + `build.sh` 产物名）；PDF 49 页 0 error；
      旧 PDF 文件名同步更新保留兼容。
- [x] 文档源码引用：`docs/程序能力与算例考核总结.md`、`docs/重启与节点流场输出说明.md`、
      `cases/beavers_joseph/README.md`、`cases/porous_fluid_phaseB/staggered_tw/README.md`。
- [x] 逐位一致回归：更名前/后二进制在 `channel_ac_2blk`（2 块 + AC，5 步）全部输出逐字节相同
      （`output_para.out` 仅版本行不同）。
- [ ] 集群副本同步 `/work/home/lijunyang/sundong/PorousTest/code/` 后 `make` 复验。
- [ ] 数据文件头字符串（`flow3d.dat`/`Ts_block_*.dat`/`bc3d.inc`）**保持旧名**；若将来要改，
      需同时给两个 util 脚本加“双读兼容”，单独一次提交。

## 并行改造待办（P0→P6，均未开始）

> 2026-10-18：`cases/parallel/` 已提供**np 无关性验证框架**（见上方「已完成」与
> `cases/parallel/README.md`）。P3–P6 每一项落地后，都必须先在本框架里加对应测试
> （固体块 / 码 12 / 码 19 / 码 14、18）并 PASS，再做算例层验证。

- [x] **P0 跨进程安全网**（2026-10-17 完成）：`check_remote_coupling_support`
      （`src/sub_multi_region.f90` 末尾）在 `Init_flow` 后遍历界面，统计远端且无 MPI 版本的
      组合（12/14/17/18/19），`MPI_Allreduce` 后 `MPI_Abort(2)` 并打印清单；np=1 直接返回。
      验证：np=1 `channel_ac_2blk` 5 步 `flow3d.dat`/`field_restart.dat` 逐字节不变；
      `fluid_solid` deck 改 `material.in=0 2` + np=2 → EXIT=2 且打印 `12(...)`；
      原 `fluid_solid`（码 16）np=2 无假阳性。手册 chap07 + 1.5.1 记录。
## 完成（2026-10-17，最新）：AC `p` halo 统一 + 缓冲区交换自检

- [x] 删除 `Umessage_send_mpi` 里的融合 `p` 拷贝；`p` 改由
      `p_send_mpi`/`p_recv_mpi` 单一路径覆盖（同进程直写 / 跨进程 Bsend），
      过滤条件 `LOWSPEED|POROUS × LOWSPEED|POROUS` 与 `U` 侧对齐。
- [x] 低速/多孔 AC 驱动条件统一为 `AC_MPI_Sync==1 .or. Total_proc==1`
      （多孔 np=1 此前完全不交换）。
- [x] 新增 `IF_Debug=3` → `check_buffer_exchange`：哨兵值逐 ghost 胞严格比对，
      非破坏性；跨进程用独立 tag 把“邻居 ghost 应有内容”发给对侧比较。
- [x] 实测三种分布全部 PASS：np=1/2 块（320 胞）、np=2/各 1 块（320 胞）、
      np=2/各 2 块（960 胞），`max|dU|=max|dp|=0`，三例 EXIT=0。
- [x] 文档：`chap04` §ac-halo/§bufchk、附录 A（`AC_MPI_Sync`、`IF_Debug=3`）、
      `93_changelog` 1.6.1；`build.sh` 0 error。

## 完成（2026-10-17）：手册 §4.7「多块离散一致性」（含曲面/非正交适用性）

- [x] `chap04` 新增 `\subsection{多块离散一致性：界面模板与通量}`（§4.7，`sec:ac-mblock`）：
      表 17（对流/粘性/角点的多块依赖）＋ 对流"三件事"（`ii=2..nx-1` 防双计、`FSH` owner
      单值通量、界面**刻意**一阶 Rusanov 的稳定理由）＋ 跨进程 `id=0` 回退与
      `AC_MPI_Sync=1` 前提 ＋ 界面粘性"等值反号"表与三条前提（共形/镜像 ghost/ghost 新鲜）
      ＋ 角点判据（无混合导数；应力张量只需 ghost-edge；真正需要角点的五类格式）
      ＋ 曲面/非正交适用性 ＋ 表 18（`cylinder_re40_half` 实测）。
- [x] 结论：**多块不需要额外粘性处理**（两层差分自然等值反号）；**O 型/极坐标网格不失效**
      （笛卡尔分量 + 常黏度 + 法向投影 ⇒ 仍二阶；实测 Cd +0.3% vs Rogers&Kwak）。
- [x] `93_changelog` **1.6.2**、`main.tex` 版本 1.6.2；`./build.sh` **0 error / 0 undefined**
      （第二遍消解前向引用，PDF 无 `??`）。
- [x] 记忆库：`activeContext` §0k、`worklog` 本条。
- [ ] 遗留（与本条相关）：若将来上 9 点/紧致/顶点型格式，**必须先扩 ghost 角点交换**
      （当前只填面 slab）。

## 待办（最新排序）

- [ ] `ac_zero_mean_p` 改为**全局**均值（`MPI_Allreduce`）：无压力 BC 的多块 AC
      目前每个块各自减本地均值，界面 `p` 会有块间偏置。
- [ ] `AC_MPI_Sync=0 & np>1`：伪步内不刷新 halo（仅启动告警），可考虑自动升级为 1。
- [ ] 接口通量注册表 `FSH` 跨进程时的回退（两侧各算一阶 Rusanov）导致 np 相关精度，
      P1/P2 的后续项。
- [x] **P3** 固体扫描数在 `np>1` 时统一（`Solid_Sync_Chunk`）—— **已完成 2026-10-18**，
      见本文件 §「✅ 已完成（2026-10-18）：P3 并行固体 GS 推进」与 `activeContext.md` §0r。
- [ ] **P4** 码 12（低速--固体）跨进程：`hl_iface_mpi`，几何克隆 `ls_por_iface_mpi`，两侧都调用。
- [ ] **P5** 码 19（可压--多孔）跨进程：`cpb_iface_mpi`，非对称 2 标量打包（气体 `p_phys`、多孔 `T_w`）。
- [ ] **P6** 码 14/18 跨进程（`ss_iface_mpi` 或复用 `Ts` 缓冲）；手册新增《并行能力矩阵》章；
      4 块/4 进程混合算例 `t07_mix4`。
- [ ] `cases/parallel/`：`mkpart.py` 生成 `partation.dat`、`run_para_verify.sh` 跑 np=1 vs np=2
      并断言 `max|d|=0`；测试 t01–t07 随各级加入。
- [x] **P3** 固体扫描数在 `np>1` 时统一（`Solid_Sync_Chunk`）—— **已完成 2026-10-18**，
      见本文件 §「✅ 已完成（2026-10-18）：P3 并行固体 GS 推进」与 `activeContext.md` §0r。
- [ ] **P4** 码 12（低速--固体）跨进程：`hl_iface_mpi`，几何克隆 `ls_por_iface_mpi`，两侧都调用。
- [ ] **P5** 码 19（可压--多孔）跨进程：`cpb_iface_mpi`，非对称 2 标量打包（气体 `p_phys`、多孔 `T_w`）。
- [ ] **P6** 码 14/18 跨进程（`ss_iface_mpi` 或复用 `Ts` 缓冲）；手册新增《并行能力矩阵》章；
      4 块/4 进程混合算例 `t07_mix4`。
- [ ] `cases/parallel/`：`mkpart.py` 生成 `partation.dat`、`run_para_verify.sh` 跑 np=1 vs np=2
      并断言 `max|d|=0`；测试 t01–t07 随各级加入。


## 已完成（2026-10-18：任务 5 旋转界面 + 任务 8 锚定回归 + 左手块保护）

- [x] 修 `cases/channel_ac_2blk_rot90/gen_case.py`：`p` 变体 **左手→右手**（块 2 `k→-z`，
      `L=(2,+1,-3)`；bc 里邻居 k 区间反序 `-2 -1`）；`n` 变体 `L=(2,-1,+3)` 不变。
- [x] 新增启动保护 `src/sub_geometry.f90::check_mesh_quality_onemesh`：最小单元体积 `<0`
      ⇒ 打印块号/最小体积/`(i,j,k)` + `MPI_Abort(3)`；`negative_vol` 时不打“OK”。本地
      故意造左手网格实测 exit=3、`Min Volume=-1.25e-3`。
- [x] 集群 5 个作业（46827434/36/37/39/40，4×10⁴ 步）全部完成：
      `p` vs 未旋转参考**逐位相同**（3.0e-16）、`n` vs 参考 4.0e-4（差值在被反转块内部）、
      np=1 vs np=2 各 2.4e-10；墙钟 4m19/4m22/6m01/2m51/2m49。
- [x] `cases/channel_ac_2blk_rot90/README.md`（新建）、`cases/lidcavity_ac_46blk/README.md`
      回归表（进程数/锚点/Ghia/界面跳变）；`analyze_rot90.py` 文档串同步 `L` 约定。
- [x] 文档：`chap04_lowspeed.tex` 新增“索引旋转（i↔j）的块与左右手约定”、
      `chap10_build_run.tex` 启动流程加体积符号检查、`93_changelog.tex` **1.6.5**；
      `./build.sh` + 第 3 遍 = 64 页、0 error、0 未定义引用、3 个 PDF 同步。
- [x] 集群/本地同步（push 全仓、pull 只取记录）；集群 `make -j4` 通过（gfortran 7.3.1）。
- [ ] 待提交：本轮改动 **尚未 `git commit`**（`src/sub_geometry.f90`、`cases/*rot90*`、
      `cases/lidcavity_ac_46blk/README.md`、docs 3 份、memory-bank 3 份）。
- [ ] 经验入库：`sbatch` 提交前**不要** `module load anaconda`（mpirun 混用 ⇒
      `NULL communicator`）；`make` 同理。

## 已完成（2026-10-18：集群部署与 Slurm 复跑）

- [x] 集群副本 `/work/home/lijunyang/sundong/MultiPhysMixSolver/`（rsync 1943 文件 / 50 MB，
      排除了 `.git/build/bin`、VTK/重启/网格大文件与 `np*/`）。
- [x] 集群 gfortran 7.3.1 clean 重编：**0 error / 2 warning / 68 s**（38 目标）。
- [x] Slurm 复跑 `cases/lidcavity_ac_3blk` np=1/2/3（1e4 步）：墙钟 27/18/8 s；
      np2/np3 对 np1 = RMS(u) 5.77e-07 / max|Δu| 1.16e-06；集群 np1 与本地 `q1_10k` 逐位相同。
- [x] 同步规则入 `.clinerules` / `constraints.md §8`；工具 `util/sync_cluster.sh`
      （push/pull）；`prep.sh` 支持步数参数；新增算例 `run.slurm`；
      手册 chap10 新增“集群作业与同步”小节 + changelog 1.6.4，`./build.sh` 0 error。
- [ ] 待办：任务 5（L 形 90° 两块的 i↔j 界面）、任务 8（4/6 块 anchoring 灵敏性）、
      任务 7（强扩展，集群上做）；`git push`。

