# Active Context — 当前任务 / 最近决策 / 下一步

> 更新时间：2026-09-26（**最新 §0f：固体/低速/多孔区域残差监视输出**）。**HEAD = `fe7f86b`**
>（`29fc24f` 功能/重构主提交之后的记忆库/IO 提交均已 push）。
> 本次全部改动（control.ec 7 组拆分 + 重启/节点输出 + 三算例迁移 + 文档/记忆库）**已提交并推送**；
> 仅 `.vscode/settings.json`（主题设置）有意未提交。
> 上一版（2026-09-18，fluid_solid 三算例）已归档为 `memory-bank/archive_activeContext_2026-09-18.md`。
> **本轮（第四轮）改动尚未提交**：`src/sub_IO.f90`、`util/check_flow3d_node.py`、两份 docs、手册与记忆库。
> **第五轮（2026-09-26）改动同样尚未提交**：`Iflag_Couple_Stop` 停机判据（6 个源文件 +
> 2 个 util 工具 + 模板 + 4 份文档 + 记忆库）。
> **第六轮（2026-09-26）改动尚未提交**：区域残差监视输出（`src/sub_debug.f90`、
> `sub_multi_region.f90`、`sub_lowspeed(_ac).f90`、`sub_porous(_ac).f90` + 手册 4 处 + 总结文档 + 记忆库）。

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
- `src/opencfd_ec3d_v1.16a.f90`：分派处按状态决定是否走交错驱动——
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
