# Active Context — 当前任务 / 最近决策 / 下一步

> 更新时间：2026-09-25（新增 §0a：集群启动崩溃修复）。**HEAD = `fe7f86b`**
>（`29fc24f` 功能/重构主提交之后的记忆库/IO 提交均已 push）。
> 本次全部改动（control.ec 7 组拆分 + 重启/节点输出 + 三算例迁移 + 文档/记忆库）**已提交并推送**；
> 仅 `.vscode/settings.json`（主题设置）有意未提交。
> 上一版（2026-09-18，fluid_solid 三算例）已归档为 `memory-bank/archive_activeContext_2026-09-18.md`。

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

**遗留**：本次改动**尚未提交**。

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
③ 本次改动尚未提交。

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
