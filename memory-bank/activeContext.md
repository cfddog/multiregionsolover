# Active Context — 当前任务 / 最近决策 / 下一步

> 更新时间：2026-09-24。**基线提交 `29fc24f`**（功能与重构主提交；已 push 到 `origin/main`，
> `git log` 的最新提交为记忆库/文档更新）。
> 本次全部改动（control.ec 7 组拆分 + 重启/节点输出 + 三算例迁移 + 文档/记忆库）**已提交并推送**；
> 仅 `.vscode/settings.json`（主题设置）有意未提交。
> 上一版（2026-09-18，fluid_solid 三算例）已归档为 `memory-bank/archive_activeContext_2026-09-18.md`。

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
5. 可选（改进）：把交错耦合的**界面量**（`fp_Tw/fp_u/fp_pw/fp_qw`）也写入重启文件，
   使 `Iflag_Couple_Scheme=1` 的续算完全零跳变（目前界面量从 `control.ec` 初值重播、由外层迭代再收敛）。

## 3. 历史上下文（仍生效的硬约束，详见 constraints.md）
- 交错耦合（`Iflag_Couple_Scheme=1`）**必须单进程**（码 12/19 无跨进程）。
- `cases/fluid_solid` 低速块必须 AC（`LS_Algorithm=3, AC_CFL=2`）；多孔预算
  `AC_Max_Iter×Porous_Chunk_Iter`；单轮耗时参考见 `constraints.md §5`。
- 面编号 1..6；接口码 11..19；`cases/fluid_solid` 拓扑 = `gas i+ (4) ↔ 对侧 j- (2)`。
- 未提交改动：`src/{sub_read_parameter,sub_modules,sub_multi_region,makefile,opencfd_ec3d_v1.16a}.f90`、
  `util/readflow3d-ver2.4a/2.5.f90`、`control.ec.template`、`docs/control.ec-说明.md`、`memory-bank/`。
