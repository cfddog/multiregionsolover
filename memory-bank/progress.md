# Progress — 已完成 / 待办

> 更新：2026-09-24。**基线提交 `29fc24f`**（本次改动已全部提交；仅 `.vscode/settings.json` 未提交）。
> 更早的开发历史见 `docs/工作日志.md`、`docs/程序能力与算例考核总结.md`、`memory-bank/worklog.md`。

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

