# Active Context — 当前任务 / 最近决策 / 下一步

> 更新时间：2026-09-18（会话中断点）
> 基线提交 `adaa43b`（分支 `main`）；**本次改动全部未提交**，在工作树中。

## 1. 当前任务
`cases/fluid_solid` 上新建**三个分段交错耦合算例**（用户将投集群测试）：

| 算例目录 | 界面码 | Block1 | Block2 |
|---|---|---|---|
| `run_case1_solid/` | 11 流体↔固体 | 不锈钢固体，j=max 等温 300 K | Ma=3 空气, Re=5000, T∞=800K, 壁面绝热 |
| `run_case2_lowspeed/` | 12 流体↔低速 | 低速空气，j=max 质量入口 G=2 kg/(m²·s), 300K | 同上 |
| `run_case3_porous/` | 19 可压↔多孔 | 不锈钢多孔 ε=0.30，j=max 冷却剂注入 G=2, 300K | 同上 |

交错调度三例相同：`Iflag_Couple_Scheme=1, Kstep_Couple_Comp=1000, Niter_Couple_Warm=2,
Kstep_Couple_Min=1, Niter_Couple_Outer=12`（= 1000,1000,500,250,125,62,31,15,7,3,1,1 步）。

**状态（2026-09-18 第二轮现场测试已完成 ✅）**：算例文件已建好；交付版 `control.ec` 参数已验证；
三例均能按交付设置启动/推进/**无 NaN** —— case1 **12 轮全跑完 99 s**；case2/case3 在 420/437 s
限时内正常推进（限时是我方测试上限，非崩溃）。case3 多孔段预算已压缩并复验通过。
**尚未写 README、未清理产物、未提交。**

## 2. 本会话的关键决策（务必遵守）
1. **不改网格**：`fluid_solid` 网格界面 = `Block1 j=1 ↔ Block2 i=101`，
   而既有耦合例程硬编码 `gas j- ↔ 对侧 j+`。为保“Block1 j=1 为对接面”的约定，
   选择**扩展代码**（而非镜像网格/重新编号），新增对 `gas i+ (4) ↔ 对侧 j- (2)` 拓扑的支持。
2. **物理公式复用**：本网格两块法向物理上都是 y ⇒ `U(2)=切向x、U(3)=法向y` 语义不变，
   只改**下标/鬼单元索引**，旧 `(2,5)` 分支代码逐字保留（向后兼容）。
3. **质量入口用等效速度**：低速/多孔为常密度（`LS_rho`），`LS_Inlet_Type=2` 通路按
   “i 面 + 网格单位面积”实现，用于 j 面会方向/量纲错 ⇒ 交付为
   `LS_V_in = Porous_V_in = G/LS_rho = 2/1.177 = 1.699 m/s`（+y），`LS_rho=1.177`（300K/1atm 空气）。
4. **Block2 的 i=1 面 4(farfield) → 2(绝热壁)**，按用户“壁面为绝热壁面”。
5. **新增发汗壁** `wall_bound_blowing`：pair42 下气体界面壁等温 + 法向喷射
   （`v_w = Bp%U(3,gi,jc_l,k)`），让冷却剂喷注真正作用于气体边界层（仅 `LS_Inlet_Type==1` 且 `v_w>0`）。
6. **新增界面诊断** `dbg_dump_interfaces`（`IF_Debug=1` 打印 face/face1/nb1/ranges/L1..L3）。
7. **case3 多孔段预算压缩（第二轮，已改交付 `control.ec`）**：原
   `AC_Max_Iter=20000 × Porous_Chunk_Iter=10` 使首轮多孔段 ≈65 min（437 s 跑不完）；
   改为 `AC_Max_Iter=5000, AC_Print=1000, AC_CFL=20, AC_beta=100, Porous_Chunk_Iter=3`
   （≈170 s/轮、`res_q`=6.4e-5，与原相当），已复验无 NaN。
8. **case2 必须用 AC 且 `AC_CFL=2`**（第二轮实测）：`AC_CFL=20, AC_beta=100` 虽让 `res_q`
   降到 3.3e-3 但**第 3 轮出现 NaN**；`LS_Algorithm=1`(SIMPLE) **第 1 步即发散**
   （`res_p`→6e31, `u`→-1e6）。⇒ 保留 `LS_Algorithm=3, AC_CFL=2, AC_beta=10, AC_Max_Iter=5000`；
   `res_q` 平台 ~1.4e-2 属该配置已知限制（界面度量逐轮下降，结果可用）。
9. **case1 判据不改**：`max|dT_w|` 平台 ~0.0996 K > `Tol_Couple_Tw=1e-2`，因 50 mm 厚固体
   +300 K 冷端未达稳态（`T_w` 每轮缓降 ~0.1 K），属物理瞬态；若要外层早停可放宽到 ~0.2 K。

## 3. 下一步（按优先级）
1. 写三个目录的 `README.md`（几何/边界表、G→V 换算、限制、复现命令、单进程要求）——**首要项**。
2. 清理算例目录运行产物（每目录约 15 MB：`flow3d.vtk` 11 MB + `flow3d.dat` 3.3 MB + 若干 log/inc），
   保留 `run_smoke.log`、`run_check_*.log`、`iface_*.dat` 作证据。
3. 文档同步：`docs/程序能力与算例考核总结.md` §1 接口表补“11/12/19 已支持 gas i+ ↔ 对侧 j- 拓扑”。
4. 提交（2 个源文件 + 3 个算例目录 + `docs/进度_fluid_solid三算例.md` + `.clinerules`/`memory-bank`）。
5. 可选（非设置问题）：改进 LS-AC 求解器本身以压低 case2 的 `res_q` 平台；
   修 `lowspeed_inlet_velocity`/`porous_inlet_velocity` 的**面法向分量 + 面积量纲**，
   使 `LS_Inlet_Type=2` 可直接按 kg/(m²·s) 输入。
6. case3 后续关注：`coolant exit |v|_max`（冷启动 10.3→15.4 m/s，入口 1.699 m/s）是否随外层收敛回落。

## 4. 上下文提醒
- 详细中断点记录（含命令、日志摘录、备份路径）见 `docs/进度_fluid_solid三算例.md`。
- 未提交改动：`src/sub_multi_region.f90`（+651/−344）、`src/opencfd_ec3d_v1.16a.f90`（+2）。
- 回滚：`git checkout -- src/ && rm -f src/opencfd-ec1.16a.out && make`。
