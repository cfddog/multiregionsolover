# porous_fluid_phaseB / staggered_tw — 分段交错耦合 (interface 19) 试跑

本目录是在 `../porous_fluid_phaseB`（接口码 19：可压通道块 1 ↔ 多孔壁块 2）
网格基础上，测试新实现的**分段交错(segregated/staggered) 外层迭代**策略：

1. **高速流段**：只推进可压块连续 `Kstep_Couple_Comp` 步；其“对接(code 19)段”
   底壁按**等温壁**（温度 `T_w(x)`，初值 `Twall_Couple_Init=300 K`）施加，
   其余 code-2 壁段仍为绝热壁（本几何整条底壁都是 code 19，无绝热段）。
   多孔块此阶段不推进。
2. 由可压解提取界面**热流 q_w(x)**（>0 进入壁面，W/m²，按 Sutherland 粘性 +
   Pr 计算）与**壁压 p_w(x)**。
3. **多孔段**：只推进多孔块（含底部冷却剂供给边界，`LS_*`/入口类型按 control.ec）；
   热端（y=0，code 19）流体出口 ghost 压力 = `p_w`，骨架 Ts 施加热流边界
   `Ts_g = Ts_i + q_w·dx/k_s_eff`，跑 `Porous_Chunk_Iter` 次 SIMPLE/LTNE 内迭代。
4. 由多孔解返回**界面温度** `T_w(x)=(Ts_i+Ts_g)/2` 作为下一高速流段的等温壁温度；
   直到 `max|ΔT_w| < Tol_Couple_Tw`（本轮未达，仅验证调度/数据流）。

## 新增控制量（control.ec namelist，缺省均保持原有逐时间步耦合不变）
```
Iflag_Couple_Scheme = 1   ! 0=原逐时间步同时耦合; 1=分段交错(仅接口码19 FLUID<->POROUS)
Kstep_Couple_Comp   = 1000! 每个高速流段推进的可压步数（用户原意 1000）
Niter_Couple_Outer  = 60  ! 外层交错最大迭代数
Porous_Chunk_Iter   = 10  ! 每个多孔段调用 porous_solver_one_block 的次数
Twall_Couple_Init   = 300.d0  ! 首次高速流段的等温壁温 [K]
Tol_Couple_Tw       = 2.d-2   ! 外层收敛判据 max|ΔT_w| [K]
```

## 代码改动（本功能）
- `src/sub_modules.f90` / `src/sub_read_parameter.f90`：上述新 namelist 控制量。
- `src/opencfd_ec3d_v1.16a.f90`：`Iflag_Couple_Scheme=1` 时主程序改走新的
  `run_staggered_fluid_porous`（否则原 `while(tt<t_end)` 时间循环不变）。
- `src/sub_multi_region.f90`：新增 `run_staggered_fluid_porous`（含
  `fill_gas_wall_ghost`、`set_porous_Ts_flux`、`T_nd_from_U` 内部例程）；
  目前支持一个共形配对：可压 **j- 面 ↔ 多孔 j+ 面**（与
  `porous_fluid_phaseB` 拓扑一致）。

## 试跑记录（本目录）
- `If_viscous=0`（无粘试调度）：无崩溃；q_w 由“等温镜像 ghost 温差”人为产生且
  恒定≈33.28 kW/m²（`wall_bound_with_Tw` 对过低 T_w 会把 ghost T 钳制到 0.5·T_int，
  因此无粘时 T_w 变化不改变 q_w），T_w 每轮单调上升
  （300→302.2→304.3→306.2→308.1 K），**不能收敛**——物理上也确实无热交换。
- `If_viscous=1`（粘性，10 步/段 短试）：q_w **随 T_w 变化**
  （34.38→35.17→35.98 kW/m²，T_w 302→304→306 K），调度与热流反馈均正常；
  说明高速段每段必须推进足够步数（如 1000）让边界层/温度剖面稳定后 q_w 才收敛。
- 默认路径回归：`check_default/` 中 `Iflag_Couple_Scheme=0`、t_end=5 短跑正常。

## 注意 / 待用户处理
1. **必须有粘性**才有物理热流：请设 `If_viscous=1`（层流/湍流按需）；且当前
   121×81 均匀通道网格对 Ma3/Re≈5.6e6/m 层流边界层远不够（首层 dy≈0.6 mm，
   需要 ~1e-5 m 量级加密）。
2. **等温/绝热分段**：若希望“前段绝热平板 + 后段多孔发汗段”，需把可压块 j- 底壁
   拆成 code-2（绝热）与 code-19（耦合）两个子面（几何上对应多孔块只铺在后段），
   本驱动会自动把 code-19 段作等温壁、code-2 段留作绝热壁。
3. 多孔热端热流目前全部进入**骨架 Ts**（经 hv 与冷却流体换热），流体相出口为
   零梯度+压力出口——若需表面孔隙分流/直接加热冷却剂需按表面孔隙率扩展。
4. 外循环通常需对 `T_w`/`q_w` 欠松弛，或“暖机大步长 → 减半小步长”精修（见下节），
   并保证每段气体推进步数足够让温度剖面稳定。
5. 实现为单块对、共形、可压 j- / 多孔 j+；多 MPI 进程时暂按单进程使用。


## “暖机 + 减间隔”两阶段用法（生成初始状态后精修）
用同一进程即可，无需文件重启：
```
Iflag_Couple_Scheme = 1
Kstep_Couple_Comp   = 1000   ! 阶段1：大步长（每轮跑满，把高速流场/壁温大致稳定）
Niter_Couple_Warm   = 3      ! 前3轮都用 1000 步暖机
Kstep_Couple_Min    = 20     ! 暖机后每轮气体步数自动减半，直到 20 步下限
Niter_Couple_Outer  = 12     ! 总外层轮数（暖机3 + 递减精修 ~7-9）
Tol_Couple_Tw       = 1.d-2  ! 最终收敛判据
```
驱动实现：第 1~`Niter_Couple_Warm` 轮气体段固定为 `Kstep_Couple_Comp` 步；
之后每轮把步数减半（下限 `Kstep_Couple_Min`），日志会打印
`outer iter k : gas chunk steps = n`。这样“先粗稳、再细修”，最终落到
`max|ΔT_w|<Tol_Couple_Tw`。

说明：若想把某轮结果作为**下一个进程/算例**的初场继续跑（如换 `Iflag_Couple_Scheme=0`
逐时间步耦合），代码会输出 `flow3d.dat`/`Ts_block_*.dat`/VTK，但当前
`Init_flow` 的重启读入路径（`Iflag_init>0` 读 `flow3d.dat`）只恢复可压区，
并把多孔区 p、Ts 重置为 `LS_P_out`/`Porous_T_ref`（见 `sub_init.f90`），
多孔骨架温场不会被带进新进程——所以推荐在**同一运行内**用上面的两阶段参数，
或后续再给多孔 Ts 增加读回功能。

