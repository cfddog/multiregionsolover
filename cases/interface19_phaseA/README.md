# interface19_phaseA — 可压↔多孔接口（码 19）阶段 A：框架与结构试跑

## 目标
实现并验证 `BLOCK_FLUID(可压) ↔ BLOCK_POROUS(多孔)` 接口（码 19）的**基本框架**：
- 接口分派：`couple_fluid_solid_interfaces` 增加 FLUID↔POROUS 分支，复用码 12 的
  跨单位换算耦合（可压无量纲 ↔ 多孔 SI，`U(1)=ρ[SI]`、`U(2..4)=u[m/s]`、`U(5)=T[K]`、
  `p[Pa]`，RHO_REF 由 Re/μ_Sutherland/Ma/T_inf 反推，速度钳制，压力锚定判断）。
- 本阶段**不含冷却剂吹出、不含骨架界面热平衡**：多孔块骨架 Ts 在界面绝热（G=0）。

## 代码改动（本阶段，均向后兼容）
1. `src/sub_multi_region.f90`：分派器把 `BLOCK_FLUID↔BLOCK_POROUS` 配对交给
   `couple_highlow_fluid_face`（多孔块与低速块同构存储 + 共用 LS_* 参数），
   反向（POROUS→FLUID 侧）跳过（每对从可压侧处理一次）。
2. `src/sub_boundary.f90`：`Boundary_condition_onemesh` 跳过多孔块
   （物理边界在多孔求解器内部施加），否则其 LS 类面码会被当作可压边界报错。
3. `src/sub_update_buffer_mpi.f90`：`Ts_send_mpi/Ts_recv_mpi` 仅在**两侧都是
   SOLID 或 POROUS**（携带骨架温度的块）时交换 Ts —— 修复流体块 Ts(恒 T_inf) 污染
   多孔骨架 ghost 的缺陷。

## 算例（基于 high_low_fluid_800K 网格）
- 可压区（Block2/3/4）：Ma=3、T_inf=800 K、Re/m=5.584e6（800 K 参考量：
  μ∞≈3.62e-5 Pa·s、ρ∞≈0.119 kg/m³、p∞≈2.73e4 Pa）。
- 原低速 Block1 → **BLOCK_POROUS**：ε=0.9、dp=0.1 m（软阻力，近似原 LS 通道流动以
  便于与码 12 对照）、hv=0（无相间换热）、骨架 k_s=13.4；Ts 初始 300 K。
- control：LS/porous 侧 ρ=0.317 kg/m³、p=2.7314e4 Pa、等效粘性 ν≈1e-2 m²/s（沿用
  high_low 网格稳定性处理）；`Iflag_vtk_onefile=0`。
- 运行：`mpirun -np 1 ./opencfd-ec1.16a.out`（t_end=150，单进程约 2 min）。

## 验收结果（t_end=150，exit 0，无 NaN）
| 块 | 码12（LS 替身）温度范围 K | 码19/PhaseA（多孔）温度范围 K |
|---|---|---|
| Block1（低速/多孔流体相） | 300.03–800.95 | 300.02–801.35 |
| Block2（可压） | 799.41–2078.98 | 799.41–2078.98 |
| Block3（可压，界面块） | 792.41–1285.90 | 792.91–1288.02 |
| Block4（可压） | 799.41–2079.09 | 799.41–2079.09 |
| 多孔骨架 Ts | —（无骨架） | **300.0（恒定）** ✔ |

- 可压三块温度范围与码 12 试跑几乎一致（Block3 差 <0.2%）；多孔流体相亦近似，
  说明耦合框架未破坏流场。
- 修复前 Ts 污染（1–800 K）已消除；`fluid_solid`（码 11）回归 exit 0、残差无回退。

## 局限 / 下一阶段
- 本阶段界面仅为"流体 ghost 交换"（剪切/接触型），未包含：
  阶段 B 冷却剂吹出（G=2 kg/(m²·s) 从多孔热端进入可压区）、
  阶段 C 骨架界面热平衡（T_w=Ts 显式交错供热）。
- 多孔块在界面的骨架 Ts 目前绝热（Ts≡300），符合阶段 A 设定。
