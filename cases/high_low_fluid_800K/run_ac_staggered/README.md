# high_low_fluid_800K / run_ac_staggered — AC + 分段交错（interface 19）试跑

本目录是 `cases/high_low_fluid_800K`（Ma=3、T∞=800 K 可压高速区 Block2/3/4 +
下方低速 Block1）网格的 **code-19 多孔版试跑**：把 Block1 按
`cases/interface19_phaseA` 改成 `BLOCK_POROUS`（ε=0.9、dp=0.1 m、hv=0、骨架
ρs=8000/Cp=500/ks=13.4，流体相沿用原 LS 参数：ρ=0.3172 kg/m³、μ=3.17e-3 Pa·s
即等效 ν≈1e-2 m²/s、LS_Inlet_Type=3），低速区求解器用 **AC**（`LS_Algorithm=3`），
外层调度用**分段交错**（`Iflag_Couple_Scheme=1`）：先跑高速区 1000 步 → 跑低速
多孔区至收敛 → 再回高速区，如此迭代；前 3 轮每轮满 1000 步“暖机”，之后每轮
气体段步数减半精修（下限 `Kstep_Couple_Min=20`）。

## 代码改动（本试跑，向后兼容）
- `src/sub_multi_region.f90::run_staggered_fluid_porous`
  1. **code-19 界面搜索推广到全部 BLOCK_FLUID 块**：原实现只查“第一个”可压块，
     而 high_low_fluid_800K 是 4 块拓扑（Block1=POROUS，Block2/3/4=FLUID），
     code-19 面在 **Block3**（j- ↔ Block1 j+），非第一个 FLUID 块 → 会误报
     “no FLUID-POROUS interface found”。现改为扫描所有 FLUID 块找
     code-19→POROUS 共形配对（仍限定 gas j- / porous j+，face_s=2/face1=5）。
  2. **多孔段按 `LS_Algorithm` 分派**：多孔段原固定调 `porous_solver_one_block`
     （SIMPLE），改调 `solver_one_block` → `LS_Algorithm=3` 时走
     `porous_ac_solver_one_block`（AC 伪时间收敛 + LTNE），SIMPLE 路径不变。

## 关键运行设置（control.ec）
```
LS_Algorithm=3, AC_beta=10 (LS_U_in=10 作 Uref 标度 → β≈1000 m²/s²),
AC_CFL=20, AC_Max_Iter=20000, AC_Tol=1e-6, If_viscous=1
Iflag_Couple_Scheme=1, Kstep_Couple_Comp=1000, Niter_Couple_Warm=3,
Niter_Couple_Outer=12, Kstep_Couple_Min=20, Tol_Couple_Tw=1e-2,
Porous_Chunk_Iter=1, Porous_Max_Iter=30000, Twall_Couple_Init=300 K
```

## 试跑结果（2026-09-09，单进程 mpirun -np 1，~2 min，exit 0、无 NaN）
调度与 AC 路径均正常：每轮打印 `outer iter k : gas chunk steps = n`，
1000→500→250→125→62→31→20→20…；多孔段 AC-LU-SGS 每块跑到
AC_Max_Iter=20000（残差 res_q≈1.2e6、res_m≈2e4，因本网格 Block1 顶部退化
单元 + 无压差驱动而收敛慢，未达 AC_Tol）；随后 AC-porous LTNE 骨架 Ts 300 余
次扫掠收敛。界面量写入 `iface_couple.dat`（x_w、T_w、q_w、p_w，19 个面元）。

| outer | 气体段步数 | max|q_w| W/m² | T_w 范围 K | max|ΔT_w| K |
|---|---|---|---|---|---|
| 1 | 1000 | 2.09e6 | 301.5–301.6 | 1.56 |
| 3 | 1000 | 9.31e5 | 302.8–303.2 | 0.69 |
| 6 | 125 | 7.98e5 | 304.0–305.0 | 0.60 |
| 9 | 20 | 7.75e5 | 305.0–306.8 | 0.58 |
| 12 | 20 | 7.61e5 | 306.1–308.5 | 0.57 |

12 轮外层结束后**未收敛**（`max|ΔT_w|=0.57 K > Tol`）。T_w 单调爬升、q_w 单调
下降（2.09e6→7.6e5 W/m²），趋势正确但放松很慢 —— 物理原因：本构型沿用
interface19_phaseA 的 **hv=0 且骨架无底端散热 BC（solid_bc.inp 缺省全绝热）**，
多孔壁只“储热”无热沉，T_w 需缓慢爬向气体近壁温才有 q_w→0。这与
`porous_fluid_phaseB/staggered_tw/README.md` 无粘试跑“T_w 单调上升不可收敛”
属同类限制，只是本运行含粘性 + 界面温度反馈使 ΔT_w 逐轮递减。

## 结论 / 后续
- ✅ AC + interface-19 分段交错在 high_low_fluid_800K 4 块网格上调度与数据流
  正常（gas 等温壁 T_w → q_w/p_w 提取 → 多孔段 AC 求解 → Ts 返回 T_w）。
- ⚠️ AC-porous 段在 20 000 伪步内未达 AC_Tol（退化网格 + 无压差），每段跑满
  AC_Max_Iter 预算；如需“真收敛”建议：(a) Block1 重建/加密后处理退化单元；
  (b) 给多孔骨架加底端热沉（solid_bc.inp 等温或 Robin 冷端）或 hv>0 冷却剂
  通流（transpiration），否则 T_w 只能以储热方式逼近气体近壁温。
- 复现：`cd run_ac_staggered && mpirun -np 1 ./opencfd-ec1.16a.out`。
