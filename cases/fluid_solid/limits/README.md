# limits —— 多孔 ε→0 / ε→1 两个极限的短调度算例集

目的是给多孔块（`BLOCK_POROUS`，码 3）的两个解析端点做**可重复的回归**：

| 端点 | 期望 | 程序内的显式处理 |
|---|---|---|
| ε = 1 | 无骨架 ⇒ 阻力 → 0，块内动量/连续方程退化为纯 N–S；`Ts≡Tf` | Ergun `K=dp²ε³/[150(1−ε)²]→∞` ⇒ `c_drag→0`；界面骨架鬼点按 `(1−ε)k_s≤1e-6` 改用零梯度 |
| ε = 0 | 不可渗固体 ⇒ 速度**严格为 0**，只剩骨架导热 | `porous.inp` 接受闭区间 `[0,1]`；阻力按 `max(ε,1e-8)`；AC 对角 `1e30`；动量后**显式冻结** `u=v=w=0`（含鬼点层） |

物理背景与公式见 `docs/程序说明/chap05_porous.tex` §“ε→0 与 ε=1 两个极限”，
算例清单见 `chap09_case_fluid_solid.tex` §“极限算例”。

## 目录与对照关系

| 目录 | 块1 类型 | 用途 / 对照 |
|---|---|---|
| `run_ref_case1` | `BLOCK_SOLID`（码 1） | case1 的短调度参照（固体导热，`Ts` 即固体温度） |
| `run_ref_case2` | `BLOCK_LOWSPEED`（码 2） | case2 的短调度参照（冷却剂通道，纯 N–S） |
| `run_por_eps1` | 多孔 ε=1.0、`h_v=2e6` | 对 `run_ref_case2`：阻力应自动为 0 |
| `run_por_eps1m6` | 多孔 ε=1e-6 | 对 `run_ref_case1`：ε→0 的连续逼近 |
| `run_por_eps0` | 多孔 ε=0 | 对 `run_ref_case1`：精确的不可渗固体极限 |

五个目录的 `Mesh3d.x`（2 块 67×101×2 / 101×201×2）、`bc3d.inp`、`material.in`、
`porous.inp`、`control.ec` 均自带；唯一的差别是块类型/ε，以及 `run_por_eps0`
的块1 `j+` 面由 case3 的冷却剂入口（码 5）改为等温壁（码 2 + `solid_bc.inp`
面 5 = 300 K）——其 `bc3d.inp` 与 `run_ref_case1` 相比**只差接口码 19 vs 11**
（块1 `j=1` 面、块2 `i=101` 面）。

### 共同短调度（相对 `run_case*` 只压缩预算）

`Kstep_Couple_Comp=100`、`Kstep_Couple_Min=100`、`Niter_Couple_Warm=10`、
`Niter_Couple_Outer=2`、`Porous_Chunk_Iter=2`、`AC_Max_Iter=5000`、
`Porous_Max_Iter=3000`、`Solid_Max_Iter=3000 / Solid_Tol=1e-8`、
`Kstep_save=100`、`Iflag_restart=-1`（不读不写重启，保证各算例初场一致）。

## 运行

```bash
cd cases/fluid_solid/limits
bash run_limits.sh run_ref_case1 run_ref_case2 run_por_eps1 run_por_eps1m6 run_por_eps0
```

`run_limits.sh` 默认使用仓库根目录下的 `bin/MultiPhysMixSolver`（按脚本位置推导，
不依赖绝对路径），每个目录先清掉上一次的 `flow3d_node.dat / *.vtk / run.log` 等产物。
2026-09-26 一轮的墙钟时间（`status.txt`）：case1 25 s、case2 3.7 min、
ε=1 9.2 min、ε=1e-6 9.6 min、ε=0 1.0 min。

## 检查与对照

```bash
# 结构自检（块数/维度/点数/变量范围，需与 Mesh3d.x 一致）
python3 ../../../util/check_flow3d_node.py .

# 逐点对照（排除界面节点行后的指标更有意义）
python3 ../../../util/compare_flow3d_node.py ../run_ref_case2 . --blocks 1
python3 ../../../util/compare_flow3d_node.py ../run_ref_case1 . --blocks 1
```

`flow3d_node.dat` 为节点值、SI 单位（`d,u,v,w,T,Ts`）。2026-09-26 那一轮的文件由
**当时的工作树二进制**写出，头记录是过渡的 2 条写法（记录2 = 每块 4 个整数
`ni,nj,nk,nvar`）；2026-10-17 修复后新写出的文件是标准的 3 条头
（块数 / 各块维数 / 变量数=6）。两个脚本对两种布局都能读。

## 实测结果（2026-09-26 一轮，块1 节点值范围）

| 算例 | u [m/s] | v [m/s] | T [K] | Ts [K] |
|---|---|---|---|---|
| `run_ref_case1`（固体） | 0 | 0 | 300 → 800 | 300 → 800 |
| `run_ref_case2`（低速） | −0.134 → 0.687 | 0.273 → 17.71 | 203.7 → 299.8 | 同 T |
| `run_por_eps1` | **−16.91 → 17.64** | −10.76 → 6.77 | 300（未变） | 300（未变） |
| `run_por_eps1m6` | ±6e−10 | −4.5e−6 → 1.425 | 300（未变） | 301.4 → 302.9 |
| `run_por_eps0` | **0（精确）** | **0（精确）** | 300 → 303.1 | 300 → 303.1 |

全部 5 例 `EXIT=0`、`run.log` 无 `NaN`；`check_flow3d_node.py` 均 `STRUCTURE OK`。

### 结论

- ✅ **ε=0 的“速度严格为零”成立**：`run_por_eps0` 块1 的 `u,v,w` 在**全部节点**
  （含鬼点层写入的边界节点行）逐点等于 0，不是“很小”而是浮点 0。
- ✅ **ε=1 的阻力确实自零**：`run_por_eps1` 没有发散/NaN，块内出现 O(10) m/s 的流动，
  说明 Darcy/Forchheimer 源项没有把速度压死（对照 `run_por_eps1m6` 的 `u≈1e-9`）。
- 🟡 **ε=1 ↔ 低速的逐点定量对照尚未通过**：`run_por_eps1` 块1 的主速度分量在 **u**
  （±17.6 m/s），而 `run_ref_case2` 同块主分量在 **v**（17.7 m/s）；`--blocks 1`
  对照给出 `max|Δu|=17.6 m/s`、`max|Δv|=28.0 m/s`、`max|ΔT|=96.3 K`，
  差值**不限于界面一排**。
- 🟡 **ε=0 ↔ 固体的温度场也未重合**：`run_por_eps0` 的 `Ts` 只到 303 K，而
  `run_ref_case1` 到 800 K（`max|ΔT|=498.6 K`）。
- 🟡 **多孔块的能量方程在这套压缩调度下几乎不推进**：三个多孔算例的流体 `T`
  都停留在初值 300 K，`Ts` 最多只变化 3 K。

### 对上述负结果的分析（待排查）

1. **界面模型不同**：码 19 的**气体侧**是“无滑移等温壁 + 可选发汗”，壁温取多孔块
   的 `Ts`；码 11/12 的界面是“温度/热流连续”。因此 `run_por_eps0`（19）的冷壁
   （300 K）不会像 `run_ref_case1`（11）那样被 2000 K 的燃气加热——温度对照
   **本质上不可比**，只能验证“块内导热 + 速度冻结”。
2. **冷却剂可能无法穿过码 19 界面**：`run_por_eps1` 的流动方向与低速参照不一致
   （u 主导 vs v 主导），提示多孔侧在 j=1 界面处的处理接近“壁”（不能流出），
   流体只能在通道内回流；而码 12 允许流-流贯穿。于是“块内方程相同”成立，
   “同边界下的数值解相同”不成立。
3. **调度太短**：`Niter_Couple_Outer=2`、`Porous_Chunk_Iter=2` 时多孔段一共只被调用
   4 次，`Ts/Tf` 来不及响应（见上表），无法用稳态判据做定量比较。

**因此本算例集的定位是**：① 回归“ε 极限的显式分支不崩、不 NaN”；② 回归“ε=0 速度
逐点为 0”；③ 提供两个端点与参照算例的**现状快照**。它**不**作为“ε=1 复现低速、
ε=0 复现固体”的定量证据。

### 后续（若要拿到定量证据）

- 用**同一种界面模型**做对照：低速块 ↔ 多孔 ε=1 都接码 12；或两个多孔块
  （ε=1 ↔ ε=1e-6）互相接码 18/19，避免“壁 vs 连续”混入；
- 把 `Niter_Couple_Outer` 提到 ≥10、`Porous_Chunk_Iter` 提到 ≥10（或让 `t_end`
  跑到稳态），确认 `Ts/Tf` 收敛后再比较；
- 若确认码 19 的多孔侧确实不允许流出，应在 `docs/程序说明/chap05_porous.tex`
  与 `chap07_coupling.tex` 明确写成“壁 + 发汗”的适用边界，而不是“多孔壁通道”。
