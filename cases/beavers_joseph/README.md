# Beavers–Joseph（BJ）两区块验证算例

标准 **Beavers–Joseph (1967)** 界面滑移问题的两块实现：
上方自由流体通道（`BLOCK_LOWSPEED`）+ 下方 Darcy–Brinkman–Forchheimer 多孔床
（`BLOCK_POROUS`），中间以显式界面（接口码 **16**，`BC_Interface_LowPorous`）
连接。压差驱动、层流（creeping），与解析解对比验证程序。

## 物理模型（与代码一致）

求解器（不可压、SI 单位）充分发展时退化为 1D：

| 区域 | 方程 |
|---|---|
| 自由流体 `0<y<h` | `μ u″ = dp/dx` |
| 多孔床 `−Hp<y<0` | `μ u″ − (ε² μ/K) u = dp/dx` |

其中 Ergun 渗透率 `K = dp_p² ε³/[150(1−ε)²]`；把床内方程写成 Brinkman 形式时
`K_eff = K/ε²`、边界层厚 `λ = √K_eff`、深床 Darcy 速度
`u_D = −(K_eff/μ) dp/dx`。

界面 y=0 处速度与剪应力连续（通过接口 ghost 交换实现），床底 y=−Hp 与
通道顶 y=h 为 no-slip 壁。对深度满足 `κHp ≫ 1`（κ=1/λ）的床，其流体区剖面
即经典 BJ 滑移解

```
du/dy|₀₊ = (u_int − u_D)/λ      (等效 BJ 滑移系数 α=1，基于 K_eff)
```

## 几何与参数（medium 网格）

| 量 | 值 |
|---|---|
| 通道高 h / 床厚 Hp / 长度 Lx / z | 1.0 / 1.2 / 8.0 / 1.0 m（准 2D） |
| ε, dp_p | 0.8, 0.5 m → K_eff=3.33e-2 m², λ=0.1826 m |
| κh, κHp | 5.48, 6.57（深床 BJ 极限成立） |
| ρ, μ | 1 kg/m³, 0.02 Pa·s |
| 驱动 | 压力入口 `LS_P_in=5e-3 Pa` → 出口 0（`LS_Inlet_Type=3`） |
| 网格 | nx=81, 流体 ny=41（dy=0.025）, 床 ny=49（dy=0.025） |
| SIMPLE | `LS_Max_Iter/Porous_Max_Iter=80`, 外层 250 步（块迭代） |

## 文件

- `gen_case.py`：生成 `Mesh3d.x`（两块）、`bc3d.inp`、`bc3d_interface.inp`
  （接口码 16）、`material.in`、`porous.inp`、`control.ec`
- `run_med/ run_coarse/ run_fine/`：三套网格的结果（flow3d_block_*.vtk）
- `plot_validation/`：
  - `bj_analytical.py` —— 有限床厚 Darcy–Brinkman 复合解析解（闭式，2×2 线性）
  - `pv_bj.py` —— 读 VTK、拟合 dp/dx、逐点对比、出图
  - `pv_vtk.py / snap_check.py / iface_check.py` —— VTK 读取与外层收敛/界面检查
- `validation_bj.png`、`validation_bj_<mesh>.png`

## 运行与验证

```bash
# 1) 生成中网格算例
python3 gen_case.py --nx 81 --nyf 41 --nyp 49 --lx 8 --tend 250 \
                    --maxiter 80 --tol 1e-6 --ksave 250
# 2) 运行
mpirun -np 1 ./opencfd-ec1.16a.out
# 3) 对比解析解（num_python 环境含 numpy/matplotlib）
conda run -n num_python python3 plot_validation/pv_bj.py
```

两求解器以“外层块迭代”耦合：每外层步各自 SIMPLE 内迭代后用
`couple_lowspeed_porous_interfaces` 交换界面 ghost，外层逐步收敛
（监测 `snap_check.py`，最后 50 步相对变化 < 4e-5）。

## 结果（x=Lx/2，中网格 dy=0.025）

| 量 | CFD | 解析 | 相对误差 |
|---|---|---|---|
| u(y) 剖面 L2 / L∞（全流域） | — | — | 4.6e-4 / 9.7e-4 |
| 界面滑移速度 u_int (m/s) | 3.2998e-3 | 3.2900e-3 | +0.30% |
| Darcy 平台 uD (m/s) | 1.019e-3 | 1.041e-3 | ~2%（深部拟合窗口依赖） |
| BJ α_eff（界面斜率反推） | 0.925 | 1.0 | 有限差分/床深效应 |

网格收敛（dy = 0.05/0.025/0.0125，见 `pv_bj.py --case run_*`）：误差随加密近二阶下降。

## 代码改动摘要

- `src/sub_multi_region.f90`：新增 `couple_lowspeed_porous_interfaces(nMesh)`
  处理 LOWSPEED–POROUS（接口码 16）共形界面，双向 ghost 填充
  （U(1..5)、B%p；不交换多孔 Ts）。
- `src/opencfd_ec3d_v1.16a.f90`：三个时间分支中在
  `couple_solid_solid_interfaces` 后调用新例程。

## 说明与局限

- 界面为平面、法向质量通量≈0（fully developed），耦合体现在切向速度
  （剪应力）与压力；入口/出口采用压力型 BC。
- ε/dp/hv 目前为整块常数；如需更一般的 3 区块（流体-多孔-流体）或非共形
  界面需扩展耦合（L1/L2/L3 映射与通量对偶）。


## run_ac —— LS_Algorithm=3(AC)双块复测(2026-09-10,粗网格)

`run_ac/`:流体块 81×21×2 + 多孔块 81×25×2(与 `run_coarse` 同网格),
`control.ec` 改 `LS_Algorithm=3` + `AC_CFL=20/AC_beta=1000(β≈1.1e-2)/
AC_Max_Iter=1`,默认时间分支每块求解后即时交换接口 16 ghost
(`opencfd_ec3d_v1.16a.f90`,`LS_Algorithm=3` 时块 Gauss–Seidel 排序),t_end=3e5 外步
(=每块 3e5 次 LU-SGS 扫掠的整耦合系统伪时间推进)。残差收敛
(res_q≈3.1e-5 / res_m≈2.3e-4,曲线见 `residual_history.txt`)。

与解析对比(`pv_bj.py --case run_ac`,x=Lx/2):

| 量 | AC | 解析(AC 拟合 dp/dx) | 相对 |
|---|---|---|---|
| dp/dx 拟合 (fluid/porous) | −5.92e-4 / −5.99e-4 Pa/m | — | 双块一致 |
| Darcy 平台 uD | 9.78e-4 m/s | 9.92e-4 | −1.5 % |
| u(y) 剖面 L2 / L∞ | 3.3e-2 / 7.8e-2 | — | (粗网格) |
| 界面滑移 u_i(线性外推) | 3.60e-3 m/s | 3.13e-3 | +15 % |
| BJ α_eff | 0.68 | 1.0(有限床/粗网格亦 <1) | — |

图 `validation_bj_ac.png`。**结论**:AC 接口16 双块可稳定收敛并复现 BJ 的
dp/dx 一致、Darcy 平台与剖面趋势;界面滑移存在 AC 系统性偏高(单元层面 ~6%,
外推口径 ~15%)——与 `porous_plug_ac/darcy` 的 AC 速度偏高 +2.7% 属同一已知
问题(见 `docs/工作日志.md` 2026-09-09/10 条目与 `docs/AC-FV-后续工作.md` §9),
需待 plug/darcy <1% 收口后复核细网格;对比用 SIMPLE 粗网格参考
`run_coarse/validation_bj_coarse.png`(u_i 误差仅 1.6%)。

复现:`cd run_ac && python3 ../gen_case.py --nx 81 --nyf 21 --nyp 25 --lx 8 &&
(按上方 control.ec 改 LS_Algorithm/AC_*/t_end)&& mpirun -np 1 ../../../src/opencfd-ec1.16a.out
&& python3 ../plot_validation/pv_bj.py --case .`。
