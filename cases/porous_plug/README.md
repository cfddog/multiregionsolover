# Porous Plug / LTNE 验证算例（首版功能自检）

单块 `BLOCK_POROUS` 多孔介质，直接借用 low-speed 的进出口 BC 类型
（`BC_Inflow/BC_Outflow`），检验首版多孔求解器的三件事：

1. **多孔流动**：Darcy（+Forchheimer）阻力项已进入 SIMPLE 动量方程
   （与压力梯度一样按体积积分后加入 `a_p`）。实测入口速度加倍时，
   Δp 与域内速度近似线性加倍。
2. **LTNE 两温度模型**：`U(5)=Tf` 与 `B%Ts` 通过 `hv` 交换；Ts 方程带
   `(1-eps)*k_s` 导热，`solid_bc.inp` 给骨架壁面热边界。
3. **输出**：`flow3d.dat`（原始量 + 后处理 vtk）、`Ts_block_1.dat/vtk`。

## 文件说明

- `Mesh3d.x`：`../lidcavity/gen_mesh.py --nx 41 --ny 21 --nz 2 --Lx 8 --Ly 1` 生成（40×20 单元，准 2D）
- `bc3d.inp`：x- 入口、x+ 压力出口、y± 壁面、z± 对称
- `material.in`：块 1 类型 = 3（porous），帧材料 `rho_s/Cp_s/k_s`
- `porous.inp`：`块号 eps dp(m) hv(W/m³K)`；hv=0 退化为冻结骨架（纯流动/单温度）
- `solid_bc.inp`：`面号 Tw(K) Qw(W/m²)`，给 porous 块的骨架热边界（复用固体热边界文件）

## 控制文件要点

- 求解参数组：`Porous_Max_Iter`、`Porous_Tol`、`Porous_alpha_Ts`、`Porous_T_ref`
- 流体物性沿用 `LS_*`（ρ、μ、k、Cp、T_ref、进出口/壁面参数、松弛）
- 对流格式 `LS_Scheme`（首版建议 1）；压力速度耦合 `LS_Algorithm`（1=SIMPLE，2=SIMPLEC）

## 运行

```bash
mpirun -np 1 ./opencfd-ec1.16a.out
```

## 验证结果（Python 绘图检查）

绘图脚本位于 `plot_validation/`，在装有 numpy/matplotlib 的环境中运行，例如：

```bash
cd cases/porous_plug
conda run -n num_python python3 plot_validation/pv_main.py darcy    validation_darcy.png          --mode darcy --eps 0.9 --dp 0.10 --mu 0.02
conda run -n num_python python3 plot_validation/pv_main.py velinlet validation_darcy_velinlet.png --mode darcy --eps 0.9 --dp 0.10 --mu 0.02 --u_in 0.01
conda run -n num_python python3 plot_validation/pv_main.py .        validation_ltne.png          --mode ltne
```

三种情形（本目录 `Mesh3d.x` 为 41×21×2 均匀网格、40×20 单元、ε=0.9、dp=0.10 m → Ergun K=4.86×10⁻³ m²，μ=0.02）：

| 配置 | 收敛 | 关键结果 |
|---|---|---|
| `darcy/`（压差驱动，LS_Inlet_Type=3, LS_P_in=0.25） | 230 步, res=1e-9 | p 沿程线性，拟合 dp/dx=−3.08e-2 Pa/m；u_max=9.17e-3、柱均 u=7.85e-3，与按拟合 dp/dx 的 Darcy–Brinkman 解误差 **+0.48%** |
| `velinlet/`（速度入口 u_in=0.01） | 643 步, res=1e-8 | 局部 Brinkman 一致（柱均 u vs 解析 +0.51%）；但 **柱均 u/u_in=0.722**（约 28% 通量缺口） |
| LTNE（本目录, hv=50, Tw=400 K, LS_U_in=0.05） | 782 步 | Tf∈[311.6,400.0]、Ts∈[380.8,400.0] K；中心线在下游迅速趋近 400 K，近壁 Ts≈400 |

结论：

- 压差驱动的多孔塞（经典 Darcy–Brinkman 定解）定量吻合（<0.5%），说明多孔阻力源项、粘度项、压力梯度在 SIMPLE 离散中实现正确；
- **速度入口 + 强多孔阻力**情形下，SIMPLE 得到的域内平均流量低于指定入口速度（u_bar/u_in≈0.72），而局部压力-速度关系仍与 Brinkman 一致 —— 说明速度入口在阻力主导下并未像预期那样“锁定”总流量（入口面质量通量与压力修正/动量耦合需要进一步排查）。因此 plug 类验证建议使用**压差驱动**入口（LS_Inlet_Type=3），或多块结构的“低速-多孔”界面。
- LTNE 温度场有界、无 NaN，骨架/流体两温度趋势合理；进一步定量验证需要与 Wakao/二维 LTNE 解析解对比（下一步）。

图：`validation_darcy.png`、`validation_darcy_velinlet.png`、`validation_ltne.png`。

## 待办（下一阶段）

- 排查“速度入口 + 强阻力”下总通量与入口速度不一致的问题（入口面通量/压力修正耦合）
- LowPorous / SolidPorous / PorousPorous 界面耦合（接口码 16/17/18 已定义）
- Ts 缓冲交换按块类型过滤（只交换 solid/porous 之间），并核对耦合与缓冲先后顺序
- 进出口处可采用“纯流体(lowspeed)-多孔-纯流体”三块结构替代直接把 BC 加在多孔块上
- 后处理输出 p/Ts 进重启文件
