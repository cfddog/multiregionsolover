# Porous LTNE 1D — 发汗冷却局部热非平衡（LTNE）定量验证算例

以 `cases/solid_1d` 的一维结构网格为模板，把"固体导热杆"替换成 **单块
`BLOCK_POROUS` 多孔介质 + 空气渗流 + 两温度（LTNE）** 的经典**一维发汗冷却**
问题，并将数值解与两温度模型的**闭式解析解**定量对比。
参考口径：中国科学技术大学董文杰博士论文中的发汗冷却一维热非平衡模型（气体冷却剂，
无相变）；论文全文无法在线获取，故按下列自洽的经典两温度定解执行（若需逐字对齐论文，
可据该文档快速修改边界参数复算）。

## 1. 物理问题与定解条件

多孔平板厚度 `L=10 mm`，坐标 `y∈[0,L]`：

- `y=0`（冷端）：冷却剂（空气）以质量流率 `G=2 kg/(m²·s)`、温度 `Tin=300 K`
  注入（速度入口）；骨架（不锈钢框架）对 300 K 环境对流换热，`h_c=31.4 W/(m²·K)`。
- `y=L`（热端）：骨架受**单位总面积**热流 `q=2e6 W/m²`；流体出口（零梯度）。

介质：`ε=0.3`，骨架 `k_s=13.4 W/(m·K)`（不锈钢），`d_p=1.0 mm`。
流体为空气（不可压常物性）：`ρ·u=G=2 kg/(m²·s)`、`cp=1007 J/(kg·K)`；
流体轴向导热相对对流可忽略（算例取 `LS_k=0`，Schumann 型模型）。
体积换热系数 `hv=2.0e6 W/(m³·K)`（Wakao–Kaguei 标定：
`Re=G·d_p/μ≈108`、`Nu=2+1.1Re^0.6 Pr^(1/3)≈18.3`、
`hv=6(1-ε)·Nu·k_f/d_p²`，见 §6 与 `porous.inp`）。

稳态两温度方程（每体积）：
```
流体:  G·cp·dTf/dy        = hv·(Ts−Tf),      Tf(0)=300
骨架:  kse·d²Ts/dy² + hv·(Tf−Ts)=0,          kse=(1−ε)·k_s
       冷端 y=0: kse·Ts'(0) = h_c·(Ts(0)−300)
       热端 y=L: kse·Ts'(L) = q        (进入骨架)
```

## 2. 解析解

消元得 `Ts''' + A·Ts'' − B·Ts' = 0`，`A=hv/(G·cp)`，`B=hv/kse`：
```
Ts = C0 + C1·exp(l1·y) + C2·exp(l2·y)
Tf = Ts − (kse/hv)·Ts'' = C0 + (kse/(G·cp))·(C1·l1·e^{l1 y} + C2·l2·e^{l2 y})
l1,2 = (−A ± √(A²+4B))/2
```
由三条定解条件（`Tf(0)=Tin`；冷端 Robin；热端定热流）解 `C0,C1,C2`
（`validate_ltne.py` 内用 3×3 线性代数闭式求解）。

## 3. 运行与验证

```bash
cd cases/porous_ltne_1d
python3 gen_case.py 401 2 2     # nx=401(Δy=25 µm), 准一维横向 2 节点；也可 ny=nz=9
mpirun -np 1 ./opencfd-ec1.16a.out      # ~1–1.5 min（Porous_Max_Iter=150000）
conda run -n num_python python3 validate_ltne.py
```

`validate_ltne.py` 读取 `Ts_block_1.dat`（Ts，质心）与 `flow3d_block_1.vtk`
（Tf `temperature` 标量，质心），按同一 `x` 列与解析解比较，输出各站位、RMS/最大
误差、能量平衡（`q ≟ G·cp·(Tf(L)−300)+h_c·(Ts(0)−300)`），并绘
`validation_ltne.png`。

## 4. 结果（nx=401 生产网格，Δy=25 µm）

| 站位 | Tf 数值 | Tf 解析 | Ts 数值 | Ts 解析 | (Ts−Tf) 数值 |
|---|---|---|---|---|---|
| y=0    | 302.6  | 302.7  | 515.5  | 516.8  | 212.8 |
| y=L/4  | 540.4  | 542.0  | 597.2  | 598.5  | 56.8  |
| y=L/2  | 696.6  | 697.8  | 770.6  | 771.7  | 73.9  |
| y=3L/4 | 925.2  | 925.9  | 1040.2 | 1041.0 | 115.0 |
| y=L    | 1287.0 | 1287.4 | 1467.9 | 1468.5 | 180.8 |

- Tf RMS 误差 **1.17 K**（最大 1.80 K），Ts RMS 误差 **1.11 K**（最大 1.40 K）；
  相对量程(≈1e3 K) 最大误差 ≈ **0.18 %**。
- 能量平衡闭合：入 `q=2e6 W/m²`，出 `G·cp·ΔTf + h_c·(Ts(0)−300)≈1.995e6 W/m²`，
  相对缺口 **0.27 %**（残余未完全收敛项）。
- 两温度显著分离（冷端 `Ts−Tf≈213 K`、热端 ≈181 K），LTNE 效应被定量复现。
- 图：`validation_ltne.png`。

## 5. 网格收敛（纵向加密，误差随 Δy ~ 一阶递减）

| nx | Δy (µm) | Tf RMS (K) | Tf max (K) | Ts max (K) | 相对误差 |
|---|---|---|---|---|---|
| 101 | 100 | 4.19 | 6.61 | 5.05 | 0.68 % |
| 201 | 50  | 2.25 | 3.49 | 2.70 | 0.36 % |
| 401 | 25  | 1.17 | 1.80 | 1.40 | 0.18 % |

数值离散中流体对流为一阶迎风（其数值扩散 O(Δy)）是误差主源，故呈一阶收敛；
解析模型忽略流体轴向导热（`LS_k=0`），与代码一致。

## 6. 物性 / hv 标定（Wakao–Kaguei）

`G=2`、`d_p=1e-3`、空气 300 K：`μ=1.846e-5 Pa·s`，`Pr≈0.707`，`k_f≈0.0262 W/mK`。
`Re_d=G·d_p/μ≈108`；`Nu=2+1.1·Re^0.6·Pr^(1/3)≈18.3`；
`h=Nu·k_f/d_p≈480 W/(m²·K)`，`a_v=6(1−ε)/d_p=4200 m²/m³`
⇒ `hv=a_v·h≈2.0e6 W/(m³·K)`（取整写入 `porous.inp`）。

## 7. 代码改动（本算例配套，均向后兼容）

- `solid_bc.inp` 支持每面骨架热边界施加于**任意物理面**（含流体入口/出口），并新增
  **对流(Robin)** 型（行格式：`face Tw Qw [h_c T_inf]`；`Tw=0,h_c>0` 触发 Robin）。
  —— `src/sub_read_parameter.f90`, `src/sub_modules.f90`, `src/sub_init.f90`,
  `src/sub_porous.f90::porous_ghost`
- **修复**：`porous_energy`（流体相）中 `h_v·Vol·(Ts−Tf)` 的 `Vol` 未赋值
  （未初始化），流-固换热项实际失效 —— 补 `vol = B%Vol(i,j,k)`。
- **修复**：骨架面热流/对流 ghost 外推改用有效导热 `kse=(1−ε)k_s`（与原方程系数一致），
  使 `Qw` 成为"单位**总**壁面热流"；原按 `k_s` 外推导致实际入骨架热流仅为 `(1−ε)Qw`。
- 影响：现有 `porous_plug` LTNE 数值更新（见其 README）；`solid_1d`/`solid_annulus`
  等固体算例回归通过。

## 8. 文件

## 9. AC 版验证(LS_Algorithm=3,2026-09-10)

`run_ac/`(nx=401)、`run_ac_g201/`、`run_ac_g101/`:`control.ec` 改 `LS_Algorithm=3`
+ `AC_*`,AC 流场(速度入口 2 m/s)约 **4.5k 伪时间步收敛**,随后自动执行
`src/sub_porous_ac.f90::ac_porous_ltne`(AC 速度重建面质量流,复用 SIMPLE 温度扫掠)
收敛 Tf/Ts(100 步窗口 ΔT<1e-3 K 停)。各网格闭式解对比:

| nx | Tf RMS (K) | rel err | Ts RMS (K) | 能量缺口 | (SIMPLE Tf RMS) |
|---|---|---|---|---|---|
| 101 | 3.97 | 0.66 % | 3.69 | 0.98 % | 4.19 |
| 201 | 2.21 | 0.35 % | 2.08 | 0.52 % | 2.25 |
| 401 | **1.22** | **0.19 %** | 1.17 | 0.28 % | 1.17 |

与 SIMPLE 结果及一阶收敛趋势一致(流体一阶迎风的数值扩散为误差主源);
图 `validation_ltne_ac.png`(401)、`validation_ltne_ac_101/201.png`,
明细见各目录 `validation_ltne_ac.txt`。复现:`cd run_ac && cp ../material.in ../porous.inp ../solid_bc.inp .
&& ../gen_case.py 401 2 2 && mpirun -np 1 ../../../src/opencfd-ec1.16a.out`。


- `gen_case.py`：生成 `Mesh3d.dat/.x`、`bc3d.inp`（i- 入口码 5，i+ 出口码 6，
  j±/k± 对称码 3）
- `material.in`：块 1 类型 3，不锈钢骨架
- `porous.inp`：`1  0.30  1.0e-3  2.0e6`
- `solid_bc.inp`：面 1（冷端）Robin `h_c=31.4, T∞=300`；面 4（热端）`Qw=2e6`
- `control.ec`：速度入口 `LS_U_in=2.0`（`LS_rho=1.0` ⇒ 质量流率 2 kg/m²·s），
  `LS_k=0`（流体轴向导热忽略，Schumann 型）
- `validate_ltne.py`、`run_grid.sh`：解析对比与网格研究脚本
