# high_low_fluid — Ma=3 / 真实气体压力入口 / 单文件 VTK 试跑记录 (2026-09-07)

## 目标设置（沿用原网格/几何，仅改来流与低速入口）
- y>0 可压区（Block 2/3/4）：来流 Ma=3，Re/m=5.584e6（Lscale=1 m → Re=5.584d6），
  来流静温 T_inf=107 K，绝热壁 Twall=-1，层流（Iflag_turbulence_model=0）。
- y<0 低速区（Block 1，类型 2）：入口为压力入口（LS_Inlet_Type=3）。

## 真实气体换算（按 Re/m=5.584e6 + Sutherland 反推）
超声速流（107 K）：μ∞=1.716e-5·(107/273.15)^1.5·(273.15+110.4)/(107+110.4)
≈7.42e-6 Pa·s；a∞≈207.35 m/s；ρ∞=Re·μ∞/(Ma·a∞·Lscale)≈0.0666 kg/m³；
p∞=ρ∞·R·T_inf≈2046 Pa（作为低速区参考压力，LS_P_in=LS_P_out=2.046d3）。

低速区空气按 **T=300 K、p=p∞≈2046 Pa**：
- μ(300 K)=1.846e-5 Pa·s；ρ=p/(R·300)=0.0238 kg/m³；ν≈7.77e-4 m²/s；
- 雷诺数测算（L=0.1 m，块高）：U=1→129、U=3→386、U=5→643
  —— 按“入口几 m/s”确为低雷诺数层流，分子粘性量级即可；
  （若误取常压 101325 Pa 的 ρ=1.177 kg/m³，则同速下 Re≈0.6–3×10⁴，需另行说明。）

## 数值稳定性说明（为何 control 中不用分子 μ=1.846e-5）
低速块 Block1 网格质量差（mesh-quality：顶部最大折转角≈90°、网格因子≈7.7e8、
最小体积 6e-9 的退化单元）。分子 μ 下（Re 即便只有数百）SIMPLE 压力修正在该网格
上病态：cap=30 时外迭代~7 步发散；cap=5 时~29 步；cap=2、LS_Max_Iter=500 时内迭代
~200 步 p 仍发到 10¹⁸⁶。因此层流 LS 求解器在当前网格需用等效（涡）粘性保证数值收敛：
`LS_mu=2.38e-4 Pa·s`（ν=μ/ρ=1e-2 m²/s，与既有稳定多块算例一致，等效 Re_LS≈O(10²)）、
`LS_k=0.342`（Pr=μCp/k≈0.7）。分子 μ 的严格物理模拟需先加密/重建 Block1 网格并
重新标定 LS 松弛参数。

## 代码改动
1. `src/sub_multi_region.f90` — `couple_highlow_fluid_face`（可压-低速接口 码12）：
   - 密度换算基准由固定 1 kg/m³ 改为按 Re/m+Sutherland 实时反推 ρ∞
     （RHO_REF = Re·μ_SI(T_inf)/(Ma·a_ref·Lscale)）；
   - 低速块 U(2..4) 存“速度”而非动量：接口 ghost 赋值去掉 ρ 因子
     （原 rho·u 仅在 ρ=1 时恰好正确）；低速→可压 ghost 的换算同步修正；
   - 仅当低速块无压力 Dirichlet 物理面时才整场平移锚定 p；压力入口
     (LS_Inlet_Type=3) / 压力出口时跳过（避免入口 ghost 反射与清零锚定冲突）。
2. `src/sub_init.f90` — 低速块冷启动 U(2..4) 直接用速度（与“存速度”约定一致）。
3. `src/sub_modules.f90` + `src/sub_read_parameter.f90` — 新增控制量
   `Iflag_vtk_onefile`（默认 0）与 `Iflag_vtk_SI`（默认 0：可压块按无量纲输出；
   =1：可压块换算为 SI 物理量纲）。
4. `src/sub_IO.f90` — `Iflag_vtk_onefile=1` 时不写 `flow3d_block_*.vtk`，
   所有块合并写入单个 VTK 非结构网格 `flow3d.vtk`
   （UNSTRUCTURED_GRID 六面体；cell data: density/velocity/temperature/pressure）；
   `Iflag_vtk_SI=1` 时把可压块速度×U∞(=Ma·a∞)、密度×ρ∞、温度×T∞、
   压力×ρ∞U∞²（ρ∞ 按 Re/m+Sutherland 反推），使两区域同为 SI 量纲可比。

## 试跑结果（单进程，t_end=150，Iflag_vtk_onefile=1，Iflag_vtk_SI=1）
- 无 NaN；低速 SIMPLE 内迭代收敛（res_p ~1.3e-5），p(1,1,1)≈2046.00 Pa；
- 单文件 `flow3d.vtk`：18406 节点 / 8760 六面体（不再写 flow3d_block_*.vtk）；
  低速区(Block1)：ρ=0.0238 kg/m³（300 K 空气 @ p≈2046 Pa）、p∈[1926,2054] Pa、
  T∈[233,300] K、|v|≤30 m/s；
  可压区（SI）：|v|∈[13.7,654] m/s（自由流≈622 m/s）、ρ∈[0.033,0.142] kg/m³、
  T∈[107,284] K、p∈[2044,5671] Pa —— 与低速区同一物理量纲可直接比较；
- 个别格点 max 残差 O(1e-3) 来自仍在发展的超声速剪切/压缩结构（非定常）。

## 复现
```bash
cd cases/high_low_fluid
./opencfd-ec1.16a.out    # Iflag_init=0 冷启动（覆盖 flow3d.dat / flow3d.vtk）
```
旧 Ma=2 设置与历次日志存档在 `_backup_ma2/`。

