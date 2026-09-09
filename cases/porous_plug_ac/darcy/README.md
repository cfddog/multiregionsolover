# porous_plug_ac/darcy — AC 求解多孔介质流动(Darcy–Brinkman plug)

`BLOCK_POROUS` 用 **AC 方式**求解(与 `cases/porous_plug/darcy` 同网格/BC/物性,SIMPLE 对照)。

## 设置
- 同 SIMPLE darcy 子例:压差驱动 `LS_Inlet_Type=3`,`LS_P_in=0.25`、`LS_P_out=0`,
  ε=0.9、d_p=0.1 m、μ=0.02,plug 长 8 m × 高 1 m(40×20 单元);`If_viscous=1`
  (AC 的粘性项需显式开启)。
- AC:`LS_Algorithm=3`、`AC_beta=2000`、`LS_U_in=0.5`(仅作速度标度:
  Uref=0.5 → **β≈500 m²/s²**,与压差尺度 √(Δp/ρ)=0.5 一致)、
  `AC_CFL=20`、`AC_Tol=1e-7`。

> **修正记录(2026-09-11)**：早期版本 `LS_U_in=0.01` 使 β 实为
> `2000·(0.01)²=0.2`(运行日志 `beta=0.2`,非文档假设的 β≈500),小 β 下
> AC-Rusanov 连续方程数值耗散过大,稳态解受端部伪影偏移,速度偏高 +2.7%。
> 修正 `LS_U_in=0.5` 后(物理 BC 不变)β≈500,结果见下,详见
> `../hstudy/README.md`(网格收敛表,80×40 验收 AC vs SIMPLE <1%)。

## 结果(β=500 修正口径)
- 收敛:5 172 次伪时间迭代到 res_q≈2.9e-9 / res_m≈9.8e-8(β=500 后迭代数反而
  由 β=0.2 的 ~99k 显著下降;CFL=20);
- u 断面均值 ≈ **7.958e-3 m/s**;拟合 dp/dx ≈ **−3.126e-2 Pa/m**
  (理想 = −Δp/L = −3.125e-2,L_eff=7.998 m → 端部无伪影);
- SIMPLE 参考(40×20):u_bar=7.847e-3、dp/dx=−3.080e-2 → AC 偏高 ≈ **+1.4%**,
  该差主要为 SIMPLE 自身粗网格偏差(其 dp/dx 距理想 −1.44%,AC 距理想 −0.02%);
- **80×40 同网格验收(见 `../hstudy/h80`)**:AC u_bar=7.9089e-3 vs SIMPLE
  7.8927e-3(+0.20%),dp/dx −3.1252e-2 vs −3.1188e-2(−0.21%)→ **<1% 达成 ✅**。

## 复现
```bash
mpirun -np 1 ../../../../src/opencfd-ec1.16a.out   # 目录内
python3 ../../porous_plug/plot_validation/pv_plug_compare.py \
    --simple ../../porous_plug/darcy --ac . --pin 0.25   # 与 SIMPLE 对比度量
```
图:`validation_plug_compare_40x20.png`(本网格 SIMPLE vs AC 并列)。

## 状态/后续
- [x] AC-porous 流动(阻力入对角)已收敛、量级正确
- [x] plug/darcy 定量收口:β 配置修正后 AC≈理想解;80×40 网格 AC vs SIMPLE
      u_bar/dp_dx 均 <1%(2026-09-11,`../hstudy`)
- [x] LTNE 双温度标量对(`src/sub_porous_ac.f90::ac_porous_ltne`)与 ltne_1d 解析验证
      (见 `cases/porous_ltne_1d/run_ac`)
