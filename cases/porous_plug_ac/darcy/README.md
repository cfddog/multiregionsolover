# porous_plug_ac/darcy — AC 求解多孔介质流动(Darcy–Brinkman plug)

`BLOCK_POROUS` 用 **AC 方式**求解(与 `cases/porous_plug/darcy` 同网格/BC/物性,SIMPLE 对照)。

## 设置
- 同 SIMPLE darcy 子例:压差驱动 `LS_Inlet_Type=3`,`LS_P_in=0.25`、`LS_P_out=0`,
  ε=0.9、d_p=0.1 m、μ=0.02,plug 长 8 m × 高 1 m(40×20 单元);`If_viscous=1`
  (AC 的粘性项需显式开启)。
- AC:`LS_Algorithm=3`、`AC_beta=2000`(Uref=√(Δp/ρ)≈0.5 → β≈500 m²/s²)、
  `AC_CFL=20`、`AC_Tol=1e-7`。

## 结果
- 收敛:99 496 次伪时间迭代到 res_q≈2.6e-8 / res_m≈1.0e-7(β 太小时收敛极慢,经验见
  `docs/工作日志.md`);
- u 断面均值 ≈ **8.06e-3 m/s**(横向变化 ~1%);拟合 dp/dx ≈ **−3.17e-2 Pa/m**;
- SIMPLE 参考:u_bar=7.85e-3、dp/dx=−3.08e-2 → AC 偏高 ≈ **+2.7%**(尚未达到 <1%,
  待出口压力 BC 口径/压差拟合一致化后复核)。

## 复现
```bash
mpirun -np 1 ../../../../src/opencfd-ec1.16a.out   # 目录内
```

## 状态/后续
- [x] AC-porous 流动(阻力入对角)已收敛、量级正确
- [ ] LTNE 双温度标量对(`src/sub_porous_ac.f90::ac_porous_ltne` 占位)与 ltne_1d 解析验证
- [ ] plug/darcy <1% 定量收口
