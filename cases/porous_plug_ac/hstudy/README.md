# porous_plug_ac/hstudy — plug/darcy 网格收敛与 SIMPLE 定量对照（P1 收口，2026-09-11）

均匀网格 20×10 / 40×20 / 80×40（Lx=8 m、Ly=1 m、z 1 层），物理参数与 BC 同
`cases/porous_plug/darcy`（ε=0.9、d_p=0.1 m → Ergun K=4.86e-3 m²、μ=0.02，
压差驱动 `LS_Inlet_Type=3`、`LS_P_in=0.25`→`LS_P_out=0`，上下无滑移壁，z 对称）。

| 子目录 | 求解器 | 说明 |
|---|---|---|
| `<h>/simple`   | SIMPLE（`LS_Algorithm=1`） | 即 `porous_plug/darcy` 参考口径（存档值一致） |
| `<h>/ac_b500`  | AC（`LS_Algorithm=3`） | **修正口径 β=500**：`LS_U_in=0.5` 作 Uref 标度（原存档误留 0.01 → β=0.2） |

## 结果（压差拟合窗口与 `pv_darcy.py` 相同：0.25–0.85 L；ū=列均值窗口平均）

| 网格 | SIMPLE ū | SIMPLE dp/dx | SIMPLE L_eff | AC(β=500) ū | AC dp/dx | AC L_eff | AC 相对 SIMPLE (ū, dpdx) |
|---|---|---|---|---|---|---|---|
| 20×10 | 7.5173e-3 | −2.89999e-2 | 8.621 m | 8.0986e-3 | −3.12856e-2 | 7.991 m | +7.73% / −7.88% |
| 40×20 | 7.8471e-3 | −3.08037e-2 | 8.116 m | 7.9583e-3 | −3.12592e-2 | 7.998 m | +1.42% / −1.48% |
| 80×40 | 7.8927e-3 | −3.11877e-2 | 8.016 m | 7.9089e-3 | −3.12523e-2 | 7.999 m | **+0.20% / −0.21% ✅** |

理想（压力钉扎两端面 x=0/8）dp/dx = −Δp/L = −3.125e-2 Pa/m。
SIMPLE 的 dp/dx 由下方以 ~二阶收敛至理想（L_eff 8.62→8.12→8.02 m）；
AC(β=500) 的 L_eff≡8.0 m、dp/dx 距理想 ≤0.06%，几乎与网格无关。

## 结论 / 根因

- 原存档 AC plug 结果（`docs/工作日志.md` 2026-09-09，ū=8.06e-3、dp/dx=−3.17e-2，
  相对 SIMPLE **+2.7%**）由 **β 配置错误**造成：`β = AC_beta·Uref²`，而 Uref 取
  `LS_U_in=0.01`（`LS_Inlet_Type=3` 下该量仅作速度标度占位）→ **β 实为 0.2**
  （运行日志 `beta=0.2`），非工作日志假设的 β≈500。
- 小 β 时 AC-Rusanov 连续方程数值耗散（特征速度 λ≈√β）相对物理项 β·u 过大，
  稳态解受其端部伪影偏移 → 速度偏高。β=500 后 AC 解复现理想压差解，收敛迭代数
  反而由 ~99 496（40×20，β=0.2）降至 5 172（β=500；80×40 为 10 148）。
- 40×20 上 AC(β=500) 与 SIMPLE 仍差 ~1.4%：该差为 **SIMPLE 自身粗网格偏差**
  （SIMPLE@40×20 距理想 −1.44%，AC 距理想 −0.02%），非 AC 误差。
- **80×40 验收（本目录 h80，两求解器同网格重跑）：ū +0.20%、dp/dx −0.21%，均 <1% ✅**，
  且两者同趋理想解 → P1 plug/darcy 定量收口达成。

图：`h40/validation_plug_compare_40x20.png`（存档网格对比）、
`h80/validation_plug_compare_80x40.png`（80×40 验收）。

## 复现

```bash
# 每格跑 SIMPLE / AC（control.ec 已含对应设置，二进制为 src 编译所得）
cd <h>/simple   && mpirun -np 1 ../../../../src/opencfd-ec1.16a.out
cd <h>/ac_b500  && mpirun -np 1 ../../../../src/opencfd-ec1.16a.out
# 对比度量（统一窗口/口径，含端部诊断）
python3 ../../porous_plug/plot_validation/pv_plug_compare.py \
    --simple h80/simple --ac h80/ac_b500 --pin 0.25
```
