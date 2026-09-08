# channel_ac — 不可压槽道:虚拟压缩法(AC-FV + LU-SGS)考核

在 `feat/lowspeed-ac-fv` 分支上新增的 **人工可压缩(Artificial Compressibility) 求解器**
(`LS_Algorithm=3`,`src/sub_lowspeed_ac.f90`)首个考核算例。

## 算例配置(与 `cases/channel/type1_vel` 完全相同的网格/边界)

- 2D 槽道 x∈[0,8]、H=1,节点 161×41×2(准 2D,两侧对称);
- ρ=1 kg/m³、μ=0.02 Pa·s、入口 u=1 m/s(速度入口)→ Re_H=50;
- i- 速度入口、i+ 压力出口 p=0、上下无滑移等温/绝热壁;
- `LS_Algorithm=3`,`If_viscous=1`,`AC_beta=1`,`AC_CFL=20`,`AC_Max_Iter=60000`,
  `AC_Tol=1e-7`(无量纲残差)。

## AC-FV 方法(本目录对应的实现)

- 状态 `q=p/ρ,u,v,w`;伪时间方程:
  `∂q/∂τ + β∇·u = 0`、`∂u/∂τ + ∇·(uu+qI) = ν∇²u`,β=AC_beta·U_ref²;
- cell-center 有限体积,直接复用多块结构网格几何(`Vol/Si/Sj/Sk/ni…`)、ghost 与
  接口/输出框架;
- 内面无粘通量:**van-Leer 限幅的 MUSCL 重构 + AC-Rusanov 迎风**;
  边界面为显式 BC 通量(壁面仅压力、入口给定状态、出口给静压);
- 粘性:沿坐标面的正交近似 Laplacian(面梯度/面距);
- 时间:**伪时间隐式 LU-SGS 对称扫掠**(定常);`Time_Method=Time_Dual_LU_SGS` 时
  自动切双时间步(动量加 BDF2 物理源项)——DTS 已接线待专项算例考核;
- 壁面 ghost:无滑移壁作**全速度镜像**(含切向),对称面仅镜像法向分量。

## 结果(`analyze_ac.py` 读 `flow3d_block_1.vtk`)

| 测站 x | u_max | 断面平均 u | 与 Poiseuille `6y(1−y)` 的 RMS |
|---|---|---|---|
| 2.025(入口段) | 1.470 | 1.000 | 0.0152(发展段,非误差) |
| 4.025 | 1.494 | 1.000 | 0.0018 |
| 6.025(充分发展) | 1.495 | 1.000 | 0.0013 |
| 7.025(充分发展) | 1.495 | 1.000 | 0.0012 |

- 充分发展区:u_max=1.495(解析 1.5,偏差 ≈0.3%);断面平均质量流量 **=1.000**(守恒);
  剖面与 Poiseuille 解析解 RMS ≈ **0.1%**;
- 收敛:39 370 次伪时间迭代到 res_q/res_m < 1e-7(`AC_CFL=20`、β=U²)。

## 复现

```bash
mpirun -np 1 ../../src/opencfd-ec1.16a.out   # 目录内运行
python3 analyze_ac.py                         # 剖面 / 守恒 / Poiseuille 对比
```

## 说明与后续

- 一阶(Rusanov 无重构)版本数值耗散过大(等效 μ_eff≈2.6μ),换 MUSCL 后达上述精度;
- 已接线但**待专项考核**:双时间步(LU-SGS-DTS)、被动标量 T 求解、压力入口与多区域
  接口(12/13/16)复跑、方腔(Ghia)等(见分支开发计划)。
