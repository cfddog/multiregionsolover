# square_duct_1977 —— 90° 弯曲方管（Taylor–Whitelaw–Yianneskis 1982 / Humphrey 1977）

参考：**Taylor, A. M. K., Whitelaw, J. H., Yianneskis, M. (1982)**, *Curved ducts with strong
secondary motion: velocity measurements of developing laminar and turbulent flow*,
J. Fluids Eng. **104**, 350–359.（层流 Rc/D=2.3 基准源于 Humphrey, Taylor & Whitelaw 1977。）

## 几何（由 `Mesh3d.dat` 坐标反推，已核验）

| 块 | 节点 | 坐标范围 | 说明 |
|---|---|---|---|
| blk1 | 41×41×41 | x∈[0, 1.8]（轴向）、y,z∈[0, 0.04] | 上游直段 45W |
| blk2 | 41×41×41 | y∈[0.112, 1.312]（轴向 30W）、x,z∈[0, 0.04] | 下游直段（转过 90°） |
| blk3 | 11×41×41 | x∈[1.8, 1.912]、y∈[0, 0.112]、z∈[0, 0.04] | **弯段块**，跨度 R_c+W/2=0.112 |

- **管宽 W = 0.04**（无量纲），断面 **40×40 单元**（Δ=W/40=1e-3）；上游段轴向 45W、下游 30W；
- 弯段跨度 0.112 = R_c + W/2 ⇒ **R_c/W = 2.3**（经典弯曲方管基准）；
- 2 个同类（`bc<0`）内连接：blk1(i+)↔blk3(i−)、blk3(i+)↔blk2(i−) ⇒ 启动打印
  `AC: shared same-class interfaces = 2` ✓；`bc3d.inp` 中无类型 0 的条件 ✓。

## 网格文件（重要）

原始 `Mesh3d.x` 实际含 **15 块**（前 4 组是其它 2D 算例的块，**最后 3 块**才是本算例：
41³、41³、11×41×41，与 `bc3d.inp` 的三块逐一对应），且记录布局是**合并式**（每块 x,y,z 在
同一条记录里）。本算例**没有** `partation.dat`，代码会走自动分块路径，而该路径按
**格式化(ASCII)** 读取 `Mesh3d.dat`，因此：

```bash
# 用 util 脚本把 15 块中的最后 3 块抽成 ASCII 版 Mesh3d.dat
python3 extract3.py          # 本目录内的脚本
# control.ec 中需设：Mesh_File_Format=1（ASCII Mesh3d.dat）
```

## 运行（多块 AC 必须"小内迭代 + 多外步"）

```
Mesh_File_Format=1
LS_rho=1.0, LS_U_in=1.0, LS_P_in=LS_P_out=0     ! 无量纲
LS_mu = W/Re = 0.04/Re                          ! Re=790 ⇒ 5.0633e-5
LS_Algorithm=3, AC_beta=1, AC_CFL=20, AC_CFLv=0.1
AC_Max_Iter=20, t_end=3000~20000, Kstep_save=1000
../../src/opencfd-ec1.16a.out
```
（`AC_Max_Iter` 大 + `t_end=1` 会冻结界面通量而产生假解，见 `docs/工作日志.md` 2026-09-15。）

## 验收判据（对应文献）

- 弯段后二次流（Dean 涡）**横向速度剖面**与测量站（文献给的若干 θ 站位）对比；
- 弯段**内外壁压力分布**、弯段前后**摩擦系数/压降**；
- 直段充分发展 Poiseuille 剖面（上游 45W 处应已充分发展，Re=790 的入口长度约 0.06Re·W≈1.9）。

首轮建议先用 **Re=790（层流）** 跑通并检查质量守恒，再按文献做定量对比。
