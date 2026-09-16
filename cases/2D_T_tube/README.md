# 2D_T_tube —— 2D T 形管（**当前不可正确运行，仅存档/**待修**）

目标：用 AC（`LS_Algorithm=3`）算 2D T 形管（一进两出），与文献对比分流比、两出口
Poiseuille 剖面与压降/损失系数 K。

## 现状（2026-09-15 核查）：**能加载、能跑，但当前设置不收敛**

### 文件约定（重要，2026-09-15 澄清）
- `bc3d.inp`：**物理边界** + **同类块（同一 BLOCK 类型）之间的对接**。内边界用负码，
  连接号与各维正负标记写在 `kb/ke` 列（`-1 -2` 形式）；转换器
  `sub_convert_inp.f90::Convert_bc` 自动求出对侧范围 `ib1..ke1`、`nb1`、`face1` 与
  连接次序 `L1..L3`，写入 `bc3d.inc` ⇒ **同类块算例只需要 `bc3d.inp`**。
- `bc3d_interface.inp`：只用于**跨类型 block**（低速↔多孔、流体↔固体等）的对接关系；
  单一块类型时该文件不参与——缺省时程序打印
  `bc3d_interface.inp not found, skip interface`，同类界面仍由 `bc3d.inc` 正常解析
  （实测：`channel_ac_2blk` 删掉该文件后仍打印 `shared same-class interfaces = 1` 并正常收敛）。
  本目录**没有**（也不需要）该文件。

### 可运行性（`run_probe.log`，`t_end=100`）
- 4 块网格读入正常，`Num_Cell=21600`（块节点 **61×61 / 61×101 / 61×101 / 101×61**）；
- 同类界面**全部解析成功**：启动打印 `AC: shared same-class interfaces = 3`；
- 无报错，能写出 `flow3d.dat` / `flow3d_block_*.vtk`。

### 但当前不收敛（同一日志）
100 外步后 `res_q≈3.0e2`、`res_m≈7.7e1`（对比 `channel_ac_2blk` 同阶段为 1e-2 量级），
场中出现 `u_min≈−3.8`（blk1）、`−7.5`（blk4），数百个 `u<0` 单元 ⇒ 属**算例设置/几何层面**
的问题，与接口文件无关。

**残差趋势（关键线索）**：第 20 外步 `res_q≈1.6e1` → 第 100 外步 `≈3.0e2`，即**增长**
（不是在消化一个大初值，而是**真不稳定**）⇒ 优先按下面的"参数"与"边界核对"两项定位。

> 该算例的 `bc3d_interface.inp/.inc` 已按文件约定删除（本算例四块同类型，对接关系全部
> 由 `bc3d.inp` 提供）。删除后实测：打印 `bc3d_interface.inp not found, skip interface`，
> 但**3 个同类界面仍全部解析**（`shared same-class interfaces = 3`），运行正常。

## 待办（按顺序）

- [ ] **边界与几何核对**：逐面核对 `bc3d.inp` 各块的物理边界码（入口 `1`/`5`、出口 `6`、
      壁 `2`、对称 `3`）与实际几何是否一致；T 形交汇处 3 个界面在拐角相交，需确认交汇点
      单元既没有重复计入也没有漏算（可参考 `channel_ac_4blk` 的角点自检做法）；
- [ ] **压力锚定**：确认至少一个出口带出口码（`BC_LS_Outlet`），使耦合压力有锚点
      （`LS_P_in=LS_P_out=0` 时压力水平由网格级锚定 `ac_have_pressure_bc` 决定）；
- [ ] **参数**：本算例 `AC_CFL=5`、`AC_CFLv=0.1`（远小于槽道算例的 20）、`AC_Flux=1`(Rusanov)、`AC_Recon=1`(MUSCL2)；
      （注：`AC_Flux=3` AUSM+ 已于 2026-09-17 从求解器中移除，详见 `docs/工作日志.md` 该日条目）
      建议先小步长短跑，确认 `res_q` 单调下降，否则先降 `AC_CFL` 定位不稳定源；
- [ ] 收敛后再验收：入口与两出口质量守恒、分流比、两出口 Poiseuille 剖面、压降/损失系数 K；
- [ ] 文献参考数据需用户指定（仓库内无 T 形管参考数据）。

## 运行时设置（多块 AC 的强制要求）

`control.ec` 已设为多块同步调度（**不要**改成 `AC_Max_Iter=60000, t_end=1`）：

```
Kstep_save=3000
t_end=3000.d0
AC_Max_Iter=20
```

**依据**：`AC_Max_Iter` 大 + `t_end=1` 会让每个块单独迭代到收敛、期间界面通量被冻结，
邻块首列被锁进假解（残差 ~1e-7 但场量非物理）。完整机理、判定链与验证数据见
`docs/工作日志.md` 2026-09-15 条目（含"补充（同日 2）"）与
`docs/AC-FV-后续工作.md` §17、§20；共形多块（2 块、2×2 四块含角点）的验收见
`cases/channel_ac_2blk/README.md`、`cases/channel_ac_4blk/README.md`。
