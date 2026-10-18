# lidcavity_ac_46blk —— 4/6 块 AC 方腔（进程数与压力锚定敏感性回归算例）

同 `cases/lidcavity_ac_3blk` 的物理问题（Re=100、42×42 单元、`U_lid=1` 顶盖驱动方腔，
全封闭、四壁无滑移、三块同为 `BLOCK_LOWSPEED`），但换成两种二维拼块与**多种块编号**：

| 前缀 | 拼块 | 每块单元 | 块数 | 内部界面 |
|---|---|---|---|---|
| `m2` | 2×1（左右两半） | 21×42 | 2 | 1 条 i 面（x=0.5） |
| `m4` | 2×2 | 21×21 | 4 | 1 条 i 面 + 1 条 j 面 |
| `m6` | 3×2 | 14×21 | 6 | 2 条 i 面 + 1 条 j 面 |

编号标签（`gen_case.py::order_list`）决定**全局 1 号块落在哪个物理子域**，也就是压力基准锚点：

| 标签 | 块 1 位置 | 用途 |
|---|---|---|
| `a` | 左下角（行主序） | 基准编号，与 3 块算例一致 |
| `b` | 右上角（完全逆序） | 锚点换到对角 ⇒ 检验 `p` 基准是否只是规范（gauge） |
| `c` | 左下角，其余逆序 | 锚点不变、访问顺序打乱 ⇒ 检验顺序无关性 |
| `d` | 左下角（列主序） | 顺序改变 ⇒ 同上 |

**这个算例抓到的 bug**：`b` 编号让界面 owner（小块号）落在**高索引一侧**，暴露了
`ac_ifc_flux_own`/`ac_ifc_flux_nb`（owner 单值通量路径）把 owner 恒当作低侧、恒写
`RAC -= FLX*A` 的错误 —— 界面变成抽吸源（`u_max` 1.83 而非 0.93）。修复见
`src/sub_lowspeed_ac.f90` 与手册 §4.5「owner 单值通量路径的符号修正」。

## 文件

| 文件 | 作用 |
|---|---|
| `gen_case.py` | 生成 `Mesh3d_{m2,m4,m6}{a,b,c,d}.x`、`bc3d_*.inp`、`material_m2/m4/m6.in` |
| `prep.sh` | 建立 13 个运行目录（输入用软链）；`./prep.sh [步数]` 写独立 `control.ec`；`FRESH=1` 清掉上次的 `field_restart.dat` 与输出 |
| `analyze_mblk.py` | 把 `flow3d_block_*.vtk` 按全局网格拼回 42×42；输出覆盖率、Ghia 比较、界面跳变，以及两目录逐点 RMS |
| `run.slurm` | 集群作业脚本（`-n` = 该目录 `partation.dat` 的进程数，自动定位三级之外的 `bin/`） |

运行产物（`flow3d*`、`field_restart.dat`、`run*.log`、`analyze_*`）不入库。

## 运行

```bash
FRESH=1 ./prep.sh 10000                 # 生成网格 + 13 个运行目录（10000 步）
cd m6a_np1 && ../../../bin/MultiPhysMixSolver > run.log 2>&1
cd ../m6a_np6 && mpirun -np 6 ../../../bin/MultiPhysMixSolver > run.log 2>&1
python3 ../analyze_mblk.py . ../m6a_np1     # 与 np=1 逐点比较
```

**必用设置**（`control.ec` 已写好）：`AC_MPI_Sync=1`、`AC_Max_Iter=1`
（理由见 `cases/lidcavity_ac_3blk/README.md`）。`FRESH=1` 不可省：目录里若留着上次
的 `field_restart.dat`，程序会从该步读回并**立刻结束**，看到的是旧结果。

## 结果（2026-10-18，修复后；wuzh02，10⁴ 外步）

进程数无关性（强扩展，同一 42×42 网格；墙钟含启动/IO）：

| 布局 | np=1 | np=2 | np=3 | np=4 | np=6 |
|---|---|---|---|---|---|
| m4 墙钟 / 加速比 | 34.8 s / 1.00 | 12.7 s / 2.73 | —— | 9.85 s / 3.53 | —— |
| m6 墙钟 / 加速比 | 39.4 s / 1.00 | 14.7 s / 2.68 | 10.3 s / 3.83 | —— | 5.81 s / 6.78 |

| 对比 | RMS(u) | max\|Δu\| | max\|Δp\| |
|---|---|---|---|
| m6 np=1 vs np=2 | 3.60e-7 | 9.41e-7 | 3.1e-6 |
| m6 np=1 vs np=3 | 3.60e-7 | 9.41e-7 | 3.1e-6 |
| m6 np=1 vs np=6 | 3.60e-7 | 9.41e-7 | 3.1e-6 |
| m4 np=1 vs np=2 | 3.12e-7 | 1.06e-6 | 2.3e-6 |
| m4 np=1 vs np=4 | 3.12e-7 | 1.06e-6 | 2.3e-6 |

np=2/3/6（以及 np=2/4）互相**逐位相同**——正是 `AC_Max_Iter=1` + 每步刷新 ghost 的设计
目标。np=1 与 np≥2 的 1e-6 差别来自两条界面路径（np=1 走 owner 单值通量 `FSH`，
np≥2 闭腔走两侧各算 `ac_face_flux_rus1_own`），与 3 块算例同一量级。

压力锚定敏感性（**只换块编号**，物理问题完全相同）：

| 对比 | RMS(u) | max\|Δu\| | mean(p−p′) | max\|Δp\| |
|---|---|---|---|---|
| m2a vs m2b | 1.28e-6 | 2.39e-6 | +3.33e-2 | 3.33e-2 |
| m4a vs m4b | 1.92e-6 | 5.37e-6 | +8.41e-3 | 8.42e-3 |
| m6a vs m6b | 2.23e-6 | 5.78e-6 | +4.53e-2 | 4.53e-2 |
| m6a vs m6c | 4.70e-9 | 2.93e-8 | +2.16e-9 | 3.35e-8 |
| m6a vs m6d | 0（逐位相同） | 0 | 0 | 0 |

结论：**速度场与压力锚点无关**（差值 1e-6 量级，纯离散/迭代残差），压力只差一个
**常数**（`max|Δp| == |mean(p−p′)|`，即整场平移）；打乱访问顺序（`c`）连这个常数
差异都没有（1e-9，纯浮点舍入）。

与 Ghia(1982) 的对照（10⁴ 步；3 块参考算例 RMS u=0.02095 / v=0.04365）：

| 布局 | RMS u(x=0.5) | RMS v(y=0.5) | u 值域 |
|---|---|---|---|
| m2a | 0.02058 | 0.04410 | [−0.2368, +0.9260] |
| m4a | 0.02082 | 0.04363 | [−0.2304, +0.9258] |
| m6a | 0.02121 | 0.04311 | [−0.2268, +0.9254] |

界面跳变（界面两侧相邻列/行的 max|Δu|）：m2 1.03e-2、m4 1.03e-2 / 5.91e-3、
m6 1.87e-2 / 7.18e-3 / 5.88e-3，与 3 块算例的 9.3e-3 同一量级，是「界面刻意一阶」的
预期误差（`h=1/42`、`∂u/∂x~0.4` ⇒ `½h|∂u/∂x|≈5e-3`），不是缺陷。

## 集群运行（wuzh02，Slurm）

```bash
bash util/sync_cluster.sh push            # 本地 -> 集群（代码/输入/脚本）
ssh wuzh02
cd ~/sundong/MultiPhysMixSolver && make -j4
cd cases/lidcavity_ac_46blk && FRESH=1 ./prep.sh 10000
for d in m2a_np1 m2b_np1 m4a_np1 m4b_np1 m6a_np1 m6b_np1 m6c_np1 m6d_np1; do (cd $d && sbatch -n 1 ../run.slurm); done
for d in m4a_np2 m6a_np2; do (cd $d && sbatch -n 2 ../run.slurm); done
(cd m6a_np3 && sbatch -n 3 ../run.slurm); (cd m4a_np4 && sbatch -n 4 ../run.slurm)
(cd m6a_np6 && sbatch -n 6 ../run.slurm)
# 分析（集群上有 anaconda：module load apps/anaconda3/2023.09）
python3 analyze_mblk.py m6a_np1 m6b_np1
```
