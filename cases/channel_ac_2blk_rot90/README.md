# channel_ac_2blk_rot90 —— 索引旋转（$i\leftrightarrow j$）的两块 AC 通道

`cases/channel_ac_2blk` 的物理问题完全相同（$L_x=8$、$H=1$、$80+80\times40$ 单元、
入口 `u=1` 插塞流、出口给定压力、上下壁无滑移、`BLOCK_LOWSPEED` + `LS_Algorithm=3`），
但**第 2 块写成转置网格**：其 `i` 沿 $y$、`j` 沿 $x$，界面因此不是“i 面接 i 面”，
而是 **blk-1 的 $i+$ 面 $\leftrightarrow$ blk-2 的 $j-$ 面**。用来检验
`Convert_bc` 的轴交换描述子（`L1/L2/L3`）与界面通量在**旋转界面 + 跨进程**下是否正确。

## 三种映射（关键字：旋向）

| 目录 | 块 2 的索引走向 | `L` | 旋向 | 说明 |
|---|---|---|---|---|
| `ref_np1` | 不转置（左右两半） | `(1,2,3)` / `(1,2,3)` | 右手 | 参考解，与 `channel_ac_2blk` 同一离散 |
| `p_np1`、`p_np2` | $i\to+y$、$k\to-z$ | blk1:`(1,2,3)` blk2:`(2,+1,-3)` | 右手 | “正序”旋转 |
| `n_np1`、`n_np2` | $i\to-y$、$k\to+z$ | blk1:`(1,2,3)` blk2:`(2,-1,+3)` | 右手 | “反序”旋转（$i$ 沿 $-y$） |

**交换两个索引方向会翻转三元组 $(x_i\times x_j)\cdot x_k$ 的旋向**，所以必须再反转
第三个方向才能保持右手。`gen_case.py` 早期版本把 `p` 变体的 `k` 也写成 $+z$
（`L=(2,+1,+3)`），结果是一个**左手块**：两块覆盖范围正确、$|V|$ 正常，但体积/面积/`dtime_mesh`
全为负，第一步就发散（`u` 冲到 $10^{10}$）。现在：

* `gen_case.py::blk2_xyz` 已改为 `k -> -z`（`p`）/ `k -> +z` 且 `i -> -y`（`n`），两种变体
  都是右手；
* `src/sub_geometry.f90::check_mesh_quality_onemesh` 启动时逐块检查最小单元体积，
  $<0$ 就打印块号、最小体积与 $(i,j,k)$，然后 `MPI_Abort(3)`，不会再“跑出乱七八糟的数”。

## 文件

| 文件 | 作用 |
|---|---|
| `gen_case.py` | 生成 `Mesh3d_{p,n}.x`、`bc3d_{p,n}.inp`（含 `L1/L2/L3`）、`material.in` |
| `prep.sh` | 建 5 个运行目录（`ref_np1 p_np1 n_np1 p_np2 n_np2`，输入软链）；`./prep.sh [步数]` 写独立 `control.ec`；`FRESH=1` 清掉 `field_restart.dat` 与输出 |
| `analyze_rot90.py` | 按块内 `i/j` 走向自动判断轴，把转置块散回 160×40 物理网格；报覆盖率、`u` 值域/列和、界面跳变，以及两目录逐点 RMS |
| `run.slurm` | 集群作业脚本（`sbatch -n` = 该目录 `partation.dat` 的进程数，自动定位三级之外的 `bin/`） |

运行产物（`flow3d*`、`field_restart.dat`、`run*.log`、`analyze_*`）不入库。

## 运行

```bash
FRESH=1 ./prep.sh 40000
cd ref_np1 && ../../../bin/MultiPhysMixSolver > run.log 2>&1     # 或 mpirun -np 1
cd ../p_np2 && mpirun -np 2 ../../../bin/MultiPhysMixSolver > run.log 2>&1
python3 ../analyze_rot90.py p_np1                # 单目录体检
python3 ../analyze_rot90.py p_np2 ref_np1        # 与参考逐点比较
```

集群（`wuzh02`）：

```bash
bash util/sync_cluster.sh push && ssh wuzh02
cd ~/sundong/MultiPhysMixSolver && make -j4            # 不要在 PATH 里带 anaconda！
cd cases/channel_ac_2blk_rot90 && FRESH=1 ./prep.sh 40000
for d in ref_np1 p_np1 n_np1; do (cd $d && sbatch -n 1 ../run.slurm); done
for d in p_np2 n_np2;       do (cd $d && sbatch -n 2 ../run.slurm); done
```

**坑**：若提交前 `module load apps/anaconda3/2023.09`（为了 `numpy`），`sbatch` 会把
anaconda 的 `mpirun`（MPICH 内核）带进作业环境，与 hpcx 链接的可执行文件不匹配，
作业在 `MPI_INIT` 之前就报 `*** on a NULL communicator`。**提交作业的 shell 里不要有
anaconda**；分析时再单独 `module load`。

## 结果（2026-10-18，修复后；wuzh02，$4\times10^4$ 步）

全场指标（`analyze_rot90.py`，160×40 网格全部覆盖且每格只写一次）：

| 目录 | u_max | u_min | 出口列和 | 界面 max\|Δu\| | 界面 max\|Δp\| |
|---|---|---|---|---|---|
| ref_np1 | 1.495712 | 0.073187 | 1.000000 | 1.671e-3 | 4.95e-3 |
| p_np1 | 1.495712 | 0.073187 | 1.000000 | 1.671e-3 | 4.95e-3 |
| n_np1 | 1.495712 | 0.073158 | 1.000000 | 1.676e-3 | 4.98e-3 |
| p_np2 | 1.495712 | 0.073187 | 1.000000 | 1.671e-3 | 4.95e-3 |
| n_np2 | 1.495712 | 0.073158 | 1.000000 | 1.676e-3 | 4.98e-3 |

逐点比较：

| 对比 | RMS(u) | max\|Δu\| | max\|Δp\| | mean(p−p′) |
|---|---|---|---|---|
| ref_np1 vs p_np1 | **3.0e-16** | 1.1e-15 | 1.3e-15 | +4.8e-17 |
| ref_np1 vs p_np2 | 2.4e-10 | 4.2e-10 | 1.1e-9 | −3.9e-10 |
| ref_np1 vs n_np1 | 4.0e-4 | 9.9e-4 | 2.4e-4 | −2.9e-7 |
| ref_np1 vs n_np2 | 4.0e-4 | 9.9e-4 | 2.4e-4 | −2.9e-7 |
| p_np1 vs p_np2 | 2.4e-10 | 4.2e-10 | 1.1e-9 | −3.9e-10 |
| n_np1 vs n_np2 | 2.4e-10 | 4.2e-10 | 1.1e-9 | −3.9e-10 |

结论：

1. **`p` 变体与未旋转参考逐位相同**（3e-16）——轴交换描述子、界面通量（含 owner 单值
   路径的符号）、`p` 的 halo 交换在旋转界面下都对；
2. **`n` 变体与参考差 4e-4**：差值不是界面跳变，而是集中在被反转的块内部
   （$x>4$ 的 RMS $5.6\times10^{-4}$，$x<4$ 只有 $4.4\times10^{-5}$），属于“同一格式在
   不同索引走向下的舍入/残差地板差”（该算例 `res_q` 只到 $10^{-5}\sim10^{-4}$），
   不是映射错误；
3. **np=1 与 np=2 各互差 2.4e-10**（`AC_Max_Iter=1` + 每步刷新 ghost 的预期水平）；
4. 窗口内质量流：出口列和 $=1.000000$（与解析 $uH=1$ 一致），界面左右列和差 1.2e-3（0.12%，
   与参考同量级）。

墙钟（`-N 1`，含启动与 VTK/重启 IO）：ref_np1 4m19s、p_np1 4m22s、n_np1 6m01s、
p_np2 2m51s、n_np2 2m49s。
