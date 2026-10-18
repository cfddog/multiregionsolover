# solid_4blk —— 4 个**互不耦合**固体块的并行 GS 推进回归算例（P3）

用最简单的拓扑把「多进程各持有多个 `BLOCK_SOLID` 时，固体 GS 扫掠必须逐进程同步」这件事
变成可判定结果，避免界面物理（码 14 局部耦合仍存在的既有问题）盖住并行逻辑。

```
rank 0 : blk-1 (41×21×3 节点, k=16.3) + blk-2 (21×21×3, k=200)
rank 1 : blk-3 (61×21×3 节点, k=16.3) + blk-4 (15×21×3, k=1)      <- partation_np2.dat
```

每个块都是**一维（x 方向）导热**：`i-` 面 400 K 等温、`i+` 面 350 K 等温，
`j±/k±` 对称 ⇒ 精确解是线性剖面

$$T(x) = 400 - 50\,(x-x_0)\quad(x_0=0\ \text{blk1/2},\ x_0=2\ \text{blk3/4}),$$

与 $k$、$y$、$z$ 无关，且均匀网格上 FVM 对线性场**精确**（误差只来自 GS 未完全收敛）。
四个块的 x 分辨率互不相同 ⇒ 各块所需扫掠数差别很大（实测 5300 / 18100 / 18700 / 5200），
正是「各进程若各自早退就会失配」的场景。

## 判据

| 检查 | 期望 |
|---|---|
| np=1 vs np=2（`Solid_Sync_Chunk=100`）| `flow3d.dat`/`Residual.dat`/`Step_mess.dat` **逐字节相同**；`Ts_block_*.dat` 差 ≲1e-6 K（见下）|
| np=1 vs np=2（`Solid_Sync_Chunk=0`，旧路径）| `Ts_block_*.dat` **逐字节相同**（块间无耦合 ⇒ 逐块早退等价）|
| 与精确线性解 | $\max|T-T_{exact}|\lesssim10^{-6}$ K |
| 并行同步行 | np=2 每个时间步恰好一行 `solid sync (parallel): sweeps= ... (identical on every MPI rank)` |

实测（2 个时间步，`Solid_Tol=1e-9`）：

| 量 | np=1 | np=2 (`Chunk=100`) | np=2 (`Chunk=0`) |
|---|---|---|---|
| 各块 sweep 数 | 18100 / 5300 / 18700 / 5200 | **全部 18700**（全局最大） | 18100 / 5300 / 18700 / 5200 |
| `Ts` vs np=1 | — | $\max|\Delta T|=3.9\times10^{-7}$ K（blk3 逐字节相同） | **逐字节相同** |
| `Ts` vs 精确解 | $\le9.0\times10^{-7}$ K | $\le9.0\times10^{-7}$ K | $\le9.0\times10^{-7}$ K |

`Chunk>0` 的差异只来自**退出时机**（全局最大扫掠数 vs 各块自己的扫掠数），量级就是
`Solid_Tol`；`Chunk=0` 逐字节相同则说明旧路径在多进程下未被改动。

> **注意**：本算例的块之间没有界面，所以 `Chunk=0` 也不会挂死（跨进程无 $T_s$ 消息）。
> 历史上真正的死锁场景是「$T_s$ 跨进程交换 + 各进程扫掠数不等」，而码 14 的
> 跨进程耦合目前被 `check_remote_coupling_support` 直接拒绝，因此这里只能验证
> 「同步推进 → 扫掠数一致 / 结果与块→进程映射无关」，无法复现死锁本身。

## 文件

| 文件 | 作用 |
|---|---|
| `gen_case.py` | 生成 `Mesh3d.x` / `bc3d.inp` / `bc3d_interface.inp` / `material.in` / `solid_bc.inp` / `partation_np2.dat` |
| `control.ec` | 计算控制（`$solid_ec` 里 `Solid_Sync_Chunk=100`；改 0 可复现旧路径）|
| `prep.sh` | 建立 `np1/`（不放分区文件）与 `np2/`（`partation.dat`）两个运行目录（输入软链）|
| `analyze_solid_4blk.py` | 逐块对比两个运行目录的 `Ts_block_*.dat` 并给出与精确线性解的最大误差 |

## 运行

```bash
python3 gen_case.py
./prep.sh
( cd np1 && ../../../bin/MultiPhysMixSolver                > run.log 2>&1 )
( cd np2 && mpirun -n 2 ../../../bin/MultiPhysMixSolver    > run.log 2>&1 )
python3 analyze_solid_4blk.py np2 np1
```

`Ts_block_N.dat` 的编号是**全局块号**（2026-10-18 修复：此前用局部块号，np>1 时两个进程
会互相覆盖同名文件）。运行产物（`*.vtk`、`field_restart.dat`、`run.log`）按 `.gitignore` 不入库。
