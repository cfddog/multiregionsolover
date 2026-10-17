# lidcavity_ac_3blk —— 3 块共形界面的 AC 方腔（跨进程界面通量回归算例）

把 `cases/lidcavity_ac/re100` 的二维顶盖驱动方腔（Re=100、42×42 单元、U_lid=1）
沿 x **等分成 3 块**（每块 14×42 单元、节点 15×43×2），两条共形 `BC_Inner(-1)`
界面位于 x=1/3、x=2/3。三块同为 `BLOCK_LOWSPEED`（`LS_Algorithm=3`，AC 人工压缩），
**全封闭域**（无压力边界，四壁无滑移），是验证「跨进程界面通量符号」的最小算例。

## 文件

| 文件 | 作用 |
|---|---|
| `gen_case.py` | 生成 `Mesh3d.x` / `bc3d.inp` / `material.in`（三块等分，界面用负码 -1） |
| `partation_np2.dat` | np=2 分区：{blk1,blk2} \| {blk3}（块数不等，专门覆盖进程负载不均） |
| `partation_np3.dat` | np=3 分区：1 / 1 / 1 |
| `prep.sh` | 建立 `np1/` `np2/` `np3/` 三个运行目录（输入用软链），np=1 不写分区文件 |
| `gen_single42.py` | 生成同网格（42×42）**单块**参考算例 `s42/`，用于「多块 vs 单块」对照 |
| `analyze_3blk.py` | 把 `flow3d_block_*.vtk` 按 +x 拼回 42×42，输出 Ghia 比较、界面跳变、与参考目录的 RMS |

运行产物（`flow3d*.vtk`、`field_restart.dat`、`run*.log`）按 `.gitignore` 不入库。

## 运行

```bash
./prep.sh                     # 生成网格 + np1/np2/np3 运行目录
cd np1 && ../../../bin/MultiPhysMixSolver                > run.log 2>&1
cd ../np2 && mpirun -np 2 ../../../bin/MultiPhysMixSolver > run.log 2>&1
cd ../np3 && mpirun -np 3 ../../../bin/MultiPhysMixSolver > run.log 2>&1
python3 ../analyze_3blk.py . ../np1        # 与本机 np1 结果逐点比较
```

**两个必需设置**（`control.ec` 已写好）：

- `AC_MPI_Sync=1`：每个伪时间步刷新 `U(1:5)` 与 `p` 的 ghost 层（闭腔且 np>1 时
  压力水平由全局 1 号块锚定，见 `src/sub_multi_region.f90::init_ac_parallel_flags`）；
- `AC_Max_Iter=1`：每个外步只走 1 个伪步，使 np=1 的「伪步刷新」与 np>1 的
  「外步末刷新」节奏一致（多块 AC 的同步调度，参见 `cases/channel_ac_2blk/README.md`）。

## 结果（2026-10-17，修复后）

10⁴ 外步（`t_end=10000`、`Kstep_save=10000`）：

| 对比 | max\|Δu\| | max\|Δv\| | max\|Δp\| |
|---|---|---|---|
| np=1 vs np=2 | 1.16e-6 | 2.4e-7 | 9.8e-7 |
| np=1 vs np=3 | 1.16e-6 | 2.4e-7 | 9.8e-7 |
| np=2 vs np=3 | 0（逐位相同） | 0 | 0 |

收敛解（3×10⁴ 步，`res_q≈3e-14`）：u∈[−0.233, 0.926]（物理界内）；
vs Ghia et al. (1982) RMS u=0.02095、v=0.04365（单块 42×42 为 0.02013 / 0.04387）。

同网格单块对照：`max|Δu|=9.3e-3`、rms=3.0e-3，全部集中在界面两侧 2~3 列，
正是「界面刻意一阶」的预期误差（h=1/42）；离界面 3 胞以外与单块**逐位一致**。

## 这个算例抓住的两个 bug（修复见 `src/sub_lowspeed_ac.f90::ac_face_flux_rus1_own`）

跨进程回退路径 `ac_face_flux_rus1_own` 里两侧用**同一个**面法向（i 面恒 +x），
因此必须一侧加、一侧减：

1. 两侧都写 `RAC += FLX*A` ⇒ Rusanov 耗散项在界面相互抵消，界面棋盘模态零阻尼；
2. `is_high` 判据写反（ghost 索引更大时被判成高侧）⇒ 中央通量两侧同向反向，
   界面变成抽吸源，u 冲到 1.68/1.83 后发 NaN。

正解：低侧（`iG>iO`、`jG>jO`、`kG>kO`，判据只与 own/ghost 索引关系有关，与块号、
进程号无关）写 `RAC -= FLX*A`，高侧写 `RAC += FLX*A`。手册 §4.5 与 `tab:ac3blk`
记录了同一内容。

## 集群运行（wuzh02，Slurm）

```bash
# 1) 本地同步整仓（规则见 .clinerules「本地 ↔ 集群同步」）
bash util/sync_cluster.sh push

# 2) 集群端：编译 + 生成三个运行目录 + 提交
ssh wuzh02
cd ~/sundong/MultiPhysMixSolver && make -j4
cd cases/lidcavity_ac_3blk && ./prep.sh 10000     # 10000 步 = 已验证基线；不带参数则软链 control.ec(1000 步)
cd np1 && sbatch -n 1 ../run.slurm
cd ../np2 && sbatch -n 2 ../run.slurm
cd ../np3 && sbatch -n 3 ../run.slurm

# 3) 结果比对（界面连续性 + np 无关性）
python3 ../analyze_3blk.py np1
python3 ../analyze_3blk.py np2 np1
python3 ../analyze_3blk.py np3 np1
```

`-n` 必须等于该目录 `partation.dat` 的进程数（np2 = `{1,2}|{3}`，np3 = `1/1/1`）；
日志 `np*/run_out.log`、`np*/run_err.log`。算例自带 `control.ec` 仍是 1000 步的冒烟配置，
长跑请用 `./prep.sh 10000`。

## 现场目录（不入库）

`np1/` `np2/` `np3/`（`prep.sh` 生成）、`s1/` `s42/` `m42/` `v1/` `v3/` `r1/` `r2/`
`q1_10k/` `q2_10k/` `q3_10k/` `q1/` `q2/` `q3/` `d1/` `d3/` `d3b/` `t1/` 均为本轮
实测目录，可随时删除后用 `prep.sh` 重建。
