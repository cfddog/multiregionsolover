# cluster_run_20261019 — 集群现场控制参数留存

`cases/solid_2blk_chain` 在集群 wuzh02 上验证“码 14 固--固耦合修复”时的现场文件：

- 作业：`sbatch -n 1 ../run.slurm`（作业号 46903753，墙钟 1 m 03 s，gfortran 7.3.1）
- 准备：`./prep.sh 2000` ⇒ 下面 `control.ec` 里 `Kstep_save=2000`、`t_end=2000.d0`、
  `Kstep_show=2000`（与仓库里的 `../control.ec` 短跑版本不同，`prep.sh <步数>` 会改写）
- 文件：
  - `control.ec` / `output_para.out`：`np1/` 下实际使用的控制参数与回显
  - `Residual.dat`：2000 步残差（固体 GS 收敛判据走 `|dTs|`，流体栏为 0）
  - `bc3d_interface.inc`：转换结果，两块 `i+/i-` 对接面 `L1 L2 L3 = 1 2 3`
- 分析结果见 `../analyze_np1.txt`（对精确解 max 误差 3.7e−9 K，界面跳变 0.625 K = 解析 ⇒ PASS）
