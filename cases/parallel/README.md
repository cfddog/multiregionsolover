# cases/parallel —— np 无关性验证框架（并行改造的回归底座）

给 `memory-bank/progress.md` 里 **P3–P6**（固体块、跨类型界面的码 12 / 19、多区域耦合的
码 14 / 18、以及 P1 之后新增的交换路径）准备的**共用**验证工具：任何一次并行相关的改动，
只要求「`np=1` 与 `np>1` 跑到同一个解」，而不要求逐位相同。

本目录**不含算例**：每个测试都复用已有算例目录里已入库的网格/边界/物性文件（软链），
只额外生成 `control.ec`（改步数与并行开关）和 `partation.dat`（不改算例目录里的版本）。

## 文件

| 文件 | 作用 |
|---|---|
| `tests.conf` | 测试清单（9 列竖线分隔），见下 |
| `run_para_verify.sh` | 驱动：建运行目录 → 生成控制参数/分区 → 跑 np=1..k → 逐点比对 → 汇总 |
| `gen_control.py` | 由算例自带 `control.ec` 生成本次运行的 `control.ec`（就地替换已有键，缺少的键按前缀**追加**到 `$flow_ec/$ac_ec/...`，base 里已有同名组时**报错**而不是追加重复组） |
| `analyze_para.py` | 断言：`flow3d_node.dat` 逐点比较（d,u,v,w,T,Ts）+ `flow3d_block_<m>.vtk` 的压力（**去常数**） |
| `run.slurm` | Slurm 模板（`MODE=slurm` 时拷进每个运行目录，`sbatch --wait -n <np>`） |
| `runs/`（gitignore） | 生成的运行目录与结果 |

外部依赖：仓库根的 `util/mkpart.py`（读写 `partation.dat`，能读 Fortran unformatted 网格）。

## 测试清单 `tests.conf`

```
name | case_dir | mesh | bc | material | nblocks | np_list | steps | overrides
```

- 第 1 列是测试名，也是 `runs/<name>/` 目录名；命令行第一个参数按**前缀**过滤。
- `np_list` 必须含 `1`（参考解）；`overrides` 是 `KEY=VAL,KEY=VAL`，交给 `gen_control.py`。
- `nblocks` 目前只作说明（真正的块数来自网格，`mkpart.py --mesh` 会校验文件格式）。

当前的两个基线测试（都是 AC，`LS_Algorithm=3`，即 P0–P2 已验证的路径）：

| 测试 | 算例 | 块 | np | 覆盖的界面 |
|---|---|---|---|---|
| `t01_two_ac` | `cases/channel_ac_2blk`（进口/出口槽道） | 2 | 1,2 | 1 条 i 面，**有**压力边界 ⇒ 压力水平由 BC 固定 |
| `t02_four_ac` | `cases/lidcavity_ac_46blk`（4 块 m4a 方腔） | 4 | 1,2,3,4 | 2 条 i 面 + 1 条 j 面，**全封闭** ⇒ 压力靠全局 1 号块锚定 |

两者都强制 `AC_Max_Iter=1`：多块 AC 的 np 无关性配方（每个外步只推 1 个伪步，见
`cases/lidcavity_ac_3blk/README.md`）。`channel_ac_2blk` 自带的 `control.ec` 是
`AC_Max_Iter=20` 的**单进程**配方，如果不覆盖，`np=2` 会差 6e-3（实测，见下「为什么」）。

## 用法

```bash
# 本地冒烟（约 4 分钟；阈值必须放宽，理由见「判据」）
STEPS=2000 TOLVEL=1e-3 TOLRMS=1e-4 TOLP=1e-3 bash run_para_verify.sh

# 完整（默认 10000 步 + 紧阈值；本地约 15 分钟，建议放集群）
bash run_para_verify.sh

# 只跑某一个测试 / 只跑某些进程数
bash run_para_verify.sh t02
NP=1,2 bash run_para_verify.sh

# 集群（wuzh02）：每个运行目录单独 sbatch --wait，作业从不落在登录节点
MODE=slurm bash run_para_verify.sh          # 建议在集群侧执行
# 或本地驱动远端：
ssh wuzh02 'cd <repo>/cases/parallel && MODE=slurm bash run_para_verify.sh'
```

环境变量：`BIN`（默认 `<repo>/bin/MultiPhysMixSolver`）、`MPIRUN`、`STEPS`、`NP`、
`TOLVEL/TOLRMS/TOLP`、`WORK`（默认 `runs`）、`MODE=run|slurm`。

集群注意事项：
- 用共享的 `cases/parallel/` 目录里的日志看进度（`runs/para_verify.log`）——
  wuzh02 有多个登录节点（login03/login05/…），`/tmp` **不共享**。
- 提交作业前不要 `module load anaconda`（规则见 `memory-bank/constraints.md`）；
  `run.slurm` 与 `cases/lidcavity_ac_3blk/run.slurm` 一致，不需要额外 module。


## 判据（为什么不是「逐位相同」）

`np=1` 与 `np>1` 的界面推进方式本来就不同：

- `np=1`：所有界面都在同一进程，`sub_lowspeed_ac.f90` 在**每个伪步**内
  `update_buffer_onemesh` + `update_p_buffer_onemesh`（块 Gauss-Seidel）；
- `np>1`：伪步内不刷（MPI 调用次数会在各进程间错配 → 死锁），只在**每个外步末**统一刷新
  （Jacobi），界面通量两侧各算一阶 Rusanov。

所以两者**只可能收敛到同一个解**，瞬态轨迹必然不同。判据因此是：

| 量 | 容差 | 说明 |
|---|---|---|
| `u,v,w` 的 max\|Δ\| | 1e-5 | 绝对（这些算例 Uref=1） |
| `u,v,w` 的 RMS\|Δ\| | 5e-6 | 残余误差的整体水平 |
| `d` | 1e-5 | 等温不可压，应当≈0（实测恒 0） |
| 压力：去常数 max\|Δp−mean(Δp)\|/max\|p\| | 5e-6 | 封闭域的压力水平是 gauge（`init_ac_parallel_flags` 用全局 1 号块锚定），所以先减掉每块的均值 |
| `T,Ts` 的相对 max\|ΔT\|/max\|T\| | 1e-6 | 为 P3（固体温度）预留 |

`np>1` 之间是**逐位相同**的（实测 np=2/3/4 的 max\|Δ\| 恒为 0），所以这一层没有容差问题。

## 结果（2026-10-18，wuzh02，gfortran 7.3.1，10000 步）

| 测试 | np | 速度 max\|Δ\| | 速度 RMS\|Δ\| | 压力 去常数(rel) | 常数偏移 | 结论 |
|---|---|---|---|---|---|---|
| `t01_two_ac` | 2 | 4.03e-6 | 2.67e-6 | 2.31e-6 | 9.72e-6 | PASS |
| `t02_four_ac` | 2 | 1.37e-6 | 7.47e-7 | 1.22e-6 | 1.51e-6 | PASS |
| `t02_four_ac` | 3 | 1.37e-6 | 7.47e-7 | 1.22e-6 | 1.51e-6 | PASS |
| `t02_four_ac` | 4 | 1.37e-6 | 7.47e-7 | 1.22e-6 | 1.51e-6 | PASS |

- 集群机时：t01 np=1 75 s / np=2 43 s；t02 np=1 42 s / np=2 43 s / np=3 43 s / np=4 10 s。
- 短跑（2000 步）的差距大一个量级以上（t01 4.1e-5、t02 8.5e-4），**只用于冒烟**，
  不要拿它判断并行正确性。
- 反面教材：t01 若不覆盖 `AC_Max_Iter=1`（沿用算例自带的 20），np=2 会差 6.6e-3，
  说明本框架确实能抓到「界面耦合节奏不一致」这类问题。

现场记录（`pull` 回来的）：`runs/<test>/np<k>/control.ec`、`partation.dat`、
`output_para.out`、`analyze_np*.txt`。

## 给 P3–P6 加测试

`tests.conf` 加一行即可，不需要改脚本：

```conf
t03_solid_ac | cases/<固体/混合算例目录> | Mesh3d_m4a.x | bc3d_m4a.inp | material_m4.in | 4 | 1,2,4 | 10000 | AC_MPI_Sync=1,AC_Max_Iter=1
```

要额外注意：

1. **块数/分区**：每个 np 都由 `util/mkpart.py --mesh` 自动均分（`split = ceil(nb/np)`）；
   若某次改造要求特定分区（例如「流体块全在 proc 0」），在 `run_np` 里改用
   `mkpart.py --split/--blocks`。
2. **多区域耦合码（14/18）**：目前 `flow3d_node.dat` 只有 d,u,v,w,T,Ts。等 P3/P6 让固体
   位移/应力进入 `output_flow_node`，同步扩展 `analyze_para.py` 的 `VARNAMES` 与断言
   （脚本按 `NVAR=6` 校验记录长度，字段数变了会直接报错而不是静默错位）。
3. **跨类型界面（码 12/19）**：把 `bc3d_interface.inp` 也软链进运行目录（在 `run_np`
   里加一行），并确认 `Iflag_bc_check` 打印的界面计数与 np 无关。
4. 新测试先 `STEPS=2000` 冒烟（阈值放宽），确认 `np>1` 之间逐位相同后再跑 10000 步定基线；
   结果写进本文件「结果」表并同步 `memory-bank/worklog.md`。
5. **闭合域（无压力边界）**：AC 的压力水平只由**全局块号最小的低速/多孔块**锚定
   （`AC_Anchor_Block`，`init_ac_parallel_flags` 用一次 `MPI_MIN` 定好，启动行打印
   `pressure anchor block=`）。加多孔块测试时注意：`porous.inp` **第一行是条目数**，
   漏掉会让块号/参数整体错位（块会退回默认 `eps=0.9/dp=1e-3`，结果看着正常但物性不对）。
   闭合域 + 多块时压力只差一个常数是允许的，`analyze_para.py` 已按“去常数”比较。

## 已知限制

- 只覆盖 AC（`LS_Algorithm=3`）路径；P4/P5 的码 12/19、P6 的码 14/18 落地后必须补测试。
- 压力取自 ASCII VTK 的单元中心值，只用于**一致性**判断（不做守恒/精度分析）。
- `flow3d_node.dat` 需要 `Iflag_flow_node=1`（`gen_control.py` 默认写）；关掉它则只剩压力比对。
- `field_restart.dat` 不参与比对（`Iflag_restart=-1`，每次都从头跑）。
