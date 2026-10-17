#!/bin/bash
# 3 块 AC 方腔的 Slurm 作业脚本（集群 wuzh02；模板源自 ~/sundong/run.slurm）。
#
# 用法：在运行目录里提交，-n 必须等于该目录 partation.dat 的进程数
#   cd np1 && sbatch -n 1 ../run.slurm
#   cd np2 && sbatch -n 2 ../run.slurm     # partation_np2: {1,2}|{3}
#   cd np3 && sbatch -n 3 ../run.slurm     # partation_np3: 1/1/1
# 提交前：cd .. && ./prep.sh 10000（生成 mesh/partation/control.ec）
#
#SBATCH -p wzacnormal
#SBATCH -N 1
#SBATCH -n 3
#SBATCH -J ac3blk
#SBATCH -o run_out.log
#SBATCH -e run_err.log

set -e
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# 运行目录是 cases/lidcavity_ac_3blk/np*/，二进制在仓库根 bin/
BIN="$(cd ../../.. && pwd)/bin/MultiPhysMixSolver"
[ -x "$BIN" ] || { echo "找不到可执行文件 $BIN（先在仓库根 make -j4）" >&2; exit 1; }
echo "job=${SLURM_JOB_ID:-local} host=$(hostname) np=${SLURM_NTASKS:-?} dir=$PWD bin=$BIN start=$(date)"
time mpirun -n "${SLURM_NTASKS:-1}" "$BIN"
echo "done $(date)"
