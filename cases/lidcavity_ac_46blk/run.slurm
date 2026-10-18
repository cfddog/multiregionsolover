#!/bin/bash
# 4/6 块 AC 方腔（任务 7 强扩展 + 任务 8 anchoring）的 Slurm 作业脚本，集群 wuzh02。
#
# 用法：在运行目录里提交，-n 必须等于该目录 partation.dat 的进程数
#   cd m6a_np1 && sbatch -n 1 ../run.slurm
#   cd m6a_np2 && sbatch -n 2 ../run.slurm
#   cd m6a_np3 && sbatch -n 3 ../run.slurm
#   cd m6a_np6 && sbatch -n 6 ../run.slurm
#   cd m4a_np2 && sbatch -n 2 ../run.slurm   （np=1 目录不带 partation.dat）
# 提交前：cd .. && ./prep.sh 10000
#
#SBATCH -p wzacnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J ac46blk
#SBATCH -o run_out.log
#SBATCH -e run_err.log

set -e
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# 运行目录是 cases/lidcavity_ac_46blk/<run>/，二进制在仓库根 bin/
BIN="$(cd ../../.. && pwd)/bin/MultiPhysMixSolver"
[ -x "$BIN" ] || { echo "找不到可执行文件 $BIN（先在仓库根 make -j4）" >&2; exit 1; }
echo "job=${SLURM_JOB_ID:-local} host=$(hostname) np=${SLURM_NTASKS:-?} dir=$PWD bin=$BIN start=$(date)"
time mpirun -n "${SLURM_NTASKS:-1}" "$BIN"
echo "done $(date)"
