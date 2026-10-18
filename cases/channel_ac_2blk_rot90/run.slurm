#!/bin/bash
# i<->j（索引 90 度旋转）两块通道的 Slurm 作业脚本（集群 wuzh02）。
#
# 用法：在运行目录里提交，-n 必须等于该目录 partation.dat 的进程数
#   cd p_np1 && sbatch -n 1 ../run.slurm
#   cd p_np2 && sbatch -n 2 ../run.slurm
#   cd n_np1 && sbatch -n 1 ../run.slurm
#   cd n_np2 && sbatch -n 2 ../run.slurm
#   cd ref_np1 && sbatch -n 1 ../run.slurm
# 提交前：cd .. && ./prep.sh 40000
#
#SBATCH -p wzacnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J acrot90
#SBATCH -o run_out.log
#SBATCH -e run_err.log

set -e
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# 运行目录是 cases/channel_ac_2blk_rot90/<run>/，二进制在仓库根 bin/
BIN="$(cd ../../.. && pwd)/bin/MultiPhysMixSolver"
[ -x "$BIN" ] || { echo "找不到可执行文件 $BIN（先在仓库根 make -j4）" >&2; exit 1; }
echo "job=${SLURM_JOB_ID:-local} host=$(hostname) np=${SLURM_NTASKS:-?} dir=$PWD bin=$BIN start=$(date)"
time mpirun -n "${SLURM_NTASKS:-1}" "$BIN"
echo "done $(date)"
