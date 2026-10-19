#!/bin/bash
# 2 块固-固串联导热（bc3d.inp 用码 14 = solid-solid，连接描述只在
# bc3d_interface.inp）的 Slurm 作业脚本（集群 wuzh02）。
#
# 只能 np1：码 14 的跨进程固--固耦合尚未实现，check_remote_coupling_support
# 在 np>1 时会中止作业。
#   cd cases/solid_2blk_chain && ./prep.sh 2000   # 生成 np1/（步数可省）
#   cd np1 && sbatch -n 1 ../run.slurm
#
#SBATCH -p wzacnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J s2chain
#SBATCH -o run_out.log
#SBATCH -e run_err.log

set -e
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# 运行目录是 cases/solid_2blk_chain/np1/，二进制在仓库根 bin/
BIN="$(cd ../../.. && pwd)/bin/MultiPhysMixSolver"
[ -x "$BIN" ] || { echo "找不到可执行文件 $BIN（先在仓库根 make -j4）" >&2; exit 1; }
echo "job=${SLURM_JOB_ID:-local} host=$(hostname) np=${SLURM_NTASKS:-?} dir=$PWD bin=$BIN start=$(date)"
time mpirun -n "${SLURM_NTASKS:-1}" "$BIN"
echo "done $(date)"
