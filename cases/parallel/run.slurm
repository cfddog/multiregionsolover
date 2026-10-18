#!/bin/bash
# cases/parallel/run.slurm —— 验证框架的 Slurm 作业脚本（集群 wuzh02）
#
# 由 run_para_verify.sh（MODE=slurm）拷进每个运行目录，并用
#   sbatch --wait --export=ALL -n <np> run.slurm
# 提交；进程数由命令行 -n 覆盖，BIN 由 --export=ALL 从调用者环境继承。
#
# 也可手工用：
#   cd runs/t01_two_ac/np2 && BIN=/path/to/bin/MultiPhysMixSolver sbatch -n 2 run.slurm
#
# 注意：提交作业前不要 `module load anaconda`（会让 mpirun/mpif90 混用，
# 见 memory-bank/constraints.md）。分析脚本可以单独用 anaconda 的 python3。
#SBATCH -p wzacnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J pverify
#SBATCH -o run_out.log
#SBATCH -e run_err.log

set -e
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# 运行目录是 cases/parallel/runs/<test>/np<k>/，仓库根在上 5 级
BIN="${BIN:-$(cd ../../../../.. && pwd)/bin/MultiPhysMixSolver}"
[ -x "$BIN" ] || { echo "找不到可执行文件 $BIN（先在仓库根 make -j4，或 export BIN=...）" >&2; exit 1; }
echo "job=${SLURM_JOB_ID:-local} host=$(hostname) np=${SLURM_NTASKS:-?} dir=$PWD bin=$BIN start=$(date)"
time mpirun -n "${SLURM_NTASKS:-1}" "$BIN"
echo "done $(date)"
