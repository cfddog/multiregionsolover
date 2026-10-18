#!/bin/bash
# 建立 np1/ np2/ 两个运行目录（共享输入用软链），并按需放置 partation.dat：
#   np=1 不放（程序自动分区 = 全部块在 0 号进程）
#   np=2 用 {blk1,blk2} | {blk3,blk4}（每进程两个固体块 —— P3 的最小并行固体拓扑）
#
# 用法：./prep.sh [总步数]
#   不给参数：软链算例自带 control.ec（t_end=4）；
#   给步数（推荐 2）：写**独立** control.ec，t_end=Kstep_save=步数、Kstep_show=步数。
set -e
cd "$(dirname "$0")"
python3 gen_case.py
NSTEP=${1:-}
for np in 1 2; do
  mkdir -p np$np
  for f in Mesh3d.x bc3d.inp bc3d_interface.inp material.in solid_bc.inp; do
    ln -sf ../$f np$np/$f
  done
  # 先删旧的（可能是软链）：直接 `> 软链` 会截断目标文件（2026-10-18 踩过）
  rm -f np$np/control.ec
  if [ -n "$NSTEP" ]; then
    sed -e "s/^  t_end=.*/  t_end=${NSTEP}.d0/" \
        -e "s/^  Kstep_save=.*/  Kstep_save=${NSTEP}/" \
        -e "s/^  Kstep_show=.*/  Kstep_show=${NSTEP}/" \
        control.ec > np$np/control.ec.tmp
    mv np$np/control.ec.tmp np$np/control.ec
  else
    ln -sf ../control.ec np$np/control.ec
  fi
  rm -f np$np/partation.dat
done
cp partation_np2.dat np2/partation.dat
echo "np1/ np2/ ready（步数 ${NSTEP:-算例自带 control.ec}）"
echo "  np1: ../../../bin/MultiPhysMixSolver"
echo "  np2: mpirun -n 2 ../../../bin/MultiPhysMixSolver"
echo "旧路径对照：在 np2/control.ec 的 \$solid_ec 里把 Solid_Sync_Chunk 改成 0"
