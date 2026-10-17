#!/bin/bash
# 建立 np1/ np2/ np3/ 三个运行目录（共享输入文件用软链），
# 并按需放置 partation.dat：np=1 不放（自动分区=全部块在 0 号进程），
# np=2 用 {1,2}|{3}，np=3 用 1/1/1。
#
# 用法：./prep.sh [总步数]
#   不给参数：三个目录软链算例自带 control.ec（t_end=1000，快速冒烟）；
#   给步数（推荐 10000）：各目录写**独立** control.ec，t_end=Kstep_save=步数、
#   Kstep_show=步数/5 —— 与已验证的 1e4 步基线（q*_10k）完全一致。
set -e
cd "$(dirname "$0")"
python3 gen_case.py
NSTEP=${1:-}
for np in 1 2 3; do
  mkdir -p np$np
  for f in Mesh3d.x bc3d.inp material.in; do
    ln -sf ../$f np$np/$f
  done
  # 先删旧的（可能是软链）：直接 `> 软链` 会截断目标文件（2026-10-18 踩过）
  rm -f np$np/control.ec
  if [ -n "$NSTEP" ]; then
    sed -e "s/^  t_end=.*/  t_end=${NSTEP}.d0/" \
        -e "s/^  Kstep_save=.*/  Kstep_save=${NSTEP}/" \
        -e "s/^  Kstep_show=.*/  Kstep_show=$((NSTEP / 5))/" \
        control.ec > np$np/control.ec.tmp
    mv np$np/control.ec.tmp np$np/control.ec
  else
    ln -sf ../control.ec np$np/control.ec
  fi
  rm -f np$np/partation.dat
done
cp partation_np2.dat np2/partation.dat
cp partation_np3.dat np3/partation.dat
echo "np1/ np2/ np3/ ready（步数 ${NSTEP:-算例自带 control.ec}）"

