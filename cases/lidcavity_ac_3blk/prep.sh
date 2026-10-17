#!/bin/bash
# 建立 np1/ np2/ np3/ 三个运行目录（共享输入文件用软链），
# 并按需放置 partation.dat：np=1 不放（自动分区=全部块在 0 号进程），
# np=2 用 {1,2}|{3}，np=3 用 1/1/1。
set -e
cd "$(dirname "$0")"
python3 gen_case.py
for np in 1 2 3; do
  mkdir -p np$np
  for f in Mesh3d.x bc3d.inp material.in control.ec; do
    ln -sf ../$f np$np/$f
  done
  rm -f np$np/partation.dat
done
cp partation_np2.dat np2/partation.dat
cp partation_np3.dat np3/partation.dat
echo "np1/ np2/ np3/ ready"
