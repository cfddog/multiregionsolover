#!/bin/bash
# 建立 np1/ 运行目录（输入文件用软链），可选覆盖总步数。
#
# 用法：./prep.sh [总步数]
#   不给参数：软链算例自带 control.ec（t_end=3）；
#   给步数  ：写独立 control.ec，t_end=Kstep_save=Kstep_show=步数。
#
# 注：码 14（固--固）跨进程耦合尚未实现（check_remote_coupling_support 会在
# np>1 时中止作业），因此本算例只准备 np1/。
set -e
cd "$(dirname "$0")"
python3 gen_case.py
NSTEP=${1:-}
mkdir -p np1
for f in Mesh3d.x bc3d.inp bc3d_interface.inp material.in solid_bc.inp; do
  ln -sf ../$f np1/$f
done
# 先删旧的（可能是软链）：直接 `> 软链` 会截断目标文件
rm -f np1/control.ec
if [ -n "$NSTEP" ]; then
  sed -e "s/^  t_end=.*/  t_end=${NSTEP}.d0/" \
      -e "s/^  Kstep_save=.*/  Kstep_save=${NSTEP}/" \
      -e "s/^  Kstep_show=.*/  Kstep_show=${NSTEP}/" \
      control.ec > np1/control.ec.tmp
  mv np1/control.ec.tmp np1/control.ec
else
  ln -sf ../control.ec np1/control.ec
fi
echo "np1/ ready（步数 ${NSTEP:-算例自带 control.ec}）"
echo "  cd np1 && ../../../bin/MultiPhysMixSolver"
echo "  python3 ../analyze_solid_2blk_chain.py np1"
