#!/bin/bash
# cases/parallel/run_para_verify.sh —— np 无关性验证框架（并行改造 P3–P6 的共用底座）
#
# 对 tests.conf 里的每个测试：为 np=1,2,... 各建一个运行目录（输入用软链），
# 生成 control.ec（gen_control.py）与 partation.dat（util/mkpart.py），跑完用
# analyze_para.py 把 np>1 的结果与 np=1 参考逐点比较并断言。
#
# 用法：
#   bash run_para_verify.sh                 # 跑 tests.conf 全部测试
#   bash run_para_verify.sh t01             # 只跑名字以 t01 开头的
#   STEPS=500 NP=1,2 bash run_para_verify.sh t02   # 覆盖步数/进程数（冒烟）
#   MODE=slurm bash run_para_verify.sh      # 集群：每个目录 sbatch --wait，
#                                           #   BIN 通过 --export=ALL 传给 run.slurm
#
# 环境变量：
#   BIN     可执行文件（默认 <repo>/bin/MultiPhysMixSolver）
#   MPIRUN  默认 mpirun
#   STEPS   覆盖 tests.conf 的步数
#   NP      只跑这些进程数（含 1，逗号分隔）
#   TOLVEL / TOLRMS / TOLP / TOLT  analyze_para.py 的阈值
#                           （默认 1e-5 / 5e-6 / 5e-6 / 1e-6，是 2026-10-18 集群实测
#                            1e4 步结果 4.0e-6/2.7e-6/2.3e-6 的约 2 倍余量；
#                            2000 步冒烟请放宽到 1e-3 / 1e-4 / 1e-3 / 1e-3）
#   WORK    运行目录根（默认 runs）
#   SKIP_RUN=1  不重跑，只对 runs/<test>/np* 里已有结果做比对（改阈值后快速复算）
# 退出码：0 = 全部 PASS，1 = 有 FAIL/运行失败，2 = 环境不满足。
set -uo pipefail
cd "$(dirname "$0")"
HERE=$PWD
ROOT=$(cd ../.. && pwd)

CONF=${CONF:-tests.conf}
WORK=${WORK:-runs}
BIN=${BIN:-$ROOT/bin/MultiPhysMixSolver}
MPIRUN=${MPIRUN:-mpirun}
MODE=${MODE:-run}
STEPS_OVR=${STEPS:-}
NP_OVR=${NP:-}
TOLVEL=${TOLVEL:-1e-5}
TOLRMS=${TOLRMS:-5e-6}
TOLP=${TOLP:-5e-6}
TOLT=${TOLT:-1e-6}
FILTER=${1:-}
export BIN
# 集群 LANG=C 时 python3 的 stdout 是 ascii，脚本里的中文提示会抛异常
export PYTHONIOENCODING=utf-8

trim() { sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'; }
say()  { printf '%s\n' "$*"; }

[ -f "$CONF" ] || { say "找不到 $CONF"; exit 2; }
[ -x "$BIN" ]  || { say "找不到可执行文件 $BIN（先在仓库根 make -j4，或用 BIN=... 指定）"; exit 2; }
command -v python3 >/dev/null || { say "需要 python3"; exit 2; }
if [ "$MODE" != slurm ] && ! command -v "$MPIRUN" >/dev/null; then
  say "找不到 $MPIRUN（np>1 需要 MPI）；或设 MODE=slurm 用 sbatch"; exit 2
fi

say "== 并行一致性验证 =="
say "   root : $ROOT"
say "   bin  : $BIN"
say "   mode : $MODE     tests: ${FILTER:-<all>}    work: $WORK"

# ---------------------------------------------------------------- 单次运行
# run_np <name> <np> <steps> <overrides> <mesh> <bc> <material> <casedir>
run_np() {
  local name=$1 np=$2 steps=$3 ovr=$4 mesh=$5 bc=$6 mat=$7 casedir=$8
  local td=$WORK/$name rd=$WORK/$name/np$np
  mkdir -p "$td"
  rm -rf "$rd"; mkdir -p "$rd"
  ln -s "$ROOT/$casedir/$mesh" "$rd/Mesh3d.x"
  ln -s "$ROOT/$casedir/$bc"   "$rd/bc3d.inp"
  ln -s "$ROOT/$casedir/$mat"  "$rd/material.in"

  local args=(--base "$ROOT/$casedir/control.ec" --out "$rd/control.ec"
              --steps "$steps" --quiet)
  if [ -n "$ovr" ]; then
    local kv
    IFS=',' read -r -a _kvs <<< "$ovr"
    for kv in "${_kvs[@]}"; do
      kv=$(printf '%s' "$kv" | trim)
      [ -n "$kv" ] && args+=(--set "$kv")
    done
  fi
  python3 gen_control.py "${args[@]}" < /dev/null || {
    say "  [fail] $name np=$np: gen_control 失败"; return 1; }
  if [ "$np" -gt 1 ]; then
    python3 "$ROOT/util/mkpart.py" --mesh "$rd/Mesh3d.x" --np "$np" \
            -o "$rd/partation.dat" --quiet < /dev/null || {
      say "  [fail] $name np=$np: mkpart 失败"; return 1; }
  fi

  local rc T0 T1
  T0=$(date +%s)
  if [ "$MODE" = slurm ]; then
    cp run.slurm "$rd/run.slurm"
    ( cd "$rd" && sbatch --wait --export=ALL -n "$np" \
        --job-name="pv_${name}_np$np" run.slurm > sbatch.out 2>&1 < /dev/null )
    rc=$?
  elif [ "$np" -eq 1 ]; then
    ( cd "$rd" && "$BIN" > run.log 2>&1 < /dev/null )
    rc=$?
  else
    ( cd "$rd" && "$MPIRUN" -n "$np" "$BIN" > run.log 2>&1 < /dev/null )
    rc=$?
  fi
  T1=$(date +%s)
  if [ $rc -ne 0 ]; then
    say "  [fail] $name np=$np: 退出码 $rc（$((T1 - T0)) s）"
    if [ -f "$rd/run.log" ]; then
      say "         run.log 末尾:"; tail -n 12 "$rd/run.log" | sed 's/^/         /'
    fi
    if [ -f "$rd/run_err.log" ]; then
      say "         run_err.log 末尾:"; tail -n 12 "$rd/run_err.log" | sed 's/^/         /'
    fi
    return 1
  fi
  printf '  [run ] %-12s np=%-2s %5ss\n' "$name" "$np" "$((T1 - T0))"
  return 0
}

# ---------------------------------------------------------------- 主流程
# 配置先整体读进内存（mapfile）：循环体里的子进程（mpirun/求解器）会继承
# 循环的 stdin，若循环直接从文件读，子进程可能把剩余的配置行吃掉
# （2026-10-18 实测：t01 跑完后 t02 直接消失）。
declare -a CONF_LINES=()
mapfile -t CONF_LINES < "$CONF"
if [ ${#CONF_LINES[@]} -eq 0 ]; then say "$CONF 为空"; exit 2; fi

npass=0; nfail=0; nskip=0
declare -a SUM=()
for raw in "${CONF_LINES[@]}"; do
  line=$(printf '%s' "$raw" | trim)
  case "$line" in ''|\#*) continue;; esac
  IFS='|' read -r -a f <<< "$line"
  name=$(printf '%s' "${f[0]:-}" | trim);  [ -z "$name" ] && continue
  casedir=$(printf '%s' "${f[1]:-}" | trim)
  mesh=$(printf '%s'    "${f[2]:-}" | trim)
  bc=$(printf '%s'      "${f[3]:-}" | trim)
  mat=$(printf '%s'     "${f[4]:-}" | trim)
  nplist=$(printf '%s'  "${f[6]:-}" | trim)
  steps=$(printf '%s'   "${f[7]:-}" | trim)
  ovr=$(printf '%s'     "${f[8]:-}" | trim)

  if [ -n "$FILTER" ]; then
    case "$name" in "$FILTER"*) ;; *) continue;; esac
  fi
  [ -n "$STEPS_OVR" ] && steps=$STEPS_OVR
  [ -n "$NP_OVR" ] && nplist=$NP_OVR

  say ""
  say "-- $name  ($casedir, np=$nplist, $steps 步)"

  # 网格必须能读出块数（mkpart 的 --mesh 路径同时验证了文件格式）
  if ! python3 "$ROOT/util/mkpart.py" --mesh "$ROOT/$casedir/$mesh" \
         --np 1 --out /dev/null --quiet >/dev/null 2>&1; then
    say "  [fail] 读不出 $casedir/$mesh 的块数"
    nfail=$((nfail + 1)); SUM+=("$name|np=?|FAIL(mesh)"); continue
  fi

  ok=1
  if [ "${SKIP_RUN:-0}" = "1" ]; then
    say "  [skip-run] SKIP_RUN=1：不重跑，只比对已有结果"
  else
    for np in ${nplist//,/ }; do
      run_np "$name" "$np" "$steps" "$ovr" "$mesh" "$bc" "$mat" "$casedir" || ok=0
    done
  fi

  ref=$WORK/$name/np1
  if [ ! -f "$ref/flow3d_node.dat" ]; then
    say "  [skip] $name: 没有 np=1 参考结果，跳过比对"
    nskip=$((nskip + 1)); SUM+=("$name|np=?|SKIP(no np1)"); continue
  fi
  for np in ${nplist//,/ }; do
    [ "$np" -eq 1 ] && continue
    rd=$WORK/$name/np$np
    out=$WORK/$name/analyze_np$np.txt
    if python3 analyze_para.py "$ref" "$rd" --tol-vel "$TOLVEL" --tol-rms "$TOLRMS" \
         --tol-p "$TOLP" --tol-T "$TOLT" --quiet > "$out" 2>&1; then
      npass=$((npass + 1)); SUM+=("$name|np=$np|PASS")
      say "  [pass] $name np=$np  $(grep -m1 '最坏值' "$out")"
    else
      nfail=$((nfail + 1)); SUM+=("$name|np=$np|FAIL")
      say "  [FAIL] $name np=$np  -> $out"
      tail -n 8 "$out" | sed 's/^/         /'
    fi
  done
  [ $ok -eq 0 ] && say "  [warn] $name 有运行失败（见上），比对结论不可信"
done

say ""
say "================ 汇总 ================"
printf '%-16s %-8s %s\n' "test" "np" "result"
for s in ${SUM[@]+"${SUM[@]}"}; do
  IFS='|' read -r a b c <<< "$s"
  printf '%-16s %-8s %s\n' "$a" "$b" "$c"
done
say "-------------------------------------"
say "PASS=$npass  FAIL=$nfail  SKIP=$nskip"
if [ $nfail -gt 0 ]; then
  say "结论：存在 FAIL（阈值/日志见 $WORK/<test>/analyze_np*.txt）"
  exit 1
fi
say "结论：全部通过（TOLVEL=$TOLVEL TOLRMS=$TOLRMS TOLP=$TOLP TOLT=$TOLT）"
exit 0
