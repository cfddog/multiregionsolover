#!/usr/bin/env python3
"""analyze_solid_4blk.py —— 4 块（互不耦合）固体导热算例的后处理。

用法
----
  python3 analyze_solid_4blk.py            # 默认比较 np2 与 np1
  python3 analyze_solid_4blk.py DIR [DIR2 ..]

对每个运行目录、每个块输出：
  * 单元数；
  * 与精确线性解 T(x)=400-50*(x-x0) 的最大误差（x0=0 给 blk1/2，x0=2 给 blk3/4）；
比较部分（给了两个以上目录时）：其余目录与**第一个**目录逐点 max|ΔT|、是否逐字节相同。
任一块对精确解误差 >1e-4 K、或目录间 max|ΔT| >1e-4 K ⇒ 退出码 1（FAIL）。

块几何/物性与 cases/solid_4blk/gen_case.py 一致：
  blk1 41x21x3  x=[0,1]  k=16.3      blk2 21x21x3  x=[0,1]  k=200
  blk3 61x21x3  x=[2,3]  k=16.3      blk4 15x21x3  x=[2,3]  k=1
"""
import os
import sys

X0 = {1: 0.0, 2: 0.0, 3: 2.0, 4: 2.0}       # block -> x origin of the linear profile
HOT, COLD, LEN = 400.0, 350.0, 1.0          # i- face, i+ face, block length


def t_exact(block, x):
    return HOT + (COLD - HOT) * (x - X0[block]) / LEN


def read_ts(path):
    rows = []
    with open(path) as f:
        for ln in f:
            s = ln.split()
            if len(s) != 4:
                continue
            try:
                rows.append(tuple(float(v) for v in s))
            except ValueError:
                pass
    return rows


def main(argv):
    dirs = [a for a in argv if not a.startswith("--")] or ["np2", "np1"]
    data = {}
    ok = True
    for d in dirs:
        data[d] = {}
        print("== %s" % d)
        for m in sorted(X0):
            fn = os.path.join(d, "Ts_block_%d.dat" % m)
            if not os.path.isfile(fn):
                print("  blk%d  MISSING %s" % (m, fn))
                ok = False
                continue
            r = read_ts(fn)
            data[d][m] = r
            err = max(abs(t - t_exact(m, x)) for (x, _, _, t) in r)
            line = "  blk%d  ncell=%4d  max|T-T_exact|=%.3e K" % (m, len(r), err)
            if err > 1.0e-4:
                line += "   <-- exceeds 1e-4 K"
                ok = False
            print(line)

    if len(dirs) > 1:
        ref = dirs[0]
        print("== comparison vs %s" % ref)
        for d in dirs[1:]:
            for m in sorted(X0):
                a, b = data.get(ref, {}).get(m), data.get(d, {}).get(m)
                if a is None or b is None:
                    print("  blk%d  MISSING data" % m)
                    ok = False
                    continue
                if len(a) != len(b):
                    print("  blk%d  cell count differs (%d vs %d)" % (m, len(a), len(b)))
                    ok = False
                    continue
                dt = max(abs(p[3] - q[3]) for p, q in zip(a, b))
                print("  %s vs %s  blk%d  max|dT|=%.3e K%s"
                      % (d, ref, m, dt, "  (identical)" if dt == 0.0 else ""))
                if dt > 1.0e-4:
                    print("      <-- exceeds 1e-4 K (projection/tolerance level is ~1e-6)")
                    ok = False
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
