#!/usr/bin/env python3
"""gen_control.py — 由算例自带 control.ec 生成本次运行的 control.ec（验证框架用）。

做法（避免出现两个同名 namelist 组）
------------------------------------
1. base 里**已有**的键：就地替换值（保留缩进，键名按传入写法，大小写不敏感匹配）。
2. base 里**没有**的键：按前缀归组，追加成一个新的 namelist 组：
     AC_*            -> $ac_ec
     LS_*            -> $lowspeed_ec
     Solid_*         -> $solid_ec
     Porous_*        -> $porous_ec
     Iflag_Couple*, Kstep_Couple*, Niter_Couple*, Tol_Couple*, Porous_Chunk* -> $couple_ec
     其余            -> $flow_ec
   若 base 里已经存在该组（例如已经是分组写法的 control.ec），**直接报错**，
   而不是追加出第二个同名组（Fortran namelist 只读第一个，会静默忽略）。
3. --steps 会默认写 `t_end=<steps>.d0`、`Kstep_save=<steps>`、`Kstep_show=<steps>`；
   默认再补 `Iflag_flow_node=1`、`Iflag_restart=-1`（关重启读写，重复跑总是从头开始）。
   显式 --set 的同名键优先。

用法
----
  python3 gen_control.py --base ../../channel_ac_2blk/control.ec \
      --out runs/t01_two_ac/np2/control.ec --steps 2000 --set AC_MPI_Sync=1
"""
import argparse
import os
import re
import sys

# 集群上 LANG=C 时 stdout 是 ascii，打印中文会抛 UnicodeEncodeError
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

GROUP_RULES = [
    ("ac_ec", r"^AC_"),
    ("lowspeed_ec", r"^LS_"),
    ("solid_ec", r"^Solid_"),
    ("porous_ec", r"^Porous_"),
    ("couple_ec", r"^(Iflag_Couple|Kstep_Couple|Niter_Couple|Tol_Couple|Porous_Chunk)"),
]
GROUP_ORDER = ["flow_ec", "ac_ec", "lowspeed_ec", "solid_ec", "porous_ec", "couple_ec"]


def group_of(key):
    for g, pat in GROUP_RULES:
        if re.match(pat, key, re.I):
            return g
    return "flow_ec"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gen_control.py",
                                 description="算例 control.ec -> 本次运行 control.ec")
    ap.add_argument("--base", required=True, help="算例自带的 control.ec")
    ap.add_argument("--out", required=True, help="要写出的 control.ec")
    ap.add_argument("--steps", type=int, help="外步数（写 t_end/Kstep_save/Kstep_show）")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL",
                    help="额外设置，可重复")
    ap.add_argument("--no-node-field", action="store_true",
                    help="不写 Iflag_flow_node=1（默认写，供 analyze_para.py 比较）")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args(argv)

    if not os.path.isfile(a.base):
        raise SystemExit("gen_control: 找不到 base control.ec: %s" % a.base)
    with open(a.base) as f:
        lines = f.read().splitlines()

    # ---- 目标键值（保持给定顺序；同名后者覆盖前者）
    want = []
    if a.steps:
        want += [("t_end", "%d.d0" % a.steps),
                 ("Kstep_save", "%d" % a.steps),
                 ("Kstep_show", "%d" % a.steps)]
    if not a.no_node_field:
        want += [("Iflag_flow_node", "1")]
    want += [("Iflag_restart", "-1")]
    for item in a.set:
        if "=" not in item:
            raise SystemExit("gen_control: --set 需要 KEY=VAL，收到 %r" % item)
        k, v = item.split("=", 1)
        want.append((k.strip(), v.strip()))
    seen = {}
    for k, v in want:
        seen[k.lower()] = (k, v)
    want = list(seen.values())

    # ---- 已存在的键：就地替换（全部出现处，保持一致）
    inplace, extra = [], []
    for k, v in want:
        pat = re.compile(r"^(\s*)" + re.escape(k) + r"(\s*)=", re.I)
        hit = False
        for n, ln in enumerate(lines):
            mt = pat.match(ln)
            if mt:
                lines[n] = "%s%s=%s" % (mt.group(1), k, v)
                hit = True
        (inplace if hit else extra).append((k, v))

    # ---- 追加分组（同名组冲突时报错）
    body = "\n".join(lines)
    groups = {}
    for k, v in extra:
        g = group_of(k)
        if re.search(r"^\s*\$\s*" + g + r"\b", body, re.I | re.M):
            raise SystemExit(
                "gen_control: base %s 里已经有 $%s 组，但键 %s 不在其中；\n"
                "  请把 %s=%s 直接写进 base 的 $%s 组，或从 --set 里去掉。" %
                (a.base, g, k, k, v, g))
        groups.setdefault(g, []).append((k, v))

    out = lines
    for g in GROUP_ORDER:
        if g not in groups:
            continue
        out.append("$%s" % g)
        for k, v in groups[g]:
            out.append("  %s=%s" % (k, v))
        out.append("$end")

    d = os.path.dirname(os.path.abspath(a.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(a.out, "w") as f:
        f.write("\n".join(out) + "\n")

    if not a.quiet:
        print("gen_control: %s -> %s" % (a.base, a.out))
        print("             就地替换 %d 个键: %s"
              % (len(inplace), ", ".join(k for k, _ in inplace) or "-"))
        for g in GROUP_ORDER:
            if g in groups:
                print("             追加 $%s: %s"
                      % (g, ", ".join("%s=%s" % kv for kv in groups[g])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
