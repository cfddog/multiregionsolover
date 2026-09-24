#!/usr/bin/env python3
"""check_flow3d_node.py — 校验 flow3d_node.dat（节点值、SI 单位）与 Mesh3d.x 是否配套。

用法：python3 check_flow3d_node.py [算例目录]        （默认当前目录）

flow3d_node.dat 由 src/sub_IO.f90::output_flow_node 写出（control.ec 的
$flow_ec: Iflag_flow_node=1）：unformatted、逐块一条记录、块序与 Mesh3d.x 相同、
每条记录含 6 个节点场： d, u, v, w, T, Ts
  T  = 流体温度 [K]（固体块为固体温度 Ts，兼容 flow3d.vtk 约定）
  Ts = 固体/骨架温度 [K]（固体与多孔为真实值；流体/低速块无固相，Ts=T）
"""
import struct
import sys
import os

d = (sys.argv[1] if len(sys.argv) > 1 else ".").rstrip("/")
mesh = os.path.join(d, "Mesh3d.x")
nodef = os.path.join(d, "flow3d_node.dat")
NVAR = 6


def read_plot3d_grid_dims(path):
    """返回 (nblock, [(ni,nj,nk), ...])（只读前两条记录）。"""
    with open(path, "rb") as f:
        rc = struct.unpack("<I", f.read(4))[0]
        nb = struct.unpack("<i", f.read(rc))[0]
        f.read(4)
        rc = struct.unpack("<I", f.read(4))[0]
        dims = struct.unpack("<%di" % (rc // 4), f.read(rc))
        f.read(4)
    return nb, [tuple(dims[3 * k:3 * k + 3]) for k in range(nb)]


nb, dims = read_plot3d_grid_dims(mesh)
print("Mesh3d.x    : nblock=%d  dims=%s" % (nb, dims))

recs = []
with open(nodef, "rb") as f:
    while True:
        raw = f.read(4)
        if len(raw) < 4:
            break
        rc = struct.unpack("<I", raw)[0]
        payload = f.read(rc)
        f.read(4)
        recs.append(struct.unpack("<%dd" % (rc // 8), payload))

print("flow3d_node : %d records, %d vars/record" % (len(recs), NVAR))
ok = len(recs) == nb
for m, vals in enumerate(recs):
    ni, nj, nk = dims[m]
    n = ni * nj * nk
    if len(vals) != NVAR * n:
        ok = False
        print("  block %d: SIZE MISMATCH got %d expect %d" % (m + 1, len(vals), NVAR * n))
        continue
    dd = vals[0 * n:1 * n]
    uu = vals[1 * n:2 * n]
    tt = vals[4 * n:5 * n]
    ts = vals[5 * n:6 * n]
    print("  block %d: npts=%-6d (=%dx%dx%d)  d[%.4g,%.4g] |u|[%.4g,%.4g]  "
          "T[%.4g,%.4g] K  Ts[%.4g,%.4g] K  max|T-Ts|=%.4g K"
          % (m + 1, n, ni, nj, nk, min(dd), max(dd),
             min(min(uu), min(vals[2 * n:3 * n]), min(vals[3 * n:4 * n])),
             max(max(uu), max(vals[2 * n:3 * n]), max(vals[3 * n:4 * n])),
             min(tt), max(tt), min(ts), max(ts),
             max(abs(a - b) for a, b in zip(tt, ts))))
print("STRUCTURE OK" if ok else "STRUCTURE MISMATCH")
