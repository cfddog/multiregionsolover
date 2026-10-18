#!/usr/bin/env python3
"""Generate a FOUR-BLOCK *decoupled* solid conduction case -> P3 regression.

Why this case exists
--------------------
It is the smallest topology that exercises the parallel (P3, 2026-10-18) solid
GS march with several SOLID blocks per process:

  rank 0 : blk-1 (41 x 21 x 3) + blk-2 (21 x 21 x 3)
  rank 1 : blk-3 (61 x 21 x 3) + blk-4 (15 x 21 x 3)      -> partation_np2.dat

The blocks are **not** coupled to each other (no code-14 interface): each one
is a plain 1-D-in-x conduction problem with two isothermal walls

    i- face : 400 K      i+ face : 350 K      j-/j+/k-/k+ : symmetry

so the exact steady state is the linear profile

    T(x) = 400 - 50 * (x - x0) / 1.0        (k, y and z independent)

on the uniform x grid of every block.  That gives three independent checks:

  1. the codes converge to an ANALYTIC solution (|T - T_exact| ~ 1e-9 K);
  2. np=1 (historical per-block march) and np=2 (chunked, Allreduce-stopped
     march) give the same Ts  ->  results are independent of the
     block->process mapping;
  3. no rank leaves the solid march early, so the per-sweep
     update_Ts_buffer_onemesh != deadlock.

Why decoupled instead of solid-solid interfaces
-----------------------------------------------
Cross-process code-14 (solid<->solid) coupling is refused by
check_remote_coupling_support, and the *local* code-14 path still shows a
pre-existing temperature jump across the interface (open issue, unrelated to
P3), which would mask an np-independence test.  With decoupled blocks the
discrete steady state is unique, so np=1 vs np=2 must agree to the GS
tolerance -- an unambiguous test of the sync logic itself.  The blocks have
different x resolutions and different k, so their local convergence histories
differ (which is exactly what the global max|dTs| stop has to reconcile).

Geometry (SI, 2-D in x-y, 3 nodes in z)
---------------------------------------
  blk-1: x in [0.0,1.0], y in [ 0.0, 0.1], k =  16.3
  blk-2: x in [0.0,1.0], y in [-0.2, 0.0], k = 200.0
  blk-3: x in [2.0,3.0], y in [ 0.0, 0.1], k =  16.3
  blk-4: x in [2.0,3.0], y in [-0.2, 0.0], k =   1.0

Writes Mesh3d.x, bc3d.inp, bc3d_interface.inp, material.in, solid_bc.inp and
partation_np2.dat.  Run with 'python3 gen_case.py' inside this directory.
"""
import struct, array

def record(f, data, fmt):
    if isinstance(data, array.array):
        raw = data.tobytes()
    else:
        raw = struct.pack(fmt, *data)
    f.write(struct.pack('I', len(raw))); f.write(raw); f.write(struct.pack('I', len(raw)))

def lin(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]

NY, NZ = 21, 3
ys_top = lin(0.0, 0.1, NY)          # blk-1 / blk-3 span
ys_bot = lin(-0.2, 0.0, NY)         # blk-2 / blk-4 span
zs = lin(0.0, 0.1, NZ)

# (x-range, y-range, conductivity, nx, name)  -- nx differs per block on purpose
BLOCKS = [((0.0, 1.0), ys_top, 16.3,  41, 'blk-1 (41 x 21 x 3, k=16.3)'),
          ((0.0, 1.0), ys_bot, 200.0, 21, 'blk-2 (21 x 21 x 3, k=200)'),
          ((2.0, 3.0), ys_top, 16.3,  61, 'blk-3 (61 x 21 x 3, k=16.3)'),
          ((2.0, 3.0), ys_bot, 1.0,   15, 'blk-4 (15 x 21 x 3, k=1.0)')]

def coords(xs_, ys_, zs_):
    xc, yc, zc = array.array('d'), array.array('d'), array.array('d')
    for k in range(len(zs_)):
        for j in range(len(ys_)):
            for i in range(len(xs_)):
                xc.append(xs_[i]); yc.append(ys_[j]); zc.append(zs_[k])
    return xc, yc, zc

with open("Mesh3d.x", "wb") as f:
    record(f, [len(BLOCKS)], 'i')
    record(f, [d for (_, _, _, NX, _) in BLOCKS for d in (NX, NY, NZ)],
           'i' * (3 * len(BLOCKS)))
    for (xr, ys_, _, NX, _) in BLOCKS:
        for a in coords(lin(xr[0], xr[1], NX), ys_, zs):
            record(f, a, '')

# ------------------------------------------------------------------ bc3d.inp
# face order is arbitrary; the reader derives i-/i+/j-/j+/k-/k+ from the ranges
def F(i1, i2, j1, j2, k1, k2, bc):
    return "%8d%8d%8d%8d%8d%8d%8d" % (i1, i2, j1, j2, k1, k2, bc)

BC_WALL, BC_SYM = 2, 3

def faces(NX):
    """6 subfaces of one solid block: isothermal i-/i+ walls, symmetry elsewhere."""
    return [F(1,  1, 1, NY, 1, NZ, BC_WALL),      # i-  -> 400 K (solid_bc)
            F(NX, NX, 1, NY, 1, NZ, BC_WALL),     # i+  -> 350 K (solid_bc)
            F(1, NX, 1, 1, 1, NZ, BC_SYM),        # j-
            F(1, NX, NY, NY, 1, NZ, BC_SYM),      # j+
            F(1, NX, 1, NY, 1, 1, BC_SYM),        # k-
            F(1, NX, 1, NY, NZ, NZ, BC_SYM)]      # k+

lines = ["! Gridgen boundary condition file", "%8d" % len(BLOCKS)]
for (_, _, _, NX, name) in BLOCKS:
    lines += ["%8d%8d%8d" % (NX, NY, NZ), name, "6"]
    lines += faces(NX)
open("bc3d.inp", "w").write("\n".join(lines) + "\n")

# ------------------------------------------------- bc3d_interface.inp
# No interfaces in this case: every block lists its 6 faces with positive codes
# and no pairing line follows (Convert_bc() only builds a pairing descriptor for
# bc<0 entries).  The file still has to exist because the reader opens it.
lines = ["! Interface connection file", "%8d" % len(BLOCKS)]
for (_, _, _, NX, name) in BLOCKS:
    lines += ["%8d%8d%8d" % (NX, NY, NZ), name, "6"]
    lines += faces(NX)
open("bc3d_interface.inp", "w").write("\n".join(lines) + "\n")

# ------------------------------------------------------------- material.in
# n_block, block types (1 = BLOCK_SOLID), then rho, Cp, k per block
mat = ["%d" % len(BLOCKS), " ".join(["1"] * len(BLOCKS))]
for (_, _, kk, _, _) in BLOCKS:
    mat.append("7800.0 460.0 %g" % kk)
open("material.in", "w").write("\n".join(mat) + "\n")

# ------------------------------------------------------------- solid_bc.inp
# per block: block no, n_face, then (face Tw Qw) lines
#   face 1 = i- : 400 K,   face 4 = i+ : 350 K  -> linear T(x) in every block
bc = ["%d" % len(BLOCKS)]
for m in range(1, len(BLOCKS) + 1):
    bc += ["%d" % m, "2",
           "%d %.1f 0.0" % (1, 400.0),
           "%d %.1f 0.0" % (4, 350.0)]
open("solid_bc.inp", "w").write("\n".join(bc) + "\n")

# ---------------------------------------------------- partation_np2.dat
# blk-1/2 -> rank 0, blk-3/4 -> rank 1 (two SOLID blocks per process)
det = ["  Total_block_number    Total_Proc_number",
       "%12d%12d" % (len(BLOCKS), 2),
       " Block_no, Proc_no"]
for m in range(1, len(BLOCKS) + 1):
    det.append("%12d%12d" % (m, 0 if m <= 2 else 1))
open("partation_np2.dat", "w").write("\n".join(det) + "\n")

print("input files written (bc3d.inp, bc3d_interface.inp, material.in, "
      "solid_bc.inp, partation_np2.dat)")
print("Mesh3d.x: %d blocks" % len(BLOCKS))
for (xr, yr, kk, NX, name) in BLOCKS:
    print("  %-30s x=[%g,%g] y=[%g,%g] k=%g" % (name, xr[0], xr[1], yr[0], yr[-1], kk))
