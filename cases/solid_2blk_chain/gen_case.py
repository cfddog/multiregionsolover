#!/usr/bin/env python3
"""Two-block SOLID<->SOLID series-conduction case (same-class internal link).

Topology (2-D in x-y, 3 nodes in z, conformal 1:1 interface):

  blk-1  x in [0.00, 0.05]  k =  16.3  (steel)   i- wall 400 K, i+ = link
  blk-2  x in [0.05, 0.10]  k = 200.0  (copper)  i- = link,    i+ wall 350 K

  j-/j+/k-/k+ : symmetry (adiabatic)

File convention (do not mix the two roles)
-----------------------------------------
  bc3d_interface.inp : declares the block CONNECTION.  A double-connected
                       (对接) face is coded -1 (BC_Inner) and is followed by
                       the pairing line  ib1 ie1 jb1 je1 kb1 ke1 nb1 .
  bc3d.inp           : declares the PHYSICAL boundary / solver type.  The same
                       face must also be -1 + pairing line for a same-class
                       pair (solid-solid): convert_inp_inc only reads a
                       continuation line after a NEGATIVE code, and
                       check_bc3d requires a negative code for same-class
                       faces.  Cross-class pairs would use 11..19 here (the
                       canonical code is derived from material.in block types,
                       see prepare_bc3d_input).

  The runtime face type used by the solvers is bc3d.inp's code (align_interface
  copies it over bc3d_interface's), while the pairing (nb1/face1/L1..L3) comes
  from bc3d_interface.inp.  Hence both files must carry the -1 link.

Exact steady state (1-D, uniform dx = 0.05/40 = 1.25e-3 m):

  R1 = L1/k1, R2 = L2/k2,  q = (400-350)/(R1+R2)
  T(x) = 400 - q x / k1          (0 <= x <= 0.05)
       = 350 + q (0.1-x)/k2      (0.05 <= x <= 0.1)
  T_iface = 400 - q R1

With these numbers q = 15071 W/m2 and T_iface = 353.77 K.  The two CELL
CENTRES next to the interface differ by the physical FVM amount

  dT_cc = q*(dx/(2 k1) + dx/(2 k2)) = 0.625 K

(not zero -- that is the discretisation of dT/dx inside each half cell), so
the pass criterion of this case is

  |T - T_exact| < 1e-6 K   everywhere
  dT_cc within 1e-3 K of the analytic 0.625 K

Why the case exists
-------------------
The local solid-solid coupling (couple_solid_solid_interfaces / the same-
process branch of Ts_send_mpi) used to write the neighbour's values into
INTERIOR cells instead of the LAP ghost layers, which produced a large
spurious jump across a solid-solid interface.  This case pins the analytic
profile and the interface jump so the regression cannot come back unnoticed.

Writes Mesh3d.x, bc3d.inp, bc3d_interface.inp, material.in, solid_bc.inp.
Run with 'python3 gen_case.py' inside this directory.
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
ys = lin(0.0, 0.1, NY)
zs = lin(0.0, 0.1, NZ)

# (x-range, conductivity, nx, name)
BLOCKS = [((0.00, 0.05), 16.3,  41, 'blk-1 (41 x 21 x 3, k=16.3)'),
          ((0.05, 0.10), 200.0, 41, 'blk-2 (41 x 21 x 3, k=200)')]

def coords(xs_, ys_, zs_):
    xc, yc, zc = array.array('d'), array.array('d'), array.array('d')
    for k in range(len(zs_)):
        for j in range(len(ys_)):
            for i in range(len(xs_)):
                xc.append(xs_[i]); yc.append(ys_[j]); zc.append(zs_[k])
    return xc, yc, zc

# ------------------------------------------------------------------ Mesh3d.x
with open("Mesh3d.x", "wb") as f:
    record(f, [len(BLOCKS)], 'i')
    record(f, [d for (_, _, NX, _) in BLOCKS for d in (NX, NY, NZ)],
           'i' * (3 * len(BLOCKS)))
    for (xr, _, NX, _) in BLOCKS:
        for a in coords(lin(xr[0], xr[1], NX), ys, zs):
            record(f, a, '')

# --------------------------------------------------- bc3d.inp / bc3d_interface.inp
def F(i1, i2, j1, j2, k1, k2, bc):
    return "%8d%8d%8d%8d%8d%8d%8d" % (i1, i2, j1, j2, k1, k2, bc)

BC_WALL, BC_SYM, BC_LINK = 2, 3, -1
BC_SOLID_SOLID = 14          # BC_Interface_SolidSolid
NXS = [NX for (_, _, NX, _) in BLOCKS]

def faces(NX, iface_face, nb, with_link):
    """6 subfaces of one block.

    bc3d_interface.inp is the only file carrying the block connection: the
    double-connected (对接) face is coded -1 (BC_Inner) there and is followed by
    its pairing line `ib1 ie1 jb1 je1 kb1 ke1 nb1`, in which the coupling
    direction is marked negative (see docs/程序说明).

    bc3d.inp only classifies the face -- it uses the explicit solid-solid type
    BC_Interface_SolidSolid (14) and has NO continuation line, because the
    neighbour block / ordering already live in bc3d_interface.inp.
    The opposite i face is the isothermal wall, the remaining four symmetry.
    """
    out = []
    for f_id in (1, 4):
        i1 = 1 if f_id == 1 else NX
        if f_id == iface_face:
            if with_link:
                out.append(F(i1, i1, 1, NY, -1, -NZ, BC_LINK))
                i1n = 1 if f_id == 4 else NXS[nb - 1]
                out.append(F(i1n, i1n, 1, NY, -1, -NZ, nb))
            else:
                # type-only line: no pairing, so no directional sign either
                out.append(F(i1, i1, 1, NY, 1, NZ, BC_SOLID_SOLID))
        else:
            out.append(F(i1, i1, 1, NY, 1, NZ, BC_WALL))
    out += [F(1, NX, 1, 1, 1, NZ, BC_SYM),      # j-
            F(1, NX, NY, NY, 1, NZ, BC_SYM),    # j+
            F(1, NX, 1, NY, 1, 1, BC_SYM),      # k-
            F(1, NX, 1, NY, NZ, NZ, BC_SYM)]    # k+
    return out

# block 1 couples on its i+ face, block 2 on its i- face
IFACE_FACE = {1: 4, 2: 1}

def write_bc_file(fname, header, with_link):
    lines = [header, "%8d" % len(BLOCKS)]
    for m, (_, _, NX, name) in enumerate(BLOCKS, start=1):
        nb = 2 if m == 1 else 1
        lines += ["%8d%8d%8d" % (NX, NY, NZ), name, "6"]
        lines += faces(NX, IFACE_FACE[m], nb, with_link)
    open(fname, "w").write("\n".join(lines) + "\n")

write_bc_file("bc3d.inp",           "! Gridgen boundary condition file", False)
write_bc_file("bc3d_interface.inp", "! Interface connection file",       True)

# ------------------------------------------------------------- material.in
mat = ["%d" % len(BLOCKS), " ".join(["1"] * len(BLOCKS))]
for (_, kk, _, _) in BLOCKS:
    mat.append("7800.0 460.0 %g" % kk)
open("material.in", "w").write("\n".join(mat) + "\n")

# ------------------------------------------------------------- solid_bc.inp
# blk-1: i- = 400 K ;  blk-2: i+ = 350 K   (only physical walls)
bc = ["%d" % len(BLOCKS),
      "1", "1", "1 %.1f 0.0" % 400.0,
      "2", "1", "4 %.1f 0.0" % 350.0]
open("solid_bc.inp", "w").write("\n".join(bc) + "\n")

print("input files written (bc3d.inp, bc3d_interface.inp, material.in, "
      "solid_bc.inp); Mesh3d.x with %d blocks" % len(BLOCKS))
