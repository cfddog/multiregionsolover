#!/usr/bin/env python3
"""2D right-angle T-junction (spec B1): FOUR conformal LOWSPEED blocks, THREE
`bc<0` inner connections -- one of them **rotated (i<->j)**.

Geometry (W = 1, all dimensionless):
    main duct : x in [-10, 10], y in [0, 1]     (upstream L1 = 10W, downstream L2 = 10W)
    branch    : x in [-0.5, 0.5], y in [1, 11]  (length 10W, width W = 1:1)

Block layout (nodes) and BC codes (2=wall, 3=symmetry, 5=inlet, 6=outlet,
-1=inner same-class connection):

        D (branch, ROTATED mesh: i along +y, j along +x, 101 x 11)
        +--------+  i+ = outlet (6);  j- / j+ = branch side walls (2)
        |        |
  A     +-- B --+     C
  x[-10,-0.5]  x[-0.5,0.5]  x[0.5,10]        C i+ = downstream outlet (6)
  96x11         11x11        96x11

  A: i- inlet(5) | i+ -> B | j- wall | j+ wall | k+- symmetry
  B: i- -> A | i+ -> C | j- wall | j+ -> D | k+- symmetry
  C: i- -> B | i+ outlet(6) | j- wall | j+ wall | k+- symmetry
  D (rotated): i- -> B | i+ outlet(6) | j- wall | j+ wall | k+- symmetry

Cell size is uniform (dx = dy = 0.1) in every block, so all three interfaces are
conformal (11 nodes across): A(i+)<->B(i-) and B(i+)<->C(i-) are non-rotated while
B(j+)<->D(i-) is a rotated (i<->j) connection -- the P3 target.  The three
interfaces also meet at the corner points (-0.5, 1) and (0.5, 1).

Re = rho*U*W/mu is set in control.ec through LS_mu (rho=1, U=1, W=1):
Re=100 -> mu=0.01, Re=500 -> mu=2e-3, Re=1000 -> mu=1e-3.

Run:  python3 gen_case.py   then   ../../src/opencfd-ec1.16a.out
"""
import struct
import os

# ---- grid (all sizes are NODE counts) ---------------------------------------
NI_A, NJ_A = 96, 11        # upstream   95 x 10 cells
NI_B, NJ_B = 11, 11        # junction   10 x 10 cells
NI_C, NJ_C = 96, 11        # downstream 95 x 10 cells
NI_D, NJ_D = 101, 11       # branch ROTATED: i along y (100 cells), j along x (10)
NZ = 2

XA0, XA1 = -10.0, -0.5
XB0, XB1 = -0.5, 0.5
XC0, XC1 = 0.5, 10.0
Y0, Y1 = 0.0, 1.0
YD0, YD1 = 1.0, 11.0

WALL, SYM, INLET, OUTLET, INNER = 2, 3, 5, 6, -1


def wrec(f, raw):
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def coords(ni, nj, nz, xfun, yfun):
    """i fastest, then j, then k (same layout as the other generated cases)."""
    x, y, z = [], [], []
    for k in range(nz):
        for j in range(nj):
            for i in range(ni):
                x.append(xfun(i, j))
                y.append(yfun(i, j))
                z.append(k * 1.0 / (nz - 1))
    return x, y, z


def lin(a, b, n, t):
    return a + (b - a) * t / (n - 1)


def main():
    d = os.path.dirname(os.path.abspath(__file__))

    # (node dims, coordinate map) -- D is ROTATED: i = y, j = x
    blocks = [(NI_A, NJ_A, lambda i, j: (lin(XA0, XA1, NI_A, i), lin(Y0, Y1, NJ_A, j))),
              (NI_B, NJ_B, lambda i, j: (lin(XB0, XB1, NI_B, i), lin(Y0, Y1, NJ_B, j))),
              (NI_C, NJ_C, lambda i, j: (lin(XC0, XC1, NI_C, i), lin(Y0, Y1, NJ_C, j))),
              (NI_D, NJ_D, lambda i, j: (lin(XB0, XB1, NJ_D, j), lin(YD0, YD1, NI_D, i)))]

    with open(os.path.join(d, 'Mesh3d.x'), 'wb') as f:
        wrec(f, struct.pack('i', len(blocks)))
        wrec(f, struct.pack('%di' % (3 * len(blocks)),
                            *[v for (ni, nj, _) in blocks for v in (ni, nj, NZ)]))
        for (ni, nj, fn) in blocks:
            for a in coords(ni, nj, NZ,
                            lambda i, j, fn=fn: fn(i, j)[0],
                            lambda i, j, fn=fn: fn(i, j)[1]):
                wrec(f, struct.pack('<%dd' % len(a), *a))

    L = ['! 2D right-angle T-junction (B1): 4 conformal blocks, 3 inner (-1) links,'
         ' B(j+)<->D(i-) rotated', '%d' % len(blocks)]

    def entry(ib, ie, jb, je, bcv, nb=None):
        """Physical face: (ib,ie,jb,je) + bc in the last column.
        Inner link: connection marker in the kb/ke columns and the neighbour's
        face ranges + block number on the following line (same style as the other
        generated cases, e.g. channel_ac_4blk)."""
        if nb is None:
            L.append('%4d %4d %4d %4d %4d %4d %7d' % (ib, ie, jb, je, 1, 2, bcv))
        else:
            L.append('%4d %4d %4d %4d %4d %4d %7d' % (ib, ie, jb, je, -1, -2, bcv))
            L.append('%4d %4d %4d %4d %4d %4d %7d' % nb)

    def block(ni, nj, name, faces):
        L.append('%d %d %d' % (ni, nj, NZ))
        L.append(name)
        L.append('%d' % len(faces))
        for e in faces:
            entry(*e)

    # ---- A upstream (i = x, j = y) ------------------------------------------
    block(NI_A, NJ_A, 'Upstream', [
        (1, 1, 1, NJ_A, INLET),                                          # i- inlet
        (NI_A, NI_A, 1, NJ_A, INNER, (1, 1, 1, NJ_B, -1, -2, 2)),        # i+ -> B i-
        (1, NI_A, 1, 1, WALL),                                           # j- wall
        (1, NI_A, NJ_A, NJ_A, WALL),                                     # j+ wall
        (1, NI_A, 1, NJ_A, SYM),                                         # k- symmetry
        (1, NI_A, 1, NJ_A, SYM)])                                        # k+ symmetry
    # ---- B junction ---------------------------------------------------------
    block(NI_B, NJ_B, 'Junction', [
        (1, 1, 1, NJ_B, INNER, (NI_A, NI_A, 1, NJ_A, -1, -2, 1)),        # i- -> A i+
        (NI_B, NI_B, 1, NJ_B, INNER, (1, 1, 1, NJ_C, -1, -2, 3)),        # i+ -> C i-
        (1, NI_B, 1, 1, WALL),                                           # j- wall
        (1, NI_B, NJ_B, NJ_B, INNER, (1, 1, 1, NJ_D, -1, -2, 4)),        # j+ -> D i- ROTATED
        (1, NI_B, 1, NJ_B, SYM),
        (1, NI_B, 1, NJ_B, SYM)])
    # ---- C downstream -------------------------------------------------------
    block(NI_C, NJ_C, 'Downstream', [
        (1, 1, 1, NJ_C, INNER, (NI_B, NI_B, 1, NJ_B, -1, -2, 2)),        # i- -> B i+
        (NI_C, NI_C, 1, NJ_C, OUTLET),                                   # i+ outlet
        (1, NI_C, 1, 1, WALL),
        (1, NI_C, NJ_C, NJ_C, WALL),
        (1, NI_C, 1, NJ_C, SYM),
        (1, NI_C, 1, NJ_C, SYM)])
    # ---- D branch (ROTATED: i = y, j = x) -----------------------------------
    block(NI_D, NJ_D, 'Branch', [
        (1, 1, 1, NJ_D, INNER, (1, NI_B, NJ_B, NJ_B, -1, -2, 2)),        # i- -> B j+ ROTATED
        (NI_D, NI_D, 1, NJ_D, OUTLET),                                   # i+ branch outlet
        (1, NI_D, 1, 1, WALL),                                           # j- wall (x=-0.5)
        (1, NI_D, NJ_D, NJ_D, WALL),                                     # j+ wall (x=+0.5)
        (1, NI_D, 1, NJ_D, SYM),
        (1, NI_D, 1, NJ_D, SYM)])

    open(os.path.join(d, 'bc3d.inp'), 'w').write('\n'.join(L) + '\n')
    open(os.path.join(d, 'material.in'), 'w').write(
        '%d\n' % len(blocks) + ' '.join(['2'] * len(blocks)) + '\n' +
        '1.0 1.0 0.0\n' * len(blocks))
    print('wrote Mesh3d.x (%d blocks), bc3d.inp, material.in' % len(blocks))
    print('  nodes: A %dx%d  B %dx%d  C %dx%d  D %dx%d (rotated)' %
          (NI_A, NJ_A, NI_B, NJ_B, NI_C, NJ_C, NI_D, NJ_D))


if __name__ == '__main__':
    main()
