#!/usr/bin/env python3
"""2-block channel: 2D Poiseuille channel split at x=Lx/2 into two conformal
blocks joined by ONE non-rotated BC_Inner (-1) interface.

Geometry / physics identical to cases/channel_ac (Lx=8, H=1, rho=1, mu=0.02,
uniform velocity inlet u=1, pressure outlet p=0, no-slip walls) so the
converged solution must reproduce the single-block result (u_max 1.495 vs
analytic 1.5) -- a clean test of the multi-block AC coupling.
"""
import struct
import os

Lx, H, DZ = 8.0, 1.0, 1.0
NX1, NX2, NY, NZ = 81, 81, 41, 2          # 40+40 cells = the 160 of channel_ac


def wrec(f, raw):
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def blk(x0, x1, dx):
    x, y, z = [], [], []
    for k in range(NZ):
        zk = k * DZ / (NZ - 1)
        for j in range(NY):
            yj = H * j / (NY - 1)
            for i in range(NX1):
                x.append(x0 + (x1 - x0) * i / (NX1 - 1))
                y.append(yj)
                z.append(zk)
    return x, y, z


def main():
    d = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(d, 'Mesh3d.x'), 'wb') as f:
        wrec(f, struct.pack('i', 2))
        wrec(f, struct.pack('6i', NX1, NY, NZ, NX2, NY, NZ))
        for x0, x1 in ((0.0, Lx / 2), (Lx / 2, Lx)):
            for a in blk(x0, x1, 0):
                wrec(f, struct.pack('<%dd' % len(a), *a))

    def face(L, ib, ie, jb, je, kb, ke, bcv, nb=None):
        L.append('%4d %4d %4d %4d %4d %4d %4d' % (ib, ie, jb, je, kb, ke, bcv))
        if nb is not None:          # BC_Inner continuation: neighbour extent + block no.
            L.append('%4d %4d %4d %4d %4d %4d %4d' % nb)

    bc = []
    bc.append('! 2-block channel, conformal -1 interface at x=Lx/2')
    bc.append('2')
    # ---- block 1 : x in [0,Lx/2] ----
    bc.append('%d %d %d' % (NX1, NY, NZ))
    bc.append('blk-1')
    bc.append('6')
    face(bc, 1, 1, 1, NY, 1, 2, 5)                                   # i- inlet
    face(bc, NX1, NX1, 1, NY, -1, -2, -1, (1, 1, 1, NY, -1, -2, 2))  # i+ -> blk2 i=1
    face(bc, 1, NX1, 1, 1, 1, 2, 2)                                  # j- wall
    face(bc, 1, NX1, NY, NY, 1, 2, 2)                                # j+ wall
    face(bc, 1, NX1, 1, NY, 1, 1, 3)                                 # k- sym
    face(bc, 1, NX1, 1, NY, 2, 2, 3)                                 # k+ sym
    # ---- block 2 : x in [Lx/2,Lx] ----
    bc.append('%d %d %d' % (NX2, NY, NZ))
    bc.append('blk-2')
    bc.append('6')
    face(bc, 1, 1, 1, NY, -1, -2, -1, (NX1, NX1, 1, NY, -1, -2, 1))  # i- -> blk1 i=41
    face(bc, NX2, NX2, 1, NY, 1, 2, 6)                               # i+ outlet
    face(bc, 1, NX2, 1, 1, 1, 2, 2)                                  # j- wall
    face(bc, 1, NX2, NY, NY, 1, 2, 2)                                # j+ wall
    face(bc, 1, NX2, 1, NY, 1, 1, 3)                                 # k- sym
    face(bc, 1, NX2, 1, NY, 2, 2, 3)                                 # k+ sym
    txt = '\n'.join(bc) + '\n'
    open(os.path.join(d, 'bc3d.inp'), 'w').write(txt)
    open(os.path.join(d, 'bc3d_interface.inp'), 'w').write(txt)
    open(os.path.join(d, 'material.in'), 'w').write('2\n2 2\n1.0 1.0 0.0\n1.0 1.0 0.0\n')
    print('wrote Mesh3d.x (%d+%d nodes), bc3d.inp, bc3d_interface.inp, material.in'
          % (NX1, NX2))


if __name__ == '__main__':
    main()
