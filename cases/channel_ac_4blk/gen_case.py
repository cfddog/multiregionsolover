#!/usr/bin/env python3
"""4-block channel: 2x2 block decomposition of the 2D Poiseuille channel.

Same physics/resolution as cases/channel_ac (Lx=8, H=1, rho=1, mu=0.02,
u_in=1, p_out=0, no-slip walls) and cases/channel_ac_2blk (80x40 cells total),
but the domain is split in BOTH x and y:

      y
   0.5 +----- blk3 -----+----- blk4 -----+   (outlet at i+ of blk4)
       |                |                |
     0 +----- blk1 -----+----- blk2 -----+   (inlet at i- of blk1/blk3)
       0                4                8   x

so FOUR conformal -1 interfaces meet at the interior point (4, 0.5) and along
the two interior lines : this is the corner/edge test for the multi-block AC
coupling (analogue of the T-junction where three interfaces meet).
"""
import struct
import os

Lx, H, DZ = 8.0, 1.0, 1.0
NI, NJ, NZ = 41, 21, 2      # per block -> 40x20 cells, 2x2 blocks = 160x40 cells


def wrec(f, raw):
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def xyz(x0, x1, y0, y1):
    x, y, z = [], [], []
    for k in range(NZ):
        zk = k * DZ / (NZ - 1)
        for j in range(NJ):
            yj = y0 + (y1 - y0) * j / (NJ - 1)
            for i in range(NI):
                x.append(x0 + (x1 - x0) * i / (NI - 1))
                y.append(yj)
                z.append(zk)
    return x, y, z


def main():
    d = os.path.dirname(os.path.abspath(__file__))
    quads = [(0.0, Lx / 2, 0.0, H / 2),      # blk1
             (Lx / 2, Lx, 0.0, H / 2),       # blk2
             (0.0, Lx / 2, H / 2, H),        # blk3
             (Lx / 2, Lx, H / 2, H)]         # blk4
    with open(os.path.join(d, 'Mesh3d.x'), 'wb') as f:
        wrec(f, struct.pack('i', 4))
        wrec(f, struct.pack('12i', *([NI, NJ, NZ] * 4)))
        for x0, x1, y0, y1 in quads:
            for a in xyz(x0, x1, y0, y1):
                wrec(f, struct.pack('<%dd' % len(a), *a))

    L = []
    L.append('! 4-block channel 2x2, conformal -1 interfaces (corner test)')
    L.append('4')

    def face(ib, ie, jb, je, kb, ke, bcv, nb=None):
        L.append('%4d %4d %4d %4d %4d %4d %4d' % (ib, ie, jb, je, kb, ke, bcv))
        if nb is not None:
            L.append('%4d %4d %4d %4d %4d %4d %4d' % nb)

    def block(name, own):
        L.append('%d %d %d' % (NI, NJ, NZ))
        L.append(name)
        L.append('6')
        for ib, ie, jb, je, kb, ke, bcv, nb in own:
            face(ib, ie, jb, je, kb, ke, bcv, nb)

    # ---- blk1 (x<4, y<0.5): inlet, i+ ->blk2, j- wall, j+ ->blk3 ----
    block('blk-1', [
        (1, 1, 1, NJ, 1, 2, 5, None),                   # i- inlet
        (NI, NI, 1, NJ, -1, -2, -1, (1, 1, 1, NJ, -1, -2, 2)),    # i+ -> blk2 i=1
        (1, NI, 1, 1, 1, 2, 2, None),                   # j- wall
        (1, NI, NJ, NJ, -1, -2, -1, (1, NI, 1, 1, -1, -2, 3)),    # j+ -> blk3 j=1
        (1, NI, 1, NJ, 1, 1, 3, None),                  # k- sym
        (1, NI, 1, NJ, 2, 2, 3, None)])                 # k+ sym
    # ---- blk2 (x>4, y<0.5): i- ->blk1, outlet, j- wall, j+ ->blk4 ----
    block('blk-2', [
        (1, 1, 1, NJ, -1, -2, -1, (NI, NI, 1, NJ, -1, -2, 1)),    # i- -> blk1 i=NI
        (NI, NI, 1, NJ, 1, 2, 6, None),                 # i+ outlet
        (1, NI, 1, 1, 1, 2, 2, None),                   # j- wall
        (1, NI, NJ, NJ, -1, -2, -1, (1, NI, 1, 1, -1, -2, 4)),    # j+ -> blk4 j=1
        (1, NI, 1, NJ, 1, 1, 3, None),
        (1, NI, 1, NJ, 2, 2, 3, None)])
    # ---- blk3 (x<4, y>0.5): inlet, i+ ->blk4, j- ->blk1, j+ wall ----
    block('blk-3', [
        (1, 1, 1, NJ, 1, 2, 5, None),                   # i- inlet
        (NI, NI, 1, NJ, -1, -2, -1, (1, 1, 1, NJ, -1, -2, 4)),    # i+ -> blk4 i=1
        (1, NI, 1, 1, -1, -2, -1, (1, NI, NJ, NJ, -1, -2, 1)),    # j- -> blk1 j=NJ
        (1, NI, NJ, NJ, 1, 2, 2, None),                 # j+ wall
        (1, NI, 1, NJ, 1, 1, 3, None),
        (1, NI, 1, NJ, 2, 2, 3, None)])
    # ---- blk4 (x>4, y>0.5): i- ->blk3, outlet, j- ->blk2, j+ wall ----
    block('blk-4', [
        (1, 1, 1, NJ, -1, -2, -1, (NI, NI, 1, NJ, -1, -2, 3)),    # i- -> blk3 i=NI
        (NI, NI, 1, NJ, 1, 2, 6, None),                 # i+ outlet
        (1, NI, 1, 1, -1, -2, -1, (1, NI, NJ, NJ, -1, -2, 2)),    # j- -> blk2 j=NJ
        (1, NI, NJ, NJ, 1, 2, 2, None),                 # j+ wall
        (1, NI, 1, NJ, 1, 1, 3, None),
        (1, NI, 1, NJ, 2, 2, 3, None)])

    txt = '\n'.join(L) + '\n'
    open(os.path.join(d, 'bc3d.inp'), 'w').write(txt)
    open(os.path.join(d, 'bc3d_interface.inp'), 'w').write(txt)
    open(os.path.join(d, 'material.in'), 'w').write(
        '4\n2 2 2 2\n' + '1.0 1.0 0.0\n' * 4)
    print('wrote Mesh3d.x (4 blocks %dx%dx%d), bc3d.inp, bc3d_interface.inp, material.in'
          % (NI, NJ, NZ))


if __name__ == '__main__':
    main()
