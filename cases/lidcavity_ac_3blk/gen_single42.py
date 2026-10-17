#!/usr/bin/env python3
"""Single-block 42x42 (43x43x2 nodes) reference mesh on the *identical* grid
used by the 3-block case (blocks of 14x42 cells, x split at 1/3 and 2/3).
Purpose: give a gold-standard single-block solution to compare the 3-block
(np=1/2/3) results against, with no discretisation or grid difference."""
import struct
import os

L, DZ = 1.0, 1.0
NI, NJ, NZ = 43, 43, 2


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
    with open(os.path.join(d, 'Mesh3d.x'), 'wb') as f:
        wrec(f, struct.pack('i', 1))
        wrec(f, struct.pack('3i', NI, NJ, NZ))
        for a in xyz(0.0, L, 0.0, L):
            wrec(f, struct.pack('<%dd' % len(a), *a))

    wall, sym = 2, 3
    L_ = []
    L_.append('! single-block 42x42 lid-driven cavity (same grid as 3-block case)')
    L_.append('1')
    L_.append('%5d%6d%6d' % (NI, NJ, NZ))
    L_.append('blk-1')
    L_.append('6')

    def face(ib, ie, jb, je, kb, ke, bcv):
        L_.append('%8d%8d%8d%8d%8d%8d%8d' % (ib, ie, jb, je, kb, ke, bcv))

    face(1, 1, 1, NJ, 1, 2, wall)
    face(NI, NI, 1, NJ, 1, 2, wall)
    face(1, NI, 1, 1, 1, 2, wall)
    face(1, NI, NJ, NJ, 1, 2, wall)     # moving lid
    face(1, NI, 1, NJ, 1, 1, sym)
    face(1, NI, 1, NJ, 2, 2, sym)
    open(os.path.join(d, 'bc3d.inp'), 'w').write('\n'.join(L_) + '\n')
    open(os.path.join(d, 'material.in'), 'w').write('1\n2\n' + '1.0 1.0 0.0\n')
    print('wrote single-block %dx%dx%d' % (NI, NJ, NZ))


if __name__ == '__main__':
    main()
