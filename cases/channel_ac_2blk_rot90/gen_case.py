#!/usr/bin/env python3
"""Two-block channel whose interface is a 90-degree rotation IN INDEX SPACE.

Physical mesh = cases/channel_ac_2blk (Lx=8, H=1, inlet u=1, outlet p=0,
no-slip walls, Re_H=50), i.e. 80+80 cells in x and 40 cells in y.  The ONLY
difference is that block 2 is written with its index axes swapped:

    blk-1:  i -> x (81 nodes, 80 cells),  j -> y (41 nodes, 40 cells)
    blk-2:  i -> y (41 nodes, 40 cells),  j -> x (81 nodes, 80 cells)

so the conformal interface (plane x = Lx/2) is

    blk-1  face 4 (i+)  <->  blk-2  face 2 (j-)

and Convert_bc must build the connection descriptor

    variant 'p' (blk-2 i along +y, its k along -z):  L = ( 2,  1, -3)
    variant 'n' (blk-2 i along -y, its k along +z):  L = ( 2, -1,  3)

instead of the usual (1,2,3).  This is exactly what an "i<->j" interface is:
an axis-aligned pair of hexahedral blocks can NOT form a physical 90-degree
bend with ONE full-face conformal interface (the two faces would have to be
perpendicular planes), so the 90 degree rotation can only live in index space.
The physics is unchanged, hence both variants MUST reproduce the non-rotated
two-block result of cases/channel_ac_2blk (u_max = 1.49546, interface column
sum du = 0.12 %, outlet column sum u = 40).

IMPORTANT (2026-10-18): the index triad of a block must be RIGHT-handed,
(x_i)x(x_j).(x_k) > 0, otherwise the cell volume is negative and the solver
blows up on the first pseudo step (it now aborts with
"left-handed block(s) in the input mesh", see sub_geometry.f90).
For variant 'p' that means blk-2 must run its k index along -z: with i -> +y
and j -> +x the product (y)x(x) is -z, so k must follow -z (which also turns
the k-connection into L3 = -3, encoded in bc3d.inp by writing the neighbour's
k range in the opposite order).  Variant 'n' (i -> -y) is right-handed with
k -> +z, hence L = (2,-1,3).

The two variants are written side by side (files *_p / *_n) and prep.sh makes
one run directory per (variant, np).
"""
import struct
import os

LX, H, DZ = 8.0, 1.0, 1.0
NX1, NY1, NZ = 81, 41, 2      # blk-1: 80 cells in x, 40 cells in y
NI2, NJ2 = 41, 81             # blk-2 (transposed): i -> y, j -> x


def wrec(f, raw):
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def blk1_xyz():
    """blk-1: i fastest along +x, j along +y."""
    x, y, z = [], [], []
    for k in range(NZ):
        zk = k * DZ / (NZ - 1)
        for j in range(NY1):
            yj = H * j / (NY1 - 1)
            for i in range(NX1):
                x.append(LX / 2 * i / (NX1 - 1))
                y.append(yj)
                z.append(zk)
    return x, y, z


def blk2_xyz(variant):
    """blk-2: i fastest along y (transposed), j along +x.

    variant 'p': i runs with +y   (i=1 -> y=0,   i=NI2 -> y=H)
    variant 'n': i runs with -y   (i=1 -> y=H,   i=NI2 -> y=0)

    Variant 'p' must also run k along -z so that the index triad stays
    RIGHT-handed:  (i)x(j) = (+y)x(+x) = -z  =>  k must follow -z
    (otherwise vol < 0 and the solver aborts on the left-handed mesh).
    """
    x, y, z = [], [], []
    for k in range(NZ):
        if variant == 'p':
            zk = DZ * (1.0 - k / (NZ - 1))
        else:
            zk = k * DZ / (NZ - 1)
        for j in range(NJ2):                     # j -> x
            xj = LX / 2 + LX / 2 * j / (NJ2 - 1)
            for i in range(NI2):                 # i -> y (transposed)
                if variant == 'p':
                    yi = H * i / (NI2 - 1)
                else:
                    yi = H * (1.0 - i / (NI2 - 1))
                x.append(xj)
                y.append(yi)
                z.append(zk)
    return x, y, z


def write_mesh(fn, variant):
    with open(fn, 'wb') as f:
        wrec(f, struct.pack('i', 2))
        wrec(f, struct.pack('6i', NX1, NY1, NZ, NI2, NJ2, NZ))
        for a in blk1_xyz():
            wrec(f, struct.pack('<%dd' % len(a), *a))
        for a in blk2_xyz(variant):
            wrec(f, struct.pack('<%dd' % len(a), *a))



def write_bc(fn, variant):
    L = []
    L.append('! 2-block channel, i<->j (90 deg index rotation) at x=4; variant '
             '%s (blk-2 i along %s, k along %s => right-handed)'
             % (variant, '+y' if variant == 'p' else '-y',
                '-z' if variant == 'p' else '+z'))
    L.append('2')

    def face(ib, ie, jb, je, kb, ke, bcv, nb=None):
        L.append('%4d %4d %4d %4d %4d %4d %4d' % (ib, ie, jb, je, kb, ke, bcv))
        if nb is not None:      # BC_Inner continuation: neighbour extent + blk no.
            L.append('%4d %4d %4d %4d %4d %4d %4d' % nb)

    # k-range marker of the PARTNER block on an interface entry.  The k range
    # in bc3d.inp is an *ordering* marker (abs+min/max, sign = direction), and
    # Convert_bc reverses the connection sign when
    # (ke-kb)*(ke1-kb1) < 0.  Variant 'p' has blk-2's k along -z, so its
    # k-connection to blk-1 (k -> +z) is reversed => L3 = -3.
    KP = (-2, -1) if variant == 'p' else (-1, -2)

    # ---- blk-1 : x in [0,4], i -> x, j -> y -------------------------------
    L.append('%d %d %d' % (NX1, NY1, NZ))
    L.append('blk-1')
    L.append('6')
    face(1, 1, 1, NY1, 1, 2, 5)                              # i- inlet (u=1)
    # i+ -> blk-2 j=1.  The target extent is blk-2's (i,j,k) sub-face, taken in
    # the direction that matches physically: variant p -> i 1..NI2 (y=0 at
    # i=1), variant n -> i NI2..1 (y=0 sits at i=NI2) => L2 = +1 / -1.
    if variant == 'p':
        face(NX1, NX1, 1, NY1, -1, -2, -1, (1, NI2, 1, 1, KP[0], KP[1], 2))
    else:
        face(NX1, NX1, 1, NY1, -1, -2, -1, (NI2, 1, 1, 1, KP[0], KP[1], 2))
    face(1, NX1, 1, 1, 1, 2, 2)                              # j- wall
    face(1, NX1, NY1, NY1, 1, 2, 2)                          # j+ wall
    face(1, NX1, 1, NY1, 1, 1, 3)                            # k- sym
    face(1, NX1, 1, NY1, 2, 2, 3)                            # k+ sym

    # ---- blk-2 : x in [4,8] but TRANSPOSED (i -> y, j -> x) ---------------
    L.append('%d %d %d' % (NI2, NJ2, NZ))
    L.append('blk-2')
    L.append('6')
    face(1, 1, 1, NJ2, 1, 2, 2)                              # i- wall
    face(NI2, NI2, 1, NJ2, 1, 2, 2)                          # i+ wall
    # j- -> blk-1 i=NX1.  For variant n the block's own i axis is reversed
    # (i increases while y decreases), so the sub-face is written reversed:
    # that is what makes Convert_bc produce the negative ordering sign L1=-2.
    if variant == 'p':
        face(1, NI2, 1, 1, -1, -2, -1, (NX1, NX1, 1, NY1, KP[0], KP[1], 1))
    else:
        face(NI2, 1, 1, 1, -1, -2, -1, (NX1, NX1, 1, NY1, KP[0], KP[1], 1))
    face(1, NI2, NJ2, NJ2, 1, 2, 6)                          # j+ outlet (p=0)
    face(1, NI2, 1, NJ2, 1, 1, 3)                            # k- sym
    face(1, NI2, 1, NJ2, 2, 2, 3)                            # k+ sym

    open(fn, 'w').write('\n'.join(L) + '\n')


def main():
    d = os.path.dirname(os.path.abspath(__file__))
    for v in ('p', 'n'):
        write_mesh(os.path.join(d, 'Mesh3d_%s.x' % v), v)
        write_bc(os.path.join(d, 'bc3d_%s.inp' % v), v)
    # same-class AC blocks: no bc3d_interface.inp is needed (as in
    # channel_ac_2blk / lidcavity_ac_3blk)
    open(os.path.join(d, 'material.in'), 'w').write(
        '2\n2 2\n1.0 1.0 0.0\n1.0 1.0 0.0\n')
    print('wrote Mesh3d_{p,n}.x (%d+%d / %d+%d nodes), bc3d_{p,n}.inp, '
          'material.in' % (NX1, NY1, NI2, NJ2))


if __name__ == '__main__':
    main()
