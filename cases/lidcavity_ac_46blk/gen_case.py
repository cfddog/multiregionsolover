#!/usr/bin/env python3
"""Lid-driven cavity (42x42 cells, Re=100, U_lid=1) split into 4 and 6 blocks.

Same physical problem as cases/lidcavity_ac/re100 and cases/lidcavity_ac_3blk
(the 3-block case splits the same 42x42 cavity into 3 equal x-slabs), but with
two lattice splits:

    m4 : 2 x 2 blocks of 21x21 cells
    m6 : 3 x 2 blocks of 14x21 cells

and, for each split, TWO BLOCK NUMBERINGS:

    a : row-major, global block 1 = the bottom-left sub-domain
    b : reversed,  global block 1 = the top-right sub-domain

Both splits are fully closed (all outer walls no-slip, the top wall is the
moving lid -- BC_Wall on face 5, see sub_lowspeed_ac::ac_set_ghost_cell), so
the pressure level is a neutral mode.  In the multi-block AC solver the datum
is set by ONE block only:

    sub_lowspeed_ac.f90:217   if(.not. have_pbc .and. B%Block_no .eq. 1)
                                 call ac_zero_mean_p(nMesh, mBlock)

i.e. the GLOBAL block 1 removes its own mean pressure and every other block
follows through the p ghost halos.  The 'b' numbering therefore anchors the
pressure in the opposite corner of the cavity.

Purpose (task 8, anchoring sensitivity): the converged VELOCITY field must be
independent of which block is the anchor (it is a pure gauge choice), while p
may differ by a constant.  Any velocity difference beyond the solver
tolerance means the p datum is leaking into the momentum equations -- e.g.
through a per-block de-meaning of p.
"""
import struct
import os

NC = 42                     # cells per direction in the whole cavity
L = 1.0
NZ = 2                      # 1 cell in z (2D)

LAYOUTS = {'m2': (2, 1), 'm4': (2, 2), 'm6': (3, 2)}
TAGS = {'m2': ('a', 'b'), 'm4': ('a', 'b'), 'm6': ('a', 'b', 'c', 'd')}


def wrec(f, raw):
    f.write(struct.pack('I', len(raw)))
    f.write(raw)
    f.write(struct.pack('I', len(raw)))


def order_list(ncx, ncy, tag):
    """Physical (ix,iy) sub-domains in the file/block order of a layout.

    a: row-major, block 1 = bottom-left   (reference numbering)
    b: reversed,  block 1 = top-right
    c: block 1 = bottom-left, the rest reversed  (anchor fixed, order scrambled)
    d: column-major, block 1 = bottom-left       (order changed, anchor fixed)
    """
    row = [(ix, iy) for iy in range(1, ncy + 1) for ix in range(1, ncx + 1)]
    if tag == 'a':
        return row
    if tag == 'b':
        return row[::-1]
    if tag == 'c':
        return [row[0]] + row[1:][::-1]
    if tag == 'd':
        return [(ix, iy) for ix in range(1, ncx + 1) for iy in range(1, ncy + 1)]
    raise SystemExit('unknown numbering tag %r' % tag)


def write_mesh(fn, ncx, ncy, tag):
    nxb, nyb = NC // ncx, NC // ncy
    lst = order_list(ncx, ncy, tag)
    with open(fn, 'wb') as f:
        wrec(f, struct.pack('i', len(lst)))
        wrec(f, struct.pack('%di' % (3 * len(lst)),
                            *[v for _ in lst for v in (nxb + 1, nyb + 1, NZ)]))
        for ix, iy in lst:
            x0, y0 = (ix - 1) * L / ncx, (iy - 1) * L / ncy
            x, y, z = [], [], []
            for k in range(NZ):
                zk = k * 0.5 / (NZ - 1)
                for j in range(nyb + 1):
                    yj = y0 + (L / ncy) * j / nyb
                    for i in range(nxb + 1):
                        x.append(x0 + (L / ncx) * i / nxb)
                        y.append(yj)
                        z.append(zk)
            for a in (x, y, z):
                wrec(f, struct.pack('<%dd' % len(a), *a))


def write_bc(fn, ncx, ncy, tag):
    nxb, nyb = NC // ncx, NC // ncy
    lst = order_list(ncx, ncy, tag)
    num = {p: n + 1 for n, p in enumerate(lst)}       # physical -> block number
    Ls = []
    Ls.append('! %dx%d-block 42x42 lid cavity, numbering %s (block 1 = %s)'
              % (ncx, ncy, tag, 'bottom-left' if tag == 'a' else 'top-right'))
    Ls.append('%d' % len(lst))

    def face(ib, ie, jb, je, kb, ke, bc, nb=None):
        Ls.append('%4d %4d %4d %4d %4d %4d %4d' % (ib, ie, jb, je, kb, ke, bc))
        if nb is not None:
            Ls.append('%4d %4d %4d %4d %4d %4d %4d' % nb)

    for ix, iy in lst:
        Ls.append('%d %d %d' % (nxb + 1, nyb + 1, NZ))
        Ls.append('blk-%d' % num[(ix, iy)])
        Ls.append('6')
        # i- (face 1)
        if ix == 1:
            face(1, 1, 1, nyb + 1, 1, 2, 2)
        else:
            face(1, 1, 1, nyb + 1, -1, -2, -1,
                 (nxb + 1, nxb + 1, 1, nyb + 1, -1, -2, num[(ix - 1, iy)]))
        # i+ (face 4)
        if ix == ncx:
            face(nxb + 1, nxb + 1, 1, nyb + 1, 1, 2, 2)
        else:
            face(nxb + 1, nxb + 1, 1, nyb + 1, -1, -2, -1,
                 (1, 1, 1, nyb + 1, -1, -2, num[(ix + 1, iy)]))
        # j- (face 2)
        if iy == 1:
            face(1, nxb + 1, 1, 1, 1, 2, 2)
        else:
            face(1, nxb + 1, 1, 1, -1, -2, -1,
                 (1, nxb + 1, nyb + 1, nyb + 1, -1, -2, num[(ix, iy - 1)]))
        # j+ (face 5) -- the top row's j+ walls carry the moving lid (LS_U_lid)
        if iy == ncy:
            face(1, nxb + 1, nyb + 1, nyb + 1, 1, 2, 2)
        else:
            face(1, nxb + 1, nyb + 1, nyb + 1, -1, -2, -1,
                 (1, nxb + 1, 1, 1, -1, -2, num[(ix, iy + 1)]))
        # k- / k+ (faces 3 / 6): symmetry in the 2D case
        face(1, nxb + 1, 1, nyb + 1, 1, 1, 3)
        face(1, nxb + 1, 1, nyb + 1, 2, 2, 3)
    open(fn, 'w').write('\n'.join(Ls) + '\n')


def main():
    d = os.path.dirname(os.path.abspath(__file__))
    for name, (ncx, ncy) in LAYOUTS.items():
        nb = ncx * ncy
        for tag in TAGS[name]:
            write_mesh(os.path.join(d, 'Mesh3d_%s%s.x' % (name, tag)), ncx, ncy, tag)
            write_bc(os.path.join(d, 'bc3d_%s%s.inp' % (name, tag)), ncx, ncy, tag)
            print('wrote Mesh3d_%s%s.x + bc3d_%s%s.inp  (%d blocks of %dx%d cells)'
                  % (name, tag, name, tag, nb, NC // ncx, NC // ncy))
        # material.in carries Total_block (sub_read_parameter::read_material_in),
        # so every block count needs its own file.  All blocks are the same
        # class (2 = BLOCK_LOWSPEED): same-class BC_Inner needs that.
        open(os.path.join(d, 'material_%s.in' % name), 'w').write(
            '%d\n%s\n' % (nb, ' '.join(['2'] * nb)) +
            '\n'.join(['1.0 1.0 0.0'] * nb) + '\n')


if __name__ == '__main__':
    main()
