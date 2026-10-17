#!/usr/bin/env python3
"""Lid-driven cavity (Re=100, AC-FV) split into THREE equal blocks in x.

    y=1  +--- lid u=1 (moving wall, j+ face of all blocks) ---+
         |                                                  |
         |      blk-1        |      blk-2      |   blk-3    |
    y=0  +------- wall -------+-----------------+-----------+
        x=0                 1/3               2/3           1

Each block is 14x42 cells (15x43x2 nodes) -> 42x42 cells in total, i.e. the
same cell count in y as the verified single-block case
cases/lidcavity_ac/re100 (40x40) and just one finer in x, so the converged
u(x=0.5)/v(y=0.5) profiles must agree with that case and with Ghia (1982)
to within the same discretisation error.

Two conformal BC_Inner (-1) interfaces (blk1 i+ <-> blk2 i=1,
blk2 i+ <-> blk3 i=1).  No bc3d_interface.inp is written: the pairing is
fully given inline in bc3d.inp, exactly like cases/channel_ac_4blk.

The domain is CLOSED (no pressure outlet) -> the AC pressure level is
anchored by global block 1 (see lowspeed_ac_solver_one_block).
"""
import struct
import os

L, DZ = 1.0, 1.0
NI, NJ, NZ = 15, 43, 2          # 14x42 cells per block, 3 blocks -> 42x42


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
    xs = [(0.0, L / 3.0), (L / 3.0, 2.0 * L / 3.0), (2.0 * L / 3.0, L)]

    with open(os.path.join(d, 'Mesh3d.x'), 'wb') as f:
        wrec(f, struct.pack('i', 3))
        wrec(f, struct.pack('9i', *([NI, NJ, NZ] * 3)))
        for x0, x1 in xs:
            for a in xyz(x0, x1, 0.0, L):
                wrec(f, struct.pack('<%dd' % len(a), *a))

    L_ = []
    L_.append('! 3-block lid-driven cavity (Re=100, AC), equal x split')
    L_.append('3')

    def face(ib, ie, jb, je, kb, ke, bcv, nb=None):
        L_.append('%4d %4d %4d %4d %4d %4d %4d' % (ib, ie, jb, je, kb, ke, bcv))
        if nb is not None:
            L_.append('%4d %4d %4d %4d %4d %4d %4d' % nb)

    def block(name, own):
        L_.append('%d %d %d' % (NI, NJ, NZ))
        L_.append(name)
        L_.append('6')
        for ib, ie, jb, je, kb, ke, bcv, nb in own:
            face(ib, ie, jb, je, kb, ke, bcv, nb)

    wall, sym = 2, 3
    # ---- blk-1 (0 <= x <= 1/3): i- wall, i+ -> blk2, lid at j+ ----
    block('blk-1', [
        (1, 1, 1, NJ, 1, 2, wall, None),
        (NI, NI, 1, NJ, -1, -2, -1, (1, 1, 1, NJ, -1, -2, 2)),      # -> blk2 i=1
        (1, NI, 1, 1, 1, 2, wall, None),
        (1, NI, NJ, NJ, 1, 2, wall, None),                          # moving lid
        (1, NI, 1, NJ, 1, 1, sym, None),
        (1, NI, 1, NJ, 2, 2, sym, None)])
    # ---- blk-2 (1/3 <= x <= 2/3): i- -> blk1, i+ -> blk3 ----
    block('blk-2', [
        (1, 1, 1, NJ, -1, -2, -1, (NI, NI, 1, NJ, -1, -2, 1)),      # -> blk1 i=NI
        (NI, NI, 1, NJ, -1, -2, -1, (1, 1, 1, NJ, -1, -2, 3)),      # -> blk3 i=1
        (1, NI, 1, 1, 1, 2, wall, None),
        (1, NI, NJ, NJ, 1, 2, wall, None),
        (1, NI, 1, NJ, 1, 1, sym, None),
        (1, NI, 1, NJ, 2, 2, sym, None)])
    # ---- blk-3 (2/3 <= x <= 1): i- -> blk2, i+ wall ----
    block('blk-3', [
        (1, 1, 1, NJ, -1, -2, -1, (NI, NI, 1, NJ, -1, -2, 2)),      # -> blk2 i=NI
        (NI, NI, 1, NJ, 1, 2, wall, None),
        (1, NI, 1, 1, 1, 2, wall, None),
        (1, NI, NJ, NJ, 1, 2, wall, None),
        (1, NI, 1, NJ, 1, 1, sym, None),
        (1, NI, 1, NJ, 2, 2, sym, None)])

    txt = '\n'.join(L_) + '\n'
    open(os.path.join(d, 'bc3d.inp'), 'w').write(txt)
    open(os.path.join(d, 'material.in'), 'w').write(
        '3\n2 2 2\n' + '1.0 1.0 0.0\n' * 3)

    # ---- partation.dat for np=2 (1+2) and np=3 (1+1+1) --------------------
    # 注意（实测约束）：每个进程按 Mesh3d.x 的块序依次调用各块的整个 AC 伪时间
    # 循环，界面 halo 是每一步的阻塞收/发。若把"夹在中间"的 blk-2 单独放在一个
    # 进程、而 blk-1/blk-3 都在另一个进程，则 blk-1 第 2 步要等 blk-2 第 2 步的
    # 数据，blk-2 第 1 步又要等 blk-3 的数据，而 blk-3 所在的进程还没走到它
    # => 死锁。故 np=2 采用 {1,2} | {3}（等价于"1 个进程 2 块、另 1 个进程 1 块"）。
    for np_, mapping in ((2, {1: 0, 2: 0, 3: 1}), (3, {1: 0, 2: 1, 3: 2})):
        lines = ['  Total_block_number    Total_Proc_number',
                 '%12d%12d' % (3, np_),
                 '  Block_no, Proc_no']
        for b in (1, 2, 3):
            lines.append('%12d%12d' % (b, mapping[b]))
        open(os.path.join(d, 'partation_np%d.dat' % np_), 'w').write(
            '\n'.join(lines) + '\n')

    print('wrote Mesh3d.x (3 blocks %dx%dx%d), bc3d.inp, material.in, '
          'partation_np{2,3}.dat' % (NI, NJ, NZ))


if __name__ == '__main__':
    main()
