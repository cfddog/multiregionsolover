#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare u(y) profiles at x = Lx/2 between saved flow snapshots
(snap_N_block_*.vtk) to check convergence of the outer block iteration."""
import sys
import glob
import numpy as np
sys.path.insert(0, 'plot_validation')
import pv_vtk

files = sorted(glob.glob('snap_*_block_1.vtk'))
prof = {}
for fn in files:
    s = int(fn.split('_')[1])
    b1 = pv_vtk.load_block(fn)
    b2 = pv_vtk.load_block('snap_%d_block_2.vtk' % s)
    i0 = int(np.argmin(np.abs(b1['xc'][0, :] - 4.0)))
    prof[s] = np.concatenate([b2['u'][:, i0], b1['u'][:, i0]])
keys = sorted(prof)
for a, b in zip(keys[:-1], keys[1:]):
    d = np.max(np.abs(prof[a] - prof[b]))
    print('snap %d vs %d : max|du|=%.3e  rel=%.4e'
          % (a, b, d, d / np.max(np.abs(prof[b]))))
