#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Print per-row u(y) and error at x=Lx/2 for two runs."""
import sys
import numpy as np
sys.path.insert(0, 'plot_validation')
import pv_vtk
import bj_analytical as bj

EPS, DP, MU = 0.8, 0.5, 0.02
HP = 1.2
which = sys.argv[1] if len(sys.argv) > 1 else 'run_fine'
P = 6.249791e-4 if 'fine' in which else 6.248674e-4
if 'coarse' in which:
    P = 6.244429e-4
b1 = pv_vtk.load_block(which + '/flow3d_block_1.vtk')
b2 = pv_vtk.load_block(which + '/flow3d_block_2.vtk')
i0 = int(np.argmin(np.abs(b1['xc'][0, :] - 4.0)))
sol = bj.solve_composite(P, EPS, DP, MU, 1.0, 1.2)
umax = np.max(bj.profile(np.linspace(-1.2, 1.0, 2001), sol))
print(which, 'nyf', b1['u'].shape[0], 'nyp', b2['u'].shape[0])
for j in range(b2['u'].shape[0]):
    y = b2['yc'][j, i0]
    if j % 6 == 0 or abs(y) < 0.25:
        print(' y=%8.4f u=%10.3e an=%10.3e err=%9.2e'
              % (y, b2['u'][j, i0], bj.u_porous(y, sol),
                 (b2['u'][j, i0] - bj.u_porous(y, sol)) / umax))
for j in range(b1['u'].shape[0]):
    y = b1['yc'][j, i0]
    if j % 6 == 0 or y < 0.25:
        print(' y=%8.4f u=%10.3e an=%10.3e err=%9.2e'
              % (y, b1['u'][j, i0], bj.u_fluid(y, sol),
                 (b1['u'][j, i0] - bj.u_fluid(y, sol)) / umax))
