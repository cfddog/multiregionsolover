#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Print coarse/fine/medium velocity profile tables at x=Lx/2."""
import sys
import numpy as np
sys.path.insert(0, 'plot_validation')
import pv_vtk
import bj_analytical as bj

EPS, DP, MU = 0.8, 0.5, 0.02
HP = 1.2
for tag, P in [('run_coarse', 6.244429e-4), ('run_med', 6.248674e-4),
               ('run_fine', 6.249791e-4)]:
    b1 = pv_vtk.load_block(tag + '/flow3d_block_1.vtk')
    b2 = pv_vtk.load_block(tag + '/flow3d_block_2.vtk')
    i0 = int(np.argmin(np.abs(b1['xc'][0, :] - 4.0)))
    sol = bj.solve_composite(P, EPS, DP, MU, 1.0, 1.2)
    ua = np.concatenate([bj.profile(b2['yc'][:, i0], sol),
                         bj.profile(b1['yc'][:, i0], sol)])
    uc = np.concatenate([b2['u'][:, i0], b1['u'][:, i0]])
    umax = np.max(bj.profile(np.linspace(-1.2, 1.0, 2001), sol))
    e = (uc - ua) / umax
    # interface & plateau estimates
    dyf = b1['yc'][1, i0] - b1['yc'][0, i0]
    u_i = b1['u'][0, i0] - (b1['u'][1, i0] - b1['u'][0, i0]) * 0.5
    yp = b2['yc'][:, i0]
    mD = (yp > -0.9 * HP) & (yp < -0.25 * HP)
    uD_c = np.mean(b2['u'][mD, i0])
    u_top = b1['u'][-1, i0]          # near no-slip top wall
    u_c = np.max(b1['u'][:, i0])
    print('%s  L2=%.3e Linf=%.3e | u_i=%.6e (an %.6e) uD=%.6e (an %.6e) '
          'u_top=%.3e u_max_f=%.3e' %
          (tag, np.sqrt(np.mean(e ** 2)), np.max(np.abs(e)),
           u_i, sol['u_int'], uD_c, sol['uD'], u_top, u_c))
