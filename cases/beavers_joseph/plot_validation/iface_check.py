#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sanity checks at the fluid/porous interface of the converged solution:
- normal velocity on the interface row ~ 0 (flat interface, no mass flux)
- pressure of the fluid bottom row equals the porous top row (continuity)
- jump of cell-centre u across the interface equals ~ du/dy * dy (smooth)"""
import sys
import numpy as np
sys.path.insert(0, 'plot_validation')
import pv_vtk

b1 = pv_vtk.load_block('flow3d_block_1.vtk')
b2 = pv_vtk.load_block('flow3d_block_2.vtk')
i0 = int(np.argmin(np.abs(b1['xc'][0, :] - 4.0)))
print('fluid bottom row  (y=%.4f): v range [%.2e, %.2e], mean v=%.2e'
      % (b1['yc'][0, i0], b1['v'][0, :].min(), b1['v'][0, :].max(),
         b1['v'][0, :].mean()))
print('porous top row    (y=%.4f): v range [%.2e, %.2e], mean v=%.2e'
      % (b2['yc'][-1, i0], b2['v'][-1, :].min(), b2['v'][-1, :].max(),
         b2['v'][-1, :].mean()))
# pressure rows next to the interface
pf = b1['p'][0, :]
pp = b2['p'][-1, :]
print('p fluid  bottom row min/max: %.6e %.6e' % (pf.min(), pf.max()))
print('p porous top    row min/max: %.6e %.6e' % (pp.min(), pp.max()))
print('max|p_f - p_p| along interface = %.3e' % np.max(np.abs(pf - pp)))
# u continuity: cell centers are dy/2 apart across the interface
du = np.abs(b1['u'][0, i0] - b2['u'][-1, i0])
print('|u_f(y=+dy/2)-u_p(y=-dy/2)| at mid = %.3e  (dy=%.4f)'
      % (du, b1['yc'][1, i0] - b1['yc'][0, i0]))
