#!/usr/bin/env python3
"""Quick wall diagnostics: Cp(theta) table + maximum speed near the wall,
potential-flow Cp for reference.  Usage: diag_wall.py [dir] [Re]"""
import sys, math
import numpy as np
R = 0.5; rho = 1.0; U = 1.0
d = sys.argv[1] if len(sys.argv) > 1 else '.'
Re = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
L = open(d + '/flow3d_block_1.vtk').read().splitlines()
def find(p):
    for n, l in enumerate(L):
        if l.startswith(p): return n
    raise SystemExit('missing ' + p)
nx, ny, nz = (int(v) for v in L[find('DIMENSIONS')].split()[1:4])
p0 = find('POINTS'); npt = int(L[p0].split()[1])
X = np.array([[float(v) for v in l.split()] for l in L[p0+1:p0+1+npt]]).reshape(nz, ny, nx, 3)
c0 = find('CELL_DATA'); ncl = int(L[c0].split()[1])
marks = [(n, l.split()[0], l.split()[1]) for n, l in enumerate(L)
         if l.startswith(('SCALARS', 'VECTORS'))]
D = {}
for k, (n0, typ, nm) in enumerate(marks):
    n1 = marks[k+1][0] if k+1 < len(marks) else len(L)
    b = [float(v) for l in L[n0+1:n1] for v in l.split() if not l.startswith('LOOKUP')]
    D[nm] = np.array(b[:ncl]) if typ == 'SCALARS' else np.array(b[:3*ncl]).reshape(ncl, 3)
vel = D['velocity'].reshape(nz-1, ny-1, nx-1, 3); pre = D['pressure'].reshape(nz-1, ny-1, nx-1)
q2 = vel[0,:,nx-2,0]**2 + vel[0,:,nx-2,1]**2
p_bern = (pre[0,:,nx-2] + 0.5*rho*(q2 - U**2)).mean()
p_mean = pre[0,:,nx-2].mean()
print('p_inf: Bernoulli=%+9.3e  plain mean=%+9.3e' % (p_bern, p_mean))
print('max |u| on i=1 layer = %.4f   global max |u| = %.4f' %
      (np.hypot(vel[0,:,0,0], vel[0,:,0,1]).max(), np.hypot(vel[0,:,:,0], vel[0,:,:,1]).max()))
print(' theta   Cp(bern)   Cp(mean)   Cp_pot')
for a in (2, 10, 20, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 178):
    j = int(np.argmin([abs(180.0*jj/(ny-1) - a) for jj in range(ny)]))
    th = math.pi*(j+0.5)/(ny-1)
    print('%6.1f  %+8.4f  %+8.4f  %+8.4f' % (math.degrees(th),
          (pre[0,j,0]-p_bern)/(0.5*rho), (pre[0,j,0]-p_mean)/(0.5*rho), 1-4*math.sin(th)**2))
