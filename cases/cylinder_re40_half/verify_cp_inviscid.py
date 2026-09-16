#!/usr/bin/env python3
"""Inviscid-cylinder verification of the AC low-speed solver on the
half-cylinder (upper half + symmetry plane) grid.

The exact (potential-flow) solution for a circular cylinder in a uniform
stream of speed U is, on the surface (theta measured from the +x axis):

    u_theta(R,theta) = -2 U sin(theta)
    Cp(theta)        = 1 - 4 sin^2(theta)

so Cp = +1 at both the front and the rear stagnation points and Cp = -3 at
theta = 90 deg.  A correct inviscid solution must also give Cd = 0
(d'Alembert) and Cl = 0.

Usage:  python3 verify_cp_inviscid.py [dir] [U] [mode]
        mode = ac (default) or simple  (only for the report header)
"""
import sys, math
import numpy as np

R = 0.5
U = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
d = sys.argv[1] if len(sys.argv) > 1 else '.'

L = open(d + '/flow3d_block_1.vtk').read().splitlines()
def find(pre):
    for n, l in enumerate(L):
        if l.startswith(pre):
            return n
    raise SystemExit('missing ' + pre)
nx, ny, nz = (int(v) for v in L[find('DIMENSIONS')].split()[1:4])
p0 = find('POINTS'); npt = int(L[p0].split()[1])
X = np.array([[float(v) for v in l.split()] for l in L[p0 + 1:p0 + 1 + npt]]).reshape(nz, ny, nx, 3)
c0 = find('CELL_DATA'); ncl = int(L[c0].split()[1])
marks = [(n, l.split()[0], l.split()[1]) for n, l in enumerate(L)
         if l.startswith(('SCALARS', 'VECTORS'))]
dat = {}
for k, (n0, typ, nm) in enumerate(marks):
    n1 = marks[k + 1][0] if k + 1 < len(marks) else len(L)
    b = [float(v) for l in L[n0 + 1:n1] for v in l.split() if not l.startswith('LOOKUP')]
    dat[nm] = np.array(b[:ncl]) if typ == 'SCALARS' else np.array(b[:3 * ncl]).reshape(ncl, 3)
vel = dat['velocity'].reshape(nz - 1, ny - 1, nx - 1, 3)
pre = dat['pressure'].reshape(nz - 1, ny - 1, nx - 1)

rho = 1.0
q2 = vel[0, :, nx - 2, 0] ** 2 + vel[0, :, nx - 2, 1] ** 2
p_inf = (pre[0, :, nx - 2] + 0.5 * rho * (q2 - U ** 2)).mean()      # Bernoulli
p_inf_mean = pre[0, :, nx - 2].mean()

# wall face centres (node i=1 plane) and the first cell-centre radius
th = np.zeros(ny - 1); Cp = np.zeros(ny - 1); Cp_old = np.zeros(ny - 1)
rc1 = np.zeros(ny - 1)
for j in range(ny - 1):
    A = X[0, j, 0, :]; B = X[0, j + 1, 0, :]
    m = 0.5 * (A + B)
    th[j] = math.atan2(m[1], m[0])
    Cp[j] = (pre[0, j, 0] - p_inf) / (0.5 * rho * U ** 2)
    Cp_old[j] = (pre[0, j, 0] - p_inf_mean) / (0.5 * rho * U ** 2)
    rc1[j] = math.hypot(*[0.25 * (X[0, j, 0, i] + X[0, j + 1, 0, i]
                                 + X[0, j, 1, i] + X[0, j + 1, 1, i]) for i in (0, 1)])
Cpa = 1.0 - 4.0 * np.sin(th) ** 2

err = Cp - Cpa
print('---- inviscid cylinder Cp verification   grid %dx%dx%d  (half domain, R_out=%.1f D)'
      % (nx, ny, nz, rc1[0] * 0 + math.hypot(X[0, 0, nx - 1, 0], X[0, 0, nx - 1, 1])))
print(' p_inf (Bernoulli) = %+9.3e    p_inf (plain mean) = %+9.3e' % (p_inf, p_inf_mean))
print(' Cp error  max|dCp| = %.4f   RMS = %.4f   (Bernoulli p_inf)' % (np.abs(err).max(), np.sqrt((err ** 2).mean())))
print(' Cp error  max|dCp| = %.4f   RMS = %.4f   (plain-mean p_inf)'
      % (np.abs(Cp_old - Cpa).max(), np.sqrt(((Cp_old - Cpa) ** 2).mean())))
deg = np.degrees(th)
for a in (0.0, 10.0, 20.0, 45.0, 80.0, 90.0, 100.0, 135.0, 160.0, 170.0, 179.9):
    k = int(np.argmin(np.abs(deg - a)))
    print('   theta=%6.1f  Cp_num=%+8.4f  Cp_exact=%+8.4f  diff=%+8.4f'
          % (deg[k], Cp[k], Cpa[k], err[k]))
# surface speed (potential flow at the first cell centre r=rc1)
ut_num = -vel[0, :, 0, 0] * np.sin(th) + vel[0, :, 0, 1] * np.cos(th)
ut_ex = -U * (1.0 + (R / rc1) ** 2) * np.sin(th)
m = np.abs(ut_ex) > 0.2 * U
print(' surface speed u_t/(2U sin) : mean=%.4f  max dev=%.4f   (analytic 1.0 inside)'
      % ((ut_num[m] / ut_ex[m]).mean(), np.abs(ut_num[m] / ut_ex[m] - 1.0).max()))
# drag / lift from the pressure only (friction is zero by construction)
Fpx = -np.sum(Cp * np.cos(th)) * (math.pi / (ny - 1))
Fpy = -np.sum(Cp * np.sin(th)) * (math.pi / (ny - 1))
print(' pressure drag  Cd_p = %.4f   (exact 0)' % Fpx)
print(' lift coefficient of the UPPER HALF Cl = %.4f   (exact 1.6667 for potential flow;'
      ' full cylinder = 0)' % (0.5 * Fpy))
# momentum-balance drag on the outer arc (also exact 0)
Fm = 0.0
for j in range(ny - 1):
    A = X[0, j, nx - 1, :]; B = X[0, j + 1, nx - 1, :]
    dA = math.hypot(B[0] - A[0], B[1] - A[1])
    m2 = 0.5 * (A + B); t2 = math.atan2(m2[1], m2[0])
    n = (math.cos(t2), math.sin(t2))
    uc = vel[0, j, nx - 2, :]
    Fm += -(rho * uc[0] * (uc[0] * n[0] + uc[1] * n[1]) + pre[0, j, nx - 2] * n[0]) * dA
print(' momentum-balance Cd on outer arc = %.4f   (exact 0)' % (2.0 * Fm))
# potential flow has NO reverse flow: u >= 0 everywhere (along the rear ray
# u_r = U(1-R^2/r^2) rises monotonically from 0 at the rear stagnation point)
ur = vel[0, 0, :, 0]
print(' min u on the rear symmetry ray (j=1) = %.4f   (exact >=0)  cells with u<0: %d'
      % (ur.min(), int((vel[0, :, :, 0] < 0).sum())))

np.savetxt(d + '/cp_wall.dat',
           np.column_stack([deg, th / math.pi, Cp, Cp_old, Cpa]),
           header='theta_deg  theta/pi  Cp_num(Bernoulli)  Cp_num(mean p_inf)  Cp_exact=1-4sin^2',
           fmt='%12.6f')
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(6.2, 4.4))
    plt.plot(deg, Cpa, 'k-', lw=2, label=r'exact  $C_p=1-4\sin^2\theta$')
    plt.plot(deg, Cp, 'o', ms=3, mfc='none', color='C3', label='AC low-speed (inviscid)')
    plt.axhline(0, color='0.7', lw=0.6)
    plt.gca().invert_xaxis()
    plt.xlabel(r'$\theta$ (deg from rear stagnation point)')
    plt.ylabel(r'$C_p$'); plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(d + '/cp_inviscid.png', dpi=140)
    print(' wrote %s/cp_inviscid.png and %s/cp_wall.dat' % (d, d))
except Exception as e:
    print(' (matplotlib not available: %s)' % e)
