#!/usr/bin/env python3
"""Post-process the half-cylinder (upper half + symmetry) Re=40 solution.

Reads flow3d_block_1.vtk (cell-centred rho,u,v,w,T,p) and reports
  * drag coefficient Cd = Cd_pressure + Cd_friction   (wall surface integral)
  * an INDEPENDENT Cd from the far-field momentum balance on the outer arc
  * separation angle, recirculation (wake) length, Cp_stag, Cp_min, Cl

Non-dimensionalisation used by the case: D=1, rho=1, U=1, mu=1/Re
  Cd = F_x / (0.5 rho U^2 D * 1) ;  only the upper half is meshed, so the
  half-domain force is doubled (symmetry) -> Cd = 4*F_half for D=rho=U=1.

Wall force (per unit depth), n = outward normal of the body (= e_r):
  dF_x = [-p_w cos(theta) - mu (du_t/dr)|_w sin(theta)] dA
with u_t = u.e_theta = -u sin(theta) + v cos(theta),  t = e_theta.

Far-field balance for the half domain (the two symmetry rays give no
x-momentum contribution because v=0 and du/dy=dv/dx=0 on them):
  F_half,x = -oint_arc [ rho u_x (u.n) + p n_x - mu du_x/dn ] dA
"""
import sys, math
import numpy as np

D = 1.0; R = 0.5 * D
rho = 1.0; U = 1.0
Re = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
mu = rho * U * D / Re
fn = (sys.argv[1] if len(sys.argv) > 1 else '.') + '/flow3d_block_1.vtk'


def read_vtk(name):
    L = open(name).read().splitlines()
    def find(pre):
        for n, l in enumerate(L):
            if l.startswith(pre):
                return n
        raise SystemExit('missing ' + pre)
    nx, ny, nz = (int(v) for v in L[find('DIMENSIONS')].split()[1:4])
    p0 = find('POINTS'); npt = int(L[p0].split()[1])
    XYZ = np.array([[float(v) for v in l.split()] for l in L[p0 + 1:p0 + 1 + npt]])
    XYZ = XYZ.reshape(nz, ny, nx, 3)
    c0 = find('CELL_DATA'); ncl = int(L[c0].split()[1])
    marks = [(n, l.split()[0], l.split()[1]) for n, l in enumerate(L)
             if l.startswith(('SCALARS', 'VECTORS'))]
    dat = {}
    for k, (n0, typ, nm) in enumerate(marks):
        n1 = marks[k + 1][0] if k + 1 < len(marks) else len(L)
        body = [float(v) for l in L[n0 + 1:n1] for v in l.split() if not l.startswith('LOOKUP')]
        dat[nm] = np.array(body[:ncl]) if typ == 'SCALARS' else np.array(body[:3 * ncl]).reshape(ncl, 3)
    return nx, ny, nz, XYZ, dat


nx, ny, nz, X, dat = read_vtk(fn)
vel = dat['velocity'].reshape(nz - 1, ny - 1, nx - 1, 3)
pre = dat['pressure'].reshape(nz - 1, ny - 1, nx - 1)

# cell-centre coordinates = average of the cell corners
xc = 0.125 * (X[0, :-1, :-1, 0] + X[0, 1:, :-1, 0] + X[0, :-1, 1:, 0] + X[0, 1:, 1:, 0]
              + X[1, :-1, :-1, 0] + X[1, 1:, :-1, 0] + X[1, :-1, 1:, 0] + X[1, 1:, 1:, 0])
yc = 0.125 * (X[0, :-1, :-1, 1] + X[0, 1:, :-1, 1] + X[0, :-1, 1:, 1] + X[0, 1:, 1:, 1]
              + X[1, :-1, :-1, 1] + X[1, 1:, :-1, 1] + X[1, :-1, 1:, 1] + X[1, 1:, 1:, 1])
rc = np.sqrt(xc ** 2 + yc ** 2)
tc = np.arctan2(yc, xc)                        # cell-centre theta in [0,pi]

# ---------------------------------------------------------------- wall force
Fp = Ff = Fp0 = Ff0 = 0.0
tau = np.zeros(ny - 1); th_w = np.zeros(ny - 1); pw = np.zeros(ny - 1)
for j in range(ny - 1):
    A = X[0, j, 0, :]; B = X[0, j + 1, 0, :]
    dA = math.hypot(B[0] - A[0], B[1] - A[1])          # times dz = 1
    mid = 0.5 * (A + B); th = math.atan2(mid[1], mid[0])
    d1 = rc[j, 0] - R; d2 = rc[j, 1] - R
    th_w[j] = th
    p0 = pre[0, j, 0]; p1 = pre[0, j, 1]
    pw[j] = p0 - d1 * (p1 - p0) / (d2 - d1)            # linear extrapolation
    u1 = -vel[0, j, 0, 0] * math.sin(tc[j, 0]) + vel[0, j, 0, 1] * math.cos(tc[j, 0])
    u2 = -vel[0, j, 1, 0] * math.sin(tc[j, 1]) + vel[0, j, 1, 1] * math.cos(tc[j, 1])
    a2 = (u1 * d2 ** 2 - u2 * d1 ** 2) / (d1 * d2 * (d2 - d1))   # u_t(0)=0 fit
    a1 = u1 / d1
    tau[j] = mu * a2
    Fp += -pw[j] * math.cos(th) * dA
    Ff += -mu * a2 * math.sin(th) * dA
    Fp0 += -p0 * math.cos(th) * dA
    Ff0 += -mu * a1 * math.sin(th) * dA

Cd_wall = 4.0 * (Fp + Ff)                              # D = rho = U = 1
Cd_p = 4.0 * Fp; Cd_f = 4.0 * Ff
Cd_wall_1st = 4.0 * (Fp0 + Ff0)

# ------------------------------------------------- far-field momentum balance
# The drag on the body can be obtained from the x-momentum balance on ANY
# closed surface around it.  For the half domain the two symmetry rays
# contribute nothing (v = du/dy = dv/dx = 0 there), so
#   F_half,x = -oint_arc [ rho u_x (u.n) + p n_x - mu du_x/dn ] dA
# Evaluating this on several radii is a strong consistency check.
def cd_mom(ico):
    """control surface = grid plane i=ico (node index); cells inside: ico-1"""
    F = 0.0
    for j in range(ny - 1):
        A = X[0, j, ico, :]; B = X[0, j + 1, ico, :]
        dA = math.hypot(B[0] - A[0], B[1] - A[1])
        mid = 0.5 * (A + B); th = math.atan2(mid[1], mid[0])
        nx1 = math.cos(th); ny1 = math.sin(th)
        uc = vel[0, j, ico - 1, :]
        un = uc[0] * nx1 + uc[1] * ny1
        dr = rc[j, ico - 1] - rc[j, ico - 2]
        dux = (vel[0, j, ico - 1, 0] - vel[0, j, ico - 2, 0]) / dr
        F += -(rho * uc[0] * un + pre[0, j, ico - 1] * nx1 - mu * dux) * dA
    return 4.0 * F

ic = nx - 2                                    # outermost cell layer
Cd_mom = cd_mom(nx - 1)                        # outer boundary
# interior control surfaces (r ~ 2, 5, 12 D) for the consistency check
r_target = [2.0, 5.0, 12.0]
Cd_arc = []
for rt in r_target:
    ico = int(np.argmin([abs(rc[0, i] - rt) for i in range(nx - 2)])) + 1
    Cd_arc.append((rc[0, ico - 1], cd_mom(ico)))

# ------------------------------------------------ free stream / Cp, Cl, wake
p_inf = pre[0, :, ic].mean()
Cp_stag = (pre[0, ny - 2, 0] - p_inf) / (0.5 * rho * U ** 2)
Cp_min = (pre[0, :, 0].min() - p_inf) / (0.5 * rho * U ** 2)
Cl = 0.0
for j in range(ny - 1):
    A = X[0, j, 0, :]; B = X[0, j + 1, 0, :]
    dA = math.hypot(B[0] - A[0], B[1] - A[1])
    mid = 0.5 * (A + B); th = math.atan2(mid[1], mid[0])
    Cl += 2.0 * (-pw[j] * math.sin(th) + tau[j] * math.cos(th)) * dA

# separation : last theta (walking back from the rear) where tau changes sign
sep = None
for j in range(ny - 3, 0, -1):
    if tau[j] * tau[j + 1] < 0:
        fr = tau[j + 1] / (tau[j + 1] - tau[j])
        th_s = th_w[j + 1] + fr * (th_w[j] - th_w[j + 1])
        sep = math.degrees(math.pi - th_s)
        break

# wake length on the symmetry ray y = 0 (cell row j = 0)
xr = rc[0, :] * np.cos(tc[0, :])
ur = vel[0, 0, :, 0]
Lw = None
for i in range(1, nx - 1):
    if ur[i - 1] < 0 <= ur[i]:
        Lw = xr[i] - R
        break

print('---- cylinder (upper half + symmetry), Re=%.1f   grid %dx%dx%d' % (Re, nx, ny, nz))
print(' Cd (wall integral, 2nd order) = %.4f   [pressure %.4f + friction %.4f]'
      % (Cd_wall, Cd_p, Cd_f))
print(' Cd (wall integral, 1st order) = %.4f' % Cd_wall_1st)
print(' Cd (far-field momentum bal., r=%.1f) = %.4f' % (rc[0, ic], Cd_mom))
for rr, cc in Cd_arc:
    print(' Cd (momentum balance on arc r=%5.1f) = %.4f' % (rr, cc))
print(' lift coefficient of the upper half Cl = %.4f (full cylinder = 0 by symmetry)'
      % (0.5 * Cl))
print(' Cp_stagnation = %.4f   Cp_min = %.4f   -Cp_rear = %.4f   p_inf = %.3e'
      % (Cp_stag, Cp_min, -(pre[0, 0, 0] - p_inf) / (0.5 * rho * U ** 2), p_inf))
print(' separation angle (from front stag.) = %s deg' % ('%.1f' % sep if sep else 'none'))
print(' wake length (reattachment, from rear) = %s D   u_min(centerline) = %.5f'
      % ('%.3f' % Lw if Lw else 'none', ur.min()))
