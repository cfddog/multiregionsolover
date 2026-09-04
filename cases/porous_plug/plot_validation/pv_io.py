import numpy as np

NXC, NYC = 40, 20          # interior cells (i=0..39, j=0..19)
LX, LY = 8.0, 1.0
DX = LX / NXC
DY = LY / NYC

def cell_centers():
    xs = (np.arange(NXC) + 0.5) * DX
    ys = (np.arange(NYC) + 0.5) * DY
    return xs, ys

def load_vtk(path):
    lines = open(path).read().splitlines()
    def read_scalar(name):
        i = lines.index('SCALARS %s double 1' % name)
        j = lines.index('LOOKUP_TABLE default', i)
        v = np.array([float(s.split()[0]) for s in lines[j+1:j+1+NXC*NYC]])
        return v.reshape((NYC, NXC))
    def read_vectors(name):
        i = lines.index('VECTORS %s double' % name)
        v = np.zeros((NYC, NXC, 3))
        for n, s in enumerate(lines[i+1:i+1+NXC*NYC]):
            t = s.split()
            j = n // NXC; i2 = n % NXC
            v[j, i2, 0] = float(t[0]); v[j, i2, 1] = float(t[1]); v[j, i2, 2] = float(t[2])
        return v
    p = read_scalar('pressure')
    T = read_scalar('temperature')
    uvw = read_vectors('velocity')
    return p, T, uvw

def load_ts(path):
    vals = []
    start = False
    for l in open(path):
        if 'DATAPACKING' in l:
            start = True; continue
        if start:
            t = l.split()
            if len(t) >= 4:
                vals.append(float(t[3]))
    a = np.array(vals)
    if a.size != NXC*NYC:
        a = a[:NXC*NYC]
    return a.reshape((NYC, NXC))

def perme_ergun(eps, dp):
    return dp**2 * eps**3 / (150.0*(1.0-eps)**2)

def brinkman_profile(y, dpdx, eps, dp, mu):
    """Steady model solved by the code (uniform, no Forchheimer, no conv):
       mu u'' - (eps^2 mu / K) u = dp/dx ,  u=0 at y=0,H.
    """
    K = perme_ergun(eps, dp)
    lam = np.sqrt(K)/eps
    uD = -dpdx * K/(eps*eps*mu)
    return uD*(1.0 - np.cosh((y-0.5*LY)/lam)/np.cosh(0.5*LY/lam)), uD, K
