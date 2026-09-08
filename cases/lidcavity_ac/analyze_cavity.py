#!/usr/bin/env python3
"""Lid-driven cavity AC-solver check: vertical u(x=0.5) and horizontal v(y=0.5)
vs Ghia et al. (1982) Tables I/II for Re=100 and Re=1000.
Usage: python3 analyze_cavity.py <dir> [<dir> ...]
"""
import sys
import numpy as np

GH = {
100: {
 'y': [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,
       0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000],
 'u': [1.00000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,
       -0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717,0.0]},
1000: {
 'y': [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,
       0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000],
 'u': [1.00000,0.65928,0.57492,0.51117,0.46604,0.33304,0.18719,0.05702,
       -0.06080,-0.10648,-0.27805,-0.38289,-0.29730,-0.22220,-0.20196,-0.18109,0.0],
 'x': [0.0,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5000,0.6172,0.7344,
       0.8516,0.9531,0.9609,0.9688,0.9766,1.0],
 'v': [0.0,0.09245,0.10091,0.15444,0.17570,0.17527,0.13690,0.06080,
       -0.05702,-0.18719,-0.33304,-0.46604,-0.51117,-0.57492,-0.65928,0.0]},
}

def load(dirn):
    lines = open(dirn+'/flow3d_block_1.vtk').read().splitlines()
    def find(pre):
        for n,ln in enumerate(lines):
            if ln.startswith(pre): return n
        raise SystemExit('nf')
    p0=find('POINTS'); pts=int(lines[p0].split()[1])
    XYZ=np.array([[float(x) for x in l.split()] for l in lines[p0+1:p0+1+pts]])
    cd=find('CELL_DATA'); cl=int(lines[cd].split()[1])
    marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(lines)
           if ln.startswith(('SCALARS','VECTORS'))]
    D={}
    for ii in range(len(marks)):
        n0,typ,nm=marks[ii]; n1=marks[ii+1][0] if ii+1<len(marks) else len(lines)
        body=[float(x) for l in lines[n0+1:n1] for x in l.split() if not l.startswith('LOOKUP')]
        if typ=='SCALARS': D[nm]=np.array(body[:cl])
        else: D[nm]=np.array(body[:3*cl]).reshape(-1,3)
    s = int(round(np.sqrt(pts/2)))          # nodes per side
    n = s-1
    Xn = XYZ.reshape(2,s,s,3)
    Yc = 0.5*(Xn[0,1:,0,1]+Xn[0,:-1,0,1])
    Xc = 0.5*(Xn[0,0,1:,0]+Xn[0,0,:-1,0])
    u = D['velocity'][:,0].reshape(n,n)      # [j,i]
    v = D['velocity'][:,1].reshape(n,n)
    return Xc, Yc, u, v

def main():
    for d in sys.argv[1:]:
        Xc,Yc,u,v = load(d.rstrip('/'))
        # Re from LS_mu in control
        mu=None
        for ln in open(d+'/control.ec'):
            if 'LS_mu' in ln and not ln.strip().startswith('!'):
                try: mu=float(ln.split('=')[1].split('d')[0].split('D')[0])
                except Exception: pass
        Re = int(round(1.0/ (mu if mu else 0.01)))
        g=GH[Re]
        ic=np.argmin(np.abs(Xc-0.5)); jc=np.argmin(np.abs(Yc-0.5))
        up = np.interp(g['y'], Yc, u[:,ic])
        eu=np.sqrt(np.mean((up-g['u'])**2))
        if 'v' in g:
            vp = np.interp(g['x'], Xc, v[jc,:])
            ev=np.sqrt(np.mean((vp-g['v'])**2))
        else:
            vp=None; ev=float('nan')
        print(f'{d:35s} Re={Re}  u-RMS={eu:.5f} v-RMS={ev:.5f}   u_min={u[:,ic].min():+.5f}')
        print('    vertical u(x=0.5):  y / sim / Ghia')
        for k in [1,3,5,7,9,10,12,14]:
            print(f'      {g["y"][k]:.4f}  {up[k]:+.5f} {g["u"][k]:+.5f}')
        if vp is not None:
            print('    horizontal v(y=0.5):  x / sim / Ghia')
            for k in [1,3,5,7,8,10,12,14]:
                print(f'      {g["x"][k]:.4f}  {vp[k]:+.5f} {g["v"][k]:+.5f}')

if __name__=='__main__':
    main()
