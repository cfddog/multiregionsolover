#!/usr/bin/env python3
import sys, numpy as np, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ghia_u import GHIA
def load(d):
    lines=open(d+'/flow3d_block_1.vtk').read().splitlines()
    n0=[n for n,l in enumerate(lines) if l.startswith('POINTS')][0]; pts=int(lines[n0].split()[1])
    A=[]; i=n0+1
    while len(A)<pts: A.append([float(x) for x in lines[i].split()]); i+=1
    X=np.array(A)
    D={}; marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(lines) if ln.startswith(('SCALARS','VECTORS'))]
    cl=int([l for l in lines if l.startswith('CELL_DATA')][0].split()[1])
    for ii in range(len(marks)):
        a,typ,nm=marks[ii]; b=marks[ii+1][0] if ii+1<len(marks) else len(lines)
        body=[float(x) for l in lines[a+1:b] for x in l.split() if not l.startswith('LOOKUP')]
        D[nm]=np.array(body[:cl]) if typ=='SCALARS' else np.array(body[:3*cl]).reshape(-1,3)
    s=int(round(np.sqrt(pts/2))); n=s-1; Xn=X.reshape(2,s,s,3)
    Yc=0.5*(Xn[0,1:,0,1]+Xn[0,:-1,0,1]); Xc=0.5*(Xn[0,0,1:,0]+Xn[0,0,:-1,0])
    u=D['velocity'][:,0].reshape(n,n); return Xc,Yc,u
for d,Re in [(sys.argv[1],int(sys.argv[2]))]:
    Xc,Yc,u=load(d); g=GHIA[Re]; ic=int(np.argmin(np.abs(Xc-0.5)))
    up=np.interp(g['y'],Yc,u[:,ic]); eu=np.sqrt(np.mean((up-g['u'])**2))
    print('Re=%d  u-RMS=%.4f  u_min_sim=%.4f (Ghia %.4f)  n=%d'%(Re,eu,u[:,ic].min(),min(g['u']),u.shape[0]))
