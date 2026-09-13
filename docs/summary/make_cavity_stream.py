#!/usr/bin/env python3
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'; F=R+'/docs/summary/figures'
def load(d):
    L=open(d+'/flow3d_block_1.vtk').read().splitlines()
    n0=[n for n,l in enumerate(L) if l.startswith('POINTS')][0]; pts=int(L[n0].split()[1]); A=[]; i=n0+1
    while len(A)<pts: A.append([float(x) for x in L[i].split()]); i+=1
    X=np.array(A)
    D={}; marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(L) if ln.startswith(('SCALARS','VECTORS'))]
    cl=int([l for l in L if l.startswith('CELL_DATA')][0].split()[1])
    for ii in range(len(marks)):
        a,t,nm=marks[ii]; b=marks[ii+1][0] if ii+1<len(marks) else len(L)
        body=[float(x) for l in L[a+1:b] for x in l.split() if not l.startswith('LOOKUP')]
        D[nm]=np.array(body[:cl]) if t=='SCALARS' else np.array(body[:3*cl]).reshape(-1,3)
    s=int(round(np.sqrt(pts/2))); n=s-1; Xn=X.reshape(2,s,s,3)
    Xc=0.5*(Xn[0,0,1:,0]+Xn[0,0,:-1,0]); Yc=0.5*(Xn[0,1:,0,1]+Xn[0,:-1,0,1])
    U=D['velocity'][:,0].reshape(n,n); V=D['velocity'][:,1].reshape(n,n)   # [j,i]
    return Xc,Yc,U,V
cases=[('/tmp/run_cav_re100',100),('/tmp/run_cav_re1000',1000),(R+'/cases/lidcavity_ac/re5000',5000)]
fig,ax=plt.subplots(1,3,figsize=(13,4.2))
for a,(d,Re) in zip(ax,cases):
    Xc,Yc,U,V=load(d); Xg,Yg=np.meshgrid(Xc,Yc); sp=np.sqrt(U*U+V*V)
    lw=2.5*sp/sp.max()
    a.streamplot(Xg,Yg,U,V,density=1.3,color=sp,cmap='viridis',linewidth=lw)
    a.set_title('cavity Re=%d streamlines'%Re); a.set_xlabel('x'); a.set_ylabel('y'); a.set_aspect('equal'); a.set_xlim(0,1); a.set_ylim(0,1)
plt.tight_layout(); plt.savefig(F+'/cavity_streamlines.png',dpi=140); print('cavity_streamlines.png saved')
for a,(d,Re) in [(None,(x[0],x[1])) for x in cases]:
    Xc,Yc,U,V=load(d); print('Re=%d  umax=%.3f vmax=%.3f'%(Re,np.abs(U).max(),np.abs(V).max()))
