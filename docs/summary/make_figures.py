#!/usr/bin/env python3
import os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
root='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
sys.path.insert(0, root+'/docs/summary'); from ghia_u import GHIA
figs=root+'/docs/summary/figures'; os.makedirs(figs,exist_ok=True)

def read_vtk(path):
    lines=open(path).read().splitlines(); D={}; marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(lines) if ln.startswith(('SCALARS','VECTORS'))]
    cd=[l for l in lines if l.startswith('CELL_DATA')][0].split()[1]; cl=int(cd)
    for ii in range(len(marks)):
        n0,typ,nm=marks[ii]; n1=marks[ii+1][0] if ii+1<len(marks) else len(lines)
        body=[float(x) for l in lines[n0+1:n1] for x in l.split() if not l.startswith('LOOKUP')]
        D[nm]=np.array(body[:cl]) if typ=='SCALARS' else np.array(body[:3*cl]).reshape(-1,3)
    pts=[float(x) for l in lines[[n for n,l in enumerate(lines) if l.startswith('POINTS')][0]+1:] for x in l.split()][: ] if False else None
    return D

def pts_xyz(path):
    lines=open(path).read().splitlines()
    n0=[n for n,l in enumerate(lines) if l.startswith('POINTS')][0]
    pts=int(lines[n0].split()[1]); A=[]
    i=n0+1
    while len(A)<pts: A.append([float(x) for x in lines[i].split()]); i+=1
    return np.array(A)

# ---- cavity Re=100/1000 (AC) vs Ghia ----
fig,ax=plt.subplots(1,2,figsize=(9,4.2))
for d,Re,ic in [('/tmp/run_cav_re100',100,19),('/tmp/run_cav_re1000',1000,39)]:
    D=read_vtk(d+'/flow3d_block_1.vtk'); X=pts_xyz(d+'/flow3d_block_1.vtk')
    s=int(round(np.sqrt(X.shape[0]/2))); n=s-1
    Xn=X.reshape(2,s,s,3)
    if Re==100: Yc=0.5*(Xn[0,1:,0,1]+Xn[0,:-1,0,1]); ic=19
    else:       Yc=0.5*(Xn[0,1:,0,1]+Xn[0,:-1,0,1]); ic=39
    u=D['velocity'][:,0].reshape(n,n)
    ax[0 if Re==100 else 1].plot(u[:,ic],Yc,'o',ms=3,label='AC CFD'); g=GHIA[Re]
    ax[0 if Re==100 else 1].plot(g['u'],g['y'],'k-',lw=1.2,label='Ghia(1982)')
    ax[0 if Re==100 else 1].set_title('cavity Re=%d'%Re); ax[0 if Re==100 else 1].set_xlabel('u at x=0.5'); ax[0 if Re==100 else 1].set_ylabel('y'); ax[0 if Re==100 else 1].legend(fontsize=8); ax[0 if Re==100 else 1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(figs+'/cavity.png',dpi=140); print('cavity.png')

# ---- channel vs Poiseuille ----
D=read_vtk('/tmp/run_chan/flow3d_block_1.vtk'); X=pts_xyz('/tmp/run_chan/flow3d_block_1.vtk')
Xr=X.reshape(2,41,161,3); Yc=0.5*(Xr[0,1:,:,1]+Xr[0,:-1,:,1])[:,0]; u=D['velocity'][:,0].reshape(40,160)
plt.figure(figsize=(4.4,4.2)); plt.plot(u[:,120],Yc,'o',ms=3,label='AC CFD (x=6)'); yy=np.linspace(0,1,200); plt.plot(6*yy*(1-yy),yy,'k-',lw=1.2,label='Poiseuille')
plt.xlabel('u'); plt.ylabel('y'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/channel.png',dpi=140); print('channel.png')

# ---- solid 1d + annulus ----
D=np.loadtxt('/tmp/run_solid_1d/Ts_block_1.dat',skiprows=3); x,T=D[:,0],D[:,3]
xu,inv=np.unique(np.round(x,10),return_inverse=True); Tx=np.array([T[inv==k].mean() for k in range(xu.size)])
plt.figure(figsize=(4.4,4.2)); plt.plot(xu,Tx,'o',ms=3,label='FVM'); plt.plot(xu,300+(800/15.1)*xu,'k-',lw=1.2,label='analytic'); plt.xlabel('x (m)'); plt.ylabel('T (K)'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/solid1d.png',dpi=140); print('solid1d.png')

D=np.loadtxt('/tmp/run_solid_annulus/Ts_block_1.dat',skiprows=3); X,Y,T=d[:,0] if False else (D[:,0],D[:,1],D[:,3]); r=np.hypot(X,Y)
i_idx=np.arange(D.shape[0])%40; ri=np.array([r[i_idx==k].mean() for k in range(40)]); Ti=np.array([T[i_idx==k].mean() for k in range(40)])
plt.figure(figsize=(4.4,4.2)); plt.plot(ri,Ti,'o',ms=3,label='FVM'); plt.plot(ri,300+(800*0.1/15.1)*np.log(0.2/ri),'k-',lw=1.2,label='analytic'); plt.xlabel('r (m)'); plt.ylabel('T (K)'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/annulus.png',dpi=140); print('annulus.png')
print('done')
