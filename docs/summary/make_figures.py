#!/usr/bin/env python3
import os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
root='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
sys.path.insert(0, root+'/docs/summary'); from ghia_u import GHIA
sys.path.insert(0, root+'/cases/porous_plug/plot_validation')
figs=root+'/docs/summary/figures'; os.makedirs(figs,exist_ok=True)
def read_vtk(p):
    L=open(p).read().splitlines(); D={}
    marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(L) if ln.startswith(('SCALARS','VECTORS'))]
    cl=int([l for l in L if l.startswith('CELL_DATA')][0].split()[1])
    for i in range(len(marks)):
        a,t,nm=marks[i]; b=marks[i+1][0] if i+1<len(marks) else len(L)
        body=[float(x) for l in L[a+1:b] for x in l.split() if not l.startswith('LOOKUP')]
        D[nm]=np.array(body[:cl]) if t=='SCALARS' else np.array(body[:3*cl]).reshape(-1,3)
    return D
def pts(p):
    L=open(p).read().splitlines(); n0=[n for n,l in enumerate(L) if l.startswith('POINTS')][0]
    npt=int(L[n0].split()[1]); A=[]; i=n0+1
    while len(A)<npt: A.append([float(x) for x in L[i].split()]); i+=1
    return np.array(A)
# cavity Re=100/1000/5000
fig,ax=plt.subplots(1,3,figsize=(12,4))
for k,(d,Re) in enumerate([('/tmp/run_cav_re100',100),('/tmp/run_cav_re1000',1000),(root+'/cases/lidcavity_ac/re5000',5000)]):
    if not os.path.exists(d+'/flow3d_block_1.vtk'): continue
    D=read_vtk(d+'/flow3d_block_1.vtk'); X=pts(d+'/flow3d_block_1.vtk')
    s=int(round(np.sqrt(X.shape[0]/2))); n=s-1; Xn=X.reshape(2,s,s,3)
    Yc=0.5*(Xn[0,1:,0,1]+Xn[0,:-1,0,1]); Xc=0.5*(Xn[0,0,1:,0]+Xn[0,0,:-1,0])
    u=D['velocity'][:,0].reshape(n,n); ic=int(np.argmin(np.abs(Xc-0.5)))
    g=GHIA[Re]; ax[k].plot(u[:,ic],Yc,'o',ms=3,label='AC CFD'); ax[k].plot(g['u'],g['y'],'k-',lw=1.2,label='Ghia(1982)')
    ax[k].set_title('cavity Re=%d'%Re); ax[k].set_xlabel('u(x=0.5)'); ax[k].set_ylabel('y'); ax[k].legend(fontsize=8); ax[k].grid(alpha=.3)
plt.tight_layout(); plt.savefig(figs+'/cavity.png',dpi=140); print('cavity.png')
# channel
D=read_vtk('/tmp/run_chan/flow3d_block_1.vtk'); X=pts('/tmp/run_chan/flow3d_block_1.vtk')
Xr=X.reshape(2,41,161,3); Yc=0.5*(Xr[0,1:,:,1]+Xr[0,:-1,:,1])[:,0]; u=D['velocity'][:,0].reshape(40,160)
plt.figure(figsize=(4.4,4.2)); plt.plot(u[:,120],Yc,'o',ms=3,label='AC CFD'); yy=np.linspace(0,1,200); plt.plot(6*yy*(1-yy),yy,'k-',lw=1.2,label='Poiseuille')
plt.xlabel('u'); plt.ylabel('y'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/channel.png',dpi=140); print('channel.png')
# solids
D=np.loadtxt('/tmp/run_solid_1d/Ts_block_1.dat',skiprows=3); x,T=D[:,0],D[:,3]
xu,inv=np.unique(np.round(x,10),return_inverse=True); Tx=np.array([T[inv==k].mean() for k in range(xu.size)])
plt.figure(figsize=(4.4,4.2)); plt.plot(xu,Tx,'o',ms=3,label='FVM'); plt.plot(xu,300+(800/15.1)*xu,'k-',lw=1.2,label='analytic'); plt.xlabel('x (m)'); plt.ylabel('T (K)'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/solid1d.png',dpi=140); print('solid1d.png')
D=np.loadtxt('/tmp/run_solid_annulus/Ts_block_1.dat',skiprows=3); X,Y,T=D[:,0],D[:,1],D[:,3]; r=np.hypot(X,Y)
ii_=np.arange(D.shape[0])%40; ri=np.array([r[ii_==k].mean() for k in range(40)]); Ti=np.array([T[ii_==k].mean() for k in range(40)])
plt.figure(figsize=(4.4,4.2)); plt.plot(ri,Ti,'o',ms=3,label='FVM'); plt.plot(ri,300+(800*0.1/15.1)*np.log(0.2/ri),'k-',lw=1.2,label='analytic'); plt.xlabel('r (m)'); plt.ylabel('T (K)'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/annulus.png',dpi=140); print('annulus.png')
# plug: u(y) for 3 eps vs Brinkman
import pv_io
plt.figure(figsize=(4.6,4.2))
for eps in (0.3,0.6,0.9):
    d=root+'/cases/porous_plug_ac/eps%02d'%int(round(eps*10))
    if not os.path.exists(d+'/flow3d_block_1.vtk'): continue
    p,T,uvw=pv_io.load_vtk(d+'/flow3d_block_1.vtk'); uu=uvw[:,:,0]; xs,ys=pv_io.cell_centers()
    coef=np.polyfit(xs[10:34],p[10,10:34],1); ub,uD,K=pv_io.brinkman_profile(ys,coef[0],eps,0.10,0.02)
    plt.plot(uu[:,20],ys,'o',ms=3); plt.plot(ub,ys,'-',lw=1,label='eps=%.1f'%eps)
plt.xlabel('u'); plt.ylabel('y'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/plug.png',dpi=140); print('plug.png')
# LTNE one case
import importlib.util
spec=importlib.util.spec_from_file_location('vl',root+'/docs/summary/validate_ltne_param.py')
# reuse readers inline
def read_flow_t(p):
    L=open(p).read().splitlines(); v=[]; i=0
    while i<len(L):
        if L[i].strip().startswith('SCALARS') and L[i].split()[1]=='temperature':
            i+=2
            while i<len(L) and not L[i].strip().startswith(('SCALARS','VECTORS','POINT_DATA','CELL_DATA','FIELD')): v+=L[i].split(); i+=1
            break
        i+=1
    return np.array(v,float)
d=root+'/cases/porous_ltne_1d/sweep/eps03_G2_q2'
if os.path.exists(d+'/Ts_block_1.dat'):
    dd=np.loadtxt(d+'/Ts_block_1.dat',skiprows=3); xx,Ts=dd[:,0],dd[:,3]; Tf=read_flow_t(d+'/flow3d_block_1.vtk')[:xx.size]
    eps,G,q,hv=0.3,2.0,2e6,2.013e6; L=0.01; ks=13.4; cp=1007.0; hc=31.4; Tinf=300.; Tin=300.; kse=(1-eps)*ks
    A=hv/(G*cp); B=hv/kse; beta=kse/(G*cp); Ds=np.sqrt(A*A+4*B); l1=(-A+Ds)/2; l2=(-A-Ds)/2
    e1,e2=np.exp(l1*xx),np.exp(l2*xx); e1L,e2L=np.exp(l1*L),np.exp(l2*L)
    M=np.array([[1,beta*l1,beta*l2],[-hc,kse*l1-hc,kse*l2-hc],[0,kse*l1*e1L,kse*l2*e2L]])
    C0,C1,C2=np.linalg.solve(M,np.array([Tin,-hc*Tinf,q])); Tsa=C0+C1*e1+C2*e2; Tfa=C0+beta*(C1*l1*e1+C2*l2*e2)
    plt.figure(figsize=(4.6,4.2)); plt.plot(xx*1e3,Tf,'s',ms=2.5,label='Tf CFD'); plt.plot(xx*1e3,Tfa,'b-',lw=1.2,label='Tf analytic'); plt.plot(xx*1e3,Ts,'o',ms=2.5,label='Ts CFD'); plt.plot(xx*1e3,Tsa,'r-',lw=1.2,label='Ts analytic'); plt.xlabel('y (mm)'); plt.ylabel('T (K)'); plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(figs+'/ltne.png',dpi=140); print('ltne.png')
print('done')
