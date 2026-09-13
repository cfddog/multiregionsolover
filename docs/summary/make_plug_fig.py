#!/usr/bin/env python3
import os,sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'; F=R+'/docs/summary/figures'
sys.path.insert(0,R+'/cases/porous_plug/plot_validation'); import pv_io
mu=0.02; dp=0.10; H=1.0
res=[]
fields={}
for eps in (0.3,0.6,0.9):
    d=R+'/cases/porous_plug_ac/eps%02d'%int(round(eps*10))
    p,T,uvw=pv_io.load_vtk(d+'/flow3d_block_1.vtk'); u=uvw[:,:,0]; v=uvw[:,:,1]
    xs,ys=pv_io.cell_centers(); mag=np.hypot(u,v)
    K=pv_io.perme_ergun(eps,dp); lam=np.sqrt(K)/eps; alpha=H/(2*lam)
    imid=np.arange(10,34)
    coef=np.polyfit(xs[imid],p[10,imid],1); dpdx=coef[0]
    ubar=u.mean(axis=0)[imid].mean()
    G_an=ubar*mu/(lam**2*(1-np.tanh(alpha)/alpha))
    p_an=p[10,0]-G_an*(xs-xs[0])
    err=np.max(np.abs(p[10,:]-p_an))/max(1e-30,np.max(np.abs(p_an-p_an[-1])))
    res.append((eps,K,G_an,-dpdx,ubar,u.max(),np.max(np.abs(p[10,:]-p_an))))
    fields[eps]=(xs,ys,u,v,p,mag)
    print('eps=%.1f K=%.3e  dp/dx CFD=%.5e  analytic=%.5e  ubar=%.4e  rel.p-err=%.3e'%(eps,K,dpdx,-G_an,ubar,err))
# Fig 1: contours (u-mag, p)
fig,ax=plt.subplots(2,3,figsize=(13,6.6))
for k,eps in enumerate((0.3,0.6,0.9)):
    xs,ys,u,v,p,mag=fields[eps]; Xg,Yg=np.meshgrid(xs,ys)
    c0=ax[0,k].pcolormesh(Xg,Yg,mag,shading='auto',cmap='viridis'); plt.colorbar(c0,ax=ax[0,k],fraction=0.046)
    ax[0,k].set_title('|u| : eps=%.1f'%eps); ax[0,k].set_xlabel('x (m)'); ax[0,k].set_ylabel('y (m)')
    c1=ax[1,k].pcolormesh(Xg,Yg,p,shading='auto',cmap='coolwarm'); plt.colorbar(c1,ax=ax[1,k],fraction=0.046)
    ax[1,k].set_title('p : eps=%.1f'%eps); ax[1,k].set_xlabel('x (m)'); ax[1,k].set_ylabel('y (m)')
plt.tight_layout(); plt.savefig(F+'/plug_fields.png',dpi=140); print('plug_fields.png saved')
# Fig 2: pressure vs analytic
fig,ax=plt.subplots(1,3,figsize=(13,3.9))
for a,eps in zip(ax,(0.3,0.6,0.9)):
    xs,ys,u,v,p,mag=fields[eps]
    K=pv_io.perme_ergun(eps,dp); lam=np.sqrt(K)/eps; alpha=H/(2*lam); imid=np.arange(10,34)
    ubar=u.mean(axis=0)[imid].mean(); G_an=ubar*mu/(lam**2*(1-np.tanh(alpha)/alpha))
    p_an=p[10,0]-G_an*(xs-xs[0])
    a.plot(xs,p[10,:],'o',ms=3,mfc='none',label='CFD p (mid-row)')
    a.plot(xs,p_an,'k-',lw=1.3,label='analytic (Darcy-Brinkman)')
    a.set_title('eps=%.1f  dp/dx CFD=%.3e, an=%.3e'%(eps,np.polyfit(xs[imid],p[10,imid],1)[0],-G_an))
    a.set_xlabel('x (m)'); a.set_ylabel('p (Pa)'); a.legend(fontsize=8); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(F+'/plug_pressure.png',dpi=140); print('plug_pressure.png saved')
