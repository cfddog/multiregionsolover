#!/usr/bin/env python3
import os,sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'; F=R+'/docs/summary/figures'
sys.path.insert(0,R+'/cases/beavers_joseph/plot_validation'); import pv_vtk, bj_analytical as bj
HP=1.2; H=1.0; DP=0.5; MU=0.02
fig,ax=plt.subplots(1,3,figsize=(13,4.0))
for a,eps in zip(ax,(0.3,0.6,0.9)):
    d=R+'/cases/beavers_joseph/eps%02d'%int(round(eps*10))
    b1=pv_vtk.load_block(d+'/flow3d_block_1.vtk'); b2=pv_vtk.load_block(d+'/flow3d_block_2.vtk')
    x=b1['xc'][0]; dx=x[1]-x[0]; L=x[-1]-x[0]+dx; xw=(x>0.25*L)&(x<0.75*L)
    jf=b1['p'].shape[0]//2; jp=b2['p'].shape[0]//2
    G=-0.5*(np.polyfit(x[xw],b1['p'][jf,xw],1)[0]+np.polyfit(x[xw],b2['p'][jp,xw],1)[0])
    sol=bj.solve_composite(G,eps,DP,MU,h=H,Hp=HP)
    i0=int(np.argmin(np.abs(x-0.5*L)))
    yc=np.concatenate([b2['yc'][:,i0],b1['yc'][:,i0]]); uc=np.concatenate([b2['u'][:,i0],b1['u'][:,i0]])
    ya=np.linspace(-HP,H,600)
    a.plot(uc,yc,'o',ms=3,mfc='none',label='CFD (x=L/2)')
    a.plot(bj.profile(ya,sol),ya,'k-',lw=1.3,label='analytic (composite)')
    a.plot(bj.u_poiseuille_noslip(np.linspace(0,H,200),G,MU,H),np.linspace(0,H,200),':',color='gray',label='no-slip Poiseuille')
    a.axhline(0,color='k',lw=0.5)
    a.set_title('BJ eps=%.1f (kappa*Hp=%.2f)'%(eps,sol['kappa']*HP)); a.set_xlabel('u (m/s)'); a.set_ylabel('y (m)')
    a.legend(fontsize=7); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(F+'/bj_profiles.png',dpi=140); print('bj_profiles.png saved')
