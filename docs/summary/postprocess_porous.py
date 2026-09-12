#!/usr/bin/env python3
"""Consolidated post-processing: Porous plug (Darcy-Brinkman) and Beavers-Joseph."""
import os, sys, glob
import numpy as np
repo='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
sys.path.insert(0, os.path.join(repo,'cases/porous_plug/plot_validation'))
sys.path.insert(0, os.path.join(repo,'cases/beavers_joseph/plot_validation'))
import pv_io, pv_vtk, bj_analytical as bj

def do_plug():
    print('===== Porous plug (pressure-driven, AC) vs Darcy-Brinkman =====')
    for eps in (0.3,0.6,0.9):
        d=os.path.join(repo,'cases/porous_plug_ac/eps%02d'%int(round(eps*10)))
        p,T,uvw=pv_io.load_vtk(os.path.join(d,'flow3d_block_1.vtk'))
        u=uvw[:,:,0]; xs,ys=pv_io.cell_centers()
        imid=np.arange(10,34)
        coef=np.polyfit(xs[imid], p[10, imid], 1); dpdx=coef[0]
        ubar=u.mean(axis=0)[imid].mean()
        ub,uD,K=pv_io.brinkman_profile(ys, dpdx, eps, 0.10, 0.02)
        lam=np.sqrt(K)/eps; alpha=0.5*pv_io.LY/lam
        ubarB=uD*(1.0-np.tanh(alpha)/alpha)
        print('  eps=%.2f  K=%.4e m2  lam=%.4e m  fitted dp/dx=%.5e Pa/m'%(eps,K,lam,dpdx))
        print('      ubar_CFD=%.5e  ubar_Brinkman=%.5e  err=%+.2f%%   u_max=%.5e'%(ubar,ubarB,100*(ubar-ubarB)/ubarB,u[:,10:34].max()))

def do_bj():
    print('===== Beavers-Joseph (two-block low-speed/porous, AC) vs composite Darcy-Brinkman =====')
    HP=1.2; H=1.0; DP=0.5; MU=0.02
    for eps in (0.3,0.6,0.9):
        d=os.path.join(repo,'cases/beavers_joseph/eps%02d'%int(round(eps*10)))
        b1=pv_vtk.load_block(os.path.join(d,'flow3d_block_1.vtk'))
        b2=pv_vtk.load_block(os.path.join(d,'flow3d_block_2.vtk'))
        x=b1['xc'][0]; dx=x[1]-x[0]; L=x[-1]-x[0]+dx
        xwin=(x>0.25*L)&(x<0.75*L); jf=b1['p'].shape[0]//2; jp=b2['p'].shape[0]//2
        cf=np.polyfit(x[xwin],b1['p'][jf,xwin],1); cp=np.polyfit(x[xwin],b2['p'][jp,xwin],1)
        G=-0.5*(cf[0]+cp[0])
        sol=bj.solve_composite(G,eps,DP,MU,h=H,Hp=HP)
        i0=int(np.argmin(np.abs(b1['xc'][0]-0.5*L)))
        uf=b1['u'][:,i0]; up=b2['u'][:,i0]; yf=b1['yc'][:,i0]; yp=b2['yc'][:,i0]
        yc=np.concatenate([yp,yf]); uc=np.concatenate([up,uf]); ua=bj.profile(yc,sol)
        umax=np.max(np.abs(bj.profile(np.linspace(-HP,H,2001),sol)))
        L2=np.sqrt(np.mean(((uc-ua)/umax)**2)); Linf=np.max(np.abs((uc-ua)/umax))
        print('  eps=%.2f  G=-dp/dx=%.5e Pa/m  kappa*Hp=%.2f  uD=%.4e u_int=%.4e'%(eps,G,sol['kappa']*HP,sol['uD'],sol['u_int']))
        print('      profile L2=%.3e  Linf=%.3e   (umax=%.4e)'%(L2,Linf,umax))

if __name__=='__main__':
    do_plug(); do_bj()
