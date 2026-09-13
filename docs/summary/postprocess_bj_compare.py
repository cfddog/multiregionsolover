#!/usr/bin/env python3
import os,sys
import numpy as np
R='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
sys.path.insert(0,R+'/cases/beavers_joseph/plot_validation')
import pv_vtk, bj_analytical as bj
HP=1.2; H=1.0; DP=0.5; MU=0.02
def L2(d,eps):
    b1=pv_vtk.load_block(d+'/flow3d_block_1.vtk'); b2=pv_vtk.load_block(d+'/flow3d_block_2.vtk')
    x=b1['xc'][0]; dx=x[1]-x[0]; L=x[-1]-x[0]+dx
    xw=(x>0.25*L)&(x<0.75*L); jf=b1['p'].shape[0]//2; jp=b2['p'].shape[0]//2
    cf=np.polyfit(x[xw],b1['p'][jf,xw],1); cp=np.polyfit(x[xw],b2['p'][jp,xw],1); G=-0.5*(cf[0]+cp[0])
    sol=bj.solve_composite(G,eps,DP,MU,h=H,Hp=HP)
    i0=int(np.argmin(np.abs(b1['xc'][0]-0.5*L)))
    uf=b1['u'][:,i0]; up=b2['u'][:,i0]; yf=b1['yc'][:,i0]; yp=b2['yc'][:,i0]
    yc=np.concatenate([yp,yf]); uc=np.concatenate([up,uf]); ua=bj.profile(yc,sol)
    umax=np.max(np.abs(bj.profile(np.linspace(-HP,H,2001),sol)))
    return np.sqrt(np.mean(((uc-ua)/umax)**2)), np.max(np.abs((uc-ua)/umax)), G, sol['kappa']*HP
print('%-10s %-6s %-10s %-10s %-10s %-8s'%('mesh','eps','L2','Linf','G(-dp/dx)','k.Hp'))
for name,base in [('coarse',R+'/cases/beavers_joseph/eps%02d'),('refined',R+'/cases/beavers_joseph/refine/eps%02d')]:
    for eps in (0.3,0.6,0.9):
        d=base%int(round(eps*10))
        if not os.path.exists(d+'/flow3d_block_2.vtk'): print(name,eps,'no vtk'); continue
        l2,li,G,k= L2(d,eps); print('%-10s %-6.1f %-10.3e %-10.3e %-10.4e %-8.2f'%(name,eps,l2,li,G,k))
