#!/usr/bin/env python3
import numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'; F=R+'/docs/summary/figures'
def T1d(path):
    d=np.loadtxt(path,skiprows=3); x,T=d[:,0],d[:,3]
    xu,inv=np.unique(np.round(x,10),return_inverse=True); Tx=np.array([T[inv==k].mean() for k in range(xu.size)])
    return xu,Tx
def ann(path):
    d=np.loadtxt(path,skiprows=3); r=np.hypot(d[:,0],d[:,1]); T=d[:,3]
    ii=np.arange(d.shape[0])%40; ri=np.array([r[ii==k].mean() for k in range(40)]); Ti=np.array([T[ii==k].mean() for k in range(40)])
    return ri,Ti
fig,ax=plt.subplots(2,2,figsize=(9.6,7.6))
# (a) 1D Dirichlet-Dirichlet
x,T=T1d('/home/sundong/Fortran_Project/OpenCFD-EC-1.16a/cases/solid_1d_2temp/Ts_block_1.dat')
Tref='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a/cases/solid_1d'
x2,T2=T1d('/tmp/run_solid_1d/Ts_block_1.dat')
ax[0,0].plot(x,T,'o',ms=3,label='CFD'); ax[0,0].plot(x,300+100*x,'k-',lw=1.3,label='analytic T=300+100x'); ax[0,0].set_title('(a) 1D: both ends fixed T (300/400 K)'); ax[0,0].set_xlabel('x (m)'); ax[0,0].set_ylabel('T (K)')
ax[0,0].legend(fontsize=8); ax[0,0].grid(alpha=.3)
print('1D Dirichlet: max err = %.3e K'%np.max(np.abs(T-(300+100*x))))
# (b) 1D Dirichlet-Neumann (flux)
ax[0,1].plot(x2,T2,'o',ms=3,label='CFD'); ax[0,1].plot(x2,300+(800/15.1)*x2,'k-',lw=1.3,label='analytic T=300+(q/k)x'); ax[0,1].set_title('(b) 1D: T(0)=300, q(L)=800 W/m2'); ax[0,1].set_xlabel('x (m)'); ax[0,1].set_ylabel('T (K)')
ax[0,1].legend(fontsize=8); ax[0,1].grid(alpha=.3)
print('1D Neumann: max err = %.3e K'%np.max(np.abs(T2-(300+(800/15.1)*x2))))
# (c) annulus Dirichlet-Dirichlet
r,Ta=ann('/home/sundong/Fortran_Project/OpenCFD-EC-1.16a/cases/solid_annulus_2temp/Ts_block_1.dat')
r2,Ta2=ann('/tmp/run_solid_annulus/Ts_block_1.dat')
an_=350-50*np.log(r/0.1)/np.log(2.0)
ax[1,0].plot(r,Ta,'o',ms=3,label='CFD'); ax[1,0].plot(r,an_,'k-',lw=1.3,label='analytic log'); ax[1,0].set_title('(c) annulus: inner/outer fixed T (350/300 K)'); ax[1,0].set_xlabel('r (m)'); ax[1,0].set_ylabel('T (K)')
ax[1,0].legend(fontsize=8); ax[1,0].grid(alpha=.3)
print('annulus Dirichlet: max err = %.3e K'%np.max(np.abs(Ta-an_)))
an2=300+(800*0.1/15.1)*np.log(0.2/r2)
ax[1,1].plot(r2,Ta2,'o',ms=3,label='CFD'); ax[1,1].plot(r2,an2,'k-',lw=1.3,label='analytic T=To+(q Ri/k)ln(Ro/r)'); ax[1,1].set_title('(d) annulus: inner q=800, outer T=300'); ax[1,1].set_xlabel('r (m)'); ax[1,1].set_ylabel('T (K)')
ax[1,1].legend(fontsize=8); ax[1,1].grid(alpha=.3)
print('annulus Neumann: max err = %.3e K'%np.max(np.abs(Ta2-an2)))
plt.tight_layout(); plt.savefig(F+'/conduction_bc.png',dpi=140); print('conduction_bc.png saved')
