#!/usr/bin/env python3
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'; F=R+'/docs/summary/figures'
def T1d(p):
    d=np.loadtxt(p,skiprows=3); x,T=d[:,0],d[:,3]
    xu,inv=np.unique(np.round(x,10),return_inverse=True); Tx=np.array([T[inv==k].mean() for k in range(xu.size)])
    return xu,Tx
def ann(p):
    d=np.loadtxt(p,skiprows=3); r=np.hypot(d[:,0],d[:,1]); T=d[:,3]
    ii=np.arange(d.shape[0])%40; ri=np.array([r[ii==k].mean() for k in range(40)]); Ti=np.array([T[ii==k].mean() for k in range(40)])
    return ri,Ti
xDD,TDD=T1d(R+'/cases/solid_1d_2temp/Ts_block_1.dat')
xDN,TDN=T1d('/tmp/run_solid_1d/Ts_block_1.dat')
rDD,TaDD=ann(R+'/cases/solid_annulus_2temp/Ts_block_1.dat')
rDN,TaDN=ann('/tmp/run_solid_annulus/Ts_block_1.dat')
fig,ax=plt.subplots(2,2,figsize=(9.8,7.8))
k=15.1
# (a) 1D T(x)
ax[0,0].plot(xDD,TDD,'o',ms=3,mfc='none',color='tab:red',label='CFD: T=300/400 K')
ax[0,0].plot(xDD,300+100*xDD,'r-',lw=1.2,label='analytic (D-D)')
ax[0,0].plot(xDN,TDN,'s',ms=3,mfc='none',color='tab:blue',label='CFD: T(0)=300, q(L)=800')
ax[0,0].plot(xDN,300+(800/k)*xDN,'b-',lw=1.2,label='analytic (D-N)')
ax[0,0].set_title('(a) 1D: T(x) for two BC types'); ax[0,0].set_xlabel('x (m)'); ax[0,0].set_ylabel('T (K)')
ax[0,0].legend(fontsize=7); ax[0,0].grid(alpha=.3)
# (b) 1D q(x)
ax[0,1].plot(xDD[1:-1],k*np.gradient(TDD,xDD)[1:-1],'o',ms=3,mfc='none',color='tab:red',label='q: D-D (const 1510)')
ax[0,1].axhline(k*100,color='r',ls='--',lw=1.1)
ax[0,1].plot(xDN[1:-1],k*np.gradient(TDN,xDN)[1:-1],'s',ms=3,mfc='none',color='tab:blue',label='q: D-N (const 800)')
ax[0,1].axhline(800,color='b',ls='--',lw=1.1)
ax[0,1].set_title('(b) 1D: heat flux q=k dT/dx (constant)'); ax[0,1].set_xlabel('x (m)'); ax[0,1].set_ylabel('q (W/m2)')
ax[0,1].legend(fontsize=7); ax[0,1].grid(alpha=.3)
# (c) annulus T(r)
ax[1,0].plot(rDD,TaDD,'o',ms=3,mfc='none',color='tab:red',label='CFD: 350/300 K')
ax[1,0].plot(rDD,350-50*np.log(rDD/0.1)/np.log(2.0),'r-',lw=1.2,label='analytic (D-D)')
ax[1,0].plot(rDN,TaDN,'s',ms=3,mfc='none',color='tab:blue',label='CFD: q_i=800, T_o=300')
ax[1,0].plot(rDN,300+(800*0.1/k)*np.log(0.2/rDN),'b-',lw=1.2,label='analytic (D-N)')
ax[1,0].set_title('(c) annulus: T(r) for two BC types'); ax[1,0].set_xlabel('r (m)'); ax[1,0].set_ylabel('T (K)')
ax[1,0].legend(fontsize=7); ax[1,0].grid(alpha=.3)
# (d) annulus q_r(r)
ax[1,1].plot(rDD,-k*np.gradient(TaDD,rDD),'o',ms=3,mfc='none',color='tab:red',label='q_r: D-D')
ax[1,1].plot(rDD,50*k/np.log(2.0)/rDD,'r-',lw=1.2)
ax[1,1].plot(rDN,-k*np.gradient(TaDN,rDN),'s',ms=3,mfc='none',color='tab:blue',label='q_r: D-N')
ax[1,1].plot(rDN,800*0.1/rDN,'b-',lw=1.2)
ax[1,1].set_title('(d) annulus: radial heat flux q_r (CFD vs analytic)'); ax[1,1].set_xlabel('r (m)'); ax[1,1].set_ylabel('q_r (W/m2)')
ax[1,1].legend(fontsize=7); ax[1,1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(F+'/conduction_bc.png',dpi=140)
print('1D D-D err %.2e'%np.max(np.abs(TDD-(300+100*xDD))))
print('1D D-N err %.2e'%np.max(np.abs(TDN-(300+(800/k)*xDN))))
print('ann D-D err %.2e'%np.max(np.abs(TaDD-(350-50*np.log(rDD/0.1)/np.log(2)))))
print('ann D-N err %.2e'%np.max(np.abs(TaDN-(300+(800*0.1/k)*np.log(0.2/rDN)))))
print('saved conduction_bc.png')
