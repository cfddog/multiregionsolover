#!/usr/bin/env python3
import os,sys,json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
root='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
sys.path.insert(0,root+'/cases/porous_plug/plot_validation'); import pv_io
figs=root+'/docs/summary/figures'; os.makedirs(figs,exist_ok=True)

def read_flow_t(p):
    L=open(p).read().splitlines(); v=[]; i=0
    while i<len(L):
        if L[i].strip().startswith('SCALARS') and L[i].split()[1]=='temperature':
            i+=2
            while i<len(L) and not L[i].strip().startswith(('SCALARS','VECTORS','POINT_DATA','CELL_DATA','FIELD')): v+=L[i].split(); i+=1
            break
        i+=1
    return np.array(v,float)

# ---------- (3) heat conduction: add heat-flux comparison ----------
D=np.loadtxt('/tmp/run_solid_1d/Ts_block_1.dat',skiprows=3); x,T=D[:,0],D[:,3]
xu,inv=np.unique(np.round(x,10),return_inverse=True); Tx=np.array([T[inv==k].mean() for k in range(xu.size)])
q=15.1*np.gradient(Tx,xu)
fig,ax=plt.subplots(1,2,figsize=(9,3.8))
ax[0].plot(xu,Tx,'o',ms=3,label='FVM'); ax[0].plot(xu,300+(800/15.1)*xu,'k-',lw=1.3,label='analytic T=300+(q/k)x')
ax[0].set_xlabel('x (m)'); ax[0].set_ylabel('T (K)'); ax[0].set_title('(a) solid 1D: T(x)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].plot(xu,q,'o',ms=3,label='CFD q=k dT/dx'); ax[1].axhline(800,color='k',ls='--',lw=1.3,label='q=800 W/m2')
ax[1].set_xlabel('x (m)'); ax[1].set_ylabel('q (W/m2)'); ax[1].set_ylim(780,820); ax[1].set_title('(b) solid 1D: heat flux'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(figs+'/solid1d_flux.png',dpi=140); print('solid1d_flux.png')

D=np.loadtxt('/tmp/run_solid_annulus/Ts_block_1.dat',skiprows=3); X,Y,T=D[:,0],D[:,1],D[:,3]; r=np.hypot(X,Y)
ii_=np.arange(D.shape[0])%40; ri=np.array([r[ii_==k].mean() for k in range(40)]); Ti=np.array([T[ii_==k].mean() for k in range(40)])
qr=15.1*np.gradient(Ti,ri)
fig,ax=plt.subplots(1,2,figsize=(9,3.8))
ax[0].plot(ri,Ti,'o',ms=3,label='FVM'); ax[0].plot(ri,300+(800*0.1/15.1)*np.log(0.2/ri),'k-',lw=1.3,label='analytic')
ax[0].set_xlabel('r (m)'); ax[0].set_ylabel('T (K)'); ax[0].set_title('(a) annulus: T(r)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].plot(ri,-qr,'o',ms=3,label='CFD q_r'); ax[1].plot(ri,800*0.1/ri,'k-',lw=1.3,label='analytic q=q_i Ri/r')
ax[1].set_xlabel('r (m)'); ax[1].set_ylabel('q_r (W/m2)'); ax[1].set_title('(b) annulus: radial heat flux'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(figs+'/annulus_flux.png',dpi=140); print('annulus_flux.png')

# ---------- (4) LTNE figures (eps / G / q), sparse markers ----------
def ltne_curve(d,eps,G,q,hv):
    dd=np.loadtxt(d+'/Ts_block_1.dat',skiprows=3); xx,Ts=dd[:,0],dd[:,3]; Tf=read_flow_t(d+'/flow3d_block_1.vtk')[:xx.size]
    L=0.01; ks=13.4; cp=1007.0; hc=31.4; Tinf=300.; Tin=300.; kse=(1-eps)*ks
    A=hv/(G*cp); B=hv/kse; beta=kse/(G*cp); Ds=np.sqrt(A*A+4*B); l1=(-A+Ds)/2; l2=(-A-Ds)/2
    e1,e2=np.exp(l1*xx),np.exp(l2*xx); e1L,e2L=np.exp(l1*L),np.exp(l2*L)
    M=np.array([[1,beta*l1,beta*l2],[-hc,kse*l1-hc,kse*l2-hc],[0,kse*l1*e1L,kse*l2*e2L]])
    C0,C1,C2=np.linalg.solve(M,np.array([Tin,-hc*Tinf,q])); Tsa=C0+C1*e1+C2*e2; Tfa=C0+beta*(C1*l1*e1+C2*l2*e2)
    return xx*1e3,Tf,Ts,Tfa,Tsa
rows=json.load(open(root+'/cases/porous_ltne_1d/sweep/cases.json'))
hvmap={(round(e,2),round(G,1)):hv for tag,e,G,q,hv,Re,Nu in rows}
def panel(ax,d,eps,G,q,hv,title):
    xm,Tf,Ts,Tfa,Tsa=ltne_curve(d,eps,G,q,hv)
    ax.plot(xm[::6],Tf[::6],'s',ms=3,mfc='none',color='tab:blue',label='Tf CFD')
    ax.plot(xm,Tfa,'b-',lw=1.3,label='Tf analytic')
    ax.plot(xm[::6],Ts[::6],'o',ms=3,mfc='none',color='tab:red',label='Ts CFD')
    ax.plot(xm,Tsa,'r-',lw=1.3,label='Ts analytic')
    ax.set_xlabel('y (mm)'); ax.set_ylabel('T (K)'); ax.set_title(title); ax.legend(fontsize=7); ax.grid(alpha=.3)
S=root+'/cases/porous_ltne_1d/sweep'
fig,ax=plt.subplots(1,3,figsize=(12.5,3.8))
for a,eps in zip(ax,(0.3,0.6,0.9)):
    panel(a,'%s/eps%02d_G2_q2'%(S,int(round(eps*10))),eps,2.0,2e6,hvmap[(round(eps,2),2.0)],'eps=%.1f (G=2, q=2e6)'%eps)
plt.tight_layout(); plt.savefig(figs+'/ltne_eps.png',dpi=140); print('ltne_eps.png')
fig,ax=plt.subplots(1,3,figsize=(12.5,3.8))
for a,G in zip(ax,(1.0,2.0,3.0)):
    panel(a,'%s/eps03_G%d_q2'%(S,int(G)),0.3,G,2e6,hvmap[(0.3,float(int(G)))],'G=%d (eps=0.3, q=2e6)'%int(G))
plt.tight_layout(); plt.savefig(figs+'/ltne_G.png',dpi=140); print('ltne_G.png')
fig,ax=plt.subplots(1,3,figsize=(12.5,3.8))
for a,qq in zip(ax,(1,2,3)):
    panel(a,'%s/eps03_G2_q%d'%(S,qq),0.3,2.0,qq*1e6,hvmap[(0.3,2.0)],'q=%de6 (eps=0.3, G=2)'%qq)
plt.tight_layout(); plt.savefig(figs+'/ltne_q.png',dpi=140); print('ltne_q.png')

# ---------- (5) plug: separate panel per porosity ----------
fig,ax=plt.subplots(1,3,figsize=(12.5,3.8))
for a,eps in zip(ax,(0.3,0.6,0.9)):
    d=root+'/cases/porous_plug_ac/eps%02d'%int(round(eps*10))
    p,T,uvw=pv_io.load_vtk(d+'/flow3d_block_1.vtk'); uu=uvw[:,:,0]; xs,ys=pv_io.cell_centers()
    coef=np.polyfit(xs[10:34],p[10,10:34],1); ub,uD,K=pv_io.brinkman_profile(ys,coef[0],eps,0.10,0.02)
    a.plot(uu[:,20],ys,'o',ms=3,mfc='none',label='CFD'); a.plot(ub,ys,'k-',lw=1.3,label='Brinkman')
    a.set_xlabel('u (m/s)'); a.set_ylabel('y (m)'); a.set_title('eps=%.1f (K=%.2e)'%(eps,K)); a.legend(fontsize=8); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig(figs+'/plug_eps.png',dpi=140); print('plug_eps.png')
print('done')
