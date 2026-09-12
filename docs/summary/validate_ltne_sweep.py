#!/usr/bin/env python3
import os, json
import numpy as np
root='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
sweep=os.path.join(root,'cases/porous_ltne_1d/sweep')
rows=json.load(open(os.path.join(sweep,'cases.json')))
L=0.010; ks=13.4; cp=1007.0; Tin=300.0; hc=31.4; Tinf=300.0
def analytic(y,eps,G,q,hv):
    kse=(1-eps)*ks; A=hv/(G*cp); B=hv/kse; beta=kse/(G*cp); D=np.sqrt(A*A+4*B)
    l1=(-A+D)/2; l2=(-A-D)/2; e1,e2=np.exp(l1*y),np.exp(l2*y); e1L,e2L=np.exp(l1*L),np.exp(l2*L)
    M=np.array([[1,beta*l1,beta*l2],[-hc,kse*l1-hc,kse*l2-hc],[0,kse*l1*e1L,kse*l2*e2L]])
    C0,C1,C2=np.linalg.solve(M,np.array([Tin,-hc*Tinf,q]))
    return C0+C1*e1+C2*e2, C0+beta*(C1*l1*e1+C2*l2*e2)
def read_Ts(p):
    d=np.loadtxt(p,skiprows=3); return d[:,0],d[:,3]
def read_flow_t(p):
    lines=open(p).read().splitlines(); vals=[]; i=0
    while i<len(lines):
        s=lines[i].strip()
        if s.startswith('SCALARS') and lines[i].split()[1]=='temperature':
            i+=2
            while i<len(lines) and not lines[i].strip().startswith(('SCALARS','VECTORS','POINT_DATA','CELL_DATA','FIELD')):
                vals+=lines[i].split(); i+=1
            break
        i+=1
    return np.array(vals,float)
print('%-16s %5s %4s %5s %10s | %8s %8s %8s'%('case','eps','G','q(M)','hv','TfRMS','TsRMS','dq%'))
for tag,eps,G,q,hv,Re,Nu in rows:
    d=os.path.join(sweep,tag)
    x,Ts=read_Ts(os.path.join(d,'Ts_block_1.dat'))
    Tf=read_flow_t(os.path.join(d,'flow3d_block_1.vtk'))[:x.size]
    Ts_a,Tf_a=analytic(x,eps,G,q,hv)
    eTf=np.sqrt(np.mean((Tf-Tf_a)**2)); eTs=np.sqrt(np.mean((Ts-Ts_a)**2))
    i0=np.argmin(np.abs(x-0)); iL=np.argmin(np.abs(x-L))
    out=G*cp*(Tf[iL]-Tin)+hc*(Ts[i0]-Tinf)
    print('%-16s %5.2f %4.1f %5.1f %10.3e | %8.4f %8.4f %+7.3f'%(tag,eps,G,q/1e6,hv,eTf,eTs,100*(q-out)/q))
