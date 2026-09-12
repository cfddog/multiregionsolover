#!/usr/bin/env python3
"""Parameterized 1D LTNE validator: closed-form two-temperature solution.
Usage: validate_ltne_param.py DIR EPS G Q HV [HC] [TINF]
"""
import sys, os
import numpy as np
DIR=sys.argv[1]; eps=float(sys.argv[2]); G=float(sys.argv[3]); q=float(sys.argv[4]); hv=float(sys.argv[5])
hc=float(sys.argv[6]) if len(sys.argv)>6 else 31.4
Tinf=float(sys.argv[7]) if len(sys.argv)>7 else 300.0
L=0.010; ks=13.4; cp=1007.0; Tin=300.0; kse=(1.0-eps)*ks
def analytic(y):
    A=hv/(G*cp); B=hv/kse; beta=kse/(G*cp); D=np.sqrt(A*A+4*B)
    l1=(-A+D)/2; l2=(-A-D)/2
    e1,e2=np.exp(l1*y),np.exp(l2*y); e1L,e2L=np.exp(l1*L),np.exp(l2*L)
    M=np.array([[1.0,beta*l1,beta*l2],[-hc,kse*l1-hc,kse*l2-hc],[0.0,kse*l1*e1L,kse*l2*e2L]])
    C0,C1,C2=np.linalg.solve(M,np.array([Tin,-hc*Tinf,q]))
    Ts=C0+C1*e1+C2*e2; Tf=C0+beta*(C1*l1*e1+C2*l2*e2)
    return Ts,Tf
def read_Ts(p):
    d=np.loadtxt(p,skiprows=3); return d[:,0],d[:,3]
def read_flow_t(p):
    lines=open(p).read().splitlines(); vals=[]; i=0
    while i<len(lines):
        s=lines[i].strip()
        if s.startswith("SCALARS") and lines[i].split()[1]=="temperature":
            i+=2
            while i<len(lines) and not lines[i].strip().startswith(("SCALARS","VECTORS","POINT_DATA","CELL_DATA","FIELD")):
                vals+=lines[i].split(); i+=1
            break
        i+=1
    return np.array(vals,float)
x,Ts=read_Ts(os.path.join(DIR,'Ts_block_1.dat'))
Tf=read_flow_t(os.path.join(DIR,'flow3d_block_1.vtk'))[:x.size]
Ts_a,Tf_a=analytic(x)
eTf=np.sqrt(np.mean((Tf-Tf_a)**2)); eTs=np.sqrt(np.mean((Ts-Ts_a)**2))
mTf=np.max(np.abs(Tf-Tf_a)); mTs=np.max(np.abs(Ts-Ts_a))
i0=np.argmin(np.abs(x-0)); iL=np.argmin(np.abs(x-L))
out=G*cp*(Tf[iL]-Tin)+hc*(Ts[i0]-Tinf)
print('eps=%.2f G=%.1f q=%.3g hv=%.4e | Tf RMS=%.4f max=%.4f | Ts RMS=%.4f max=%.4f | energy gap=%.3f%%'
      %(eps,G,q,hv,eTf,mTf,eTs,mTs,100*(q-out)/q))
print('  Tf(L)=%.2f (an %.2f)  Ts(0)=%.2f (an %.2f)  Ts(L)-Tf(L)=%.2f'
      %(Tf[iL],Tf_a[iL],Ts[i0],Ts_a[i0],Ts[iL]-Tf[iL]))
