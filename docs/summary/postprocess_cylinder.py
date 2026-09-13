#!/usr/bin/env python3
import sys, math
import numpy as np
dir=sys.argv[1]; mu=float(sys.argv[2]); U=1.0; R_in=0.5; nx=121; ny=61; dz=1.0
def read(fn):
    L=open(fn).read().splitlines(); n0=[n for n,l in enumerate(L) if l.startswith('POINTS')][0]
    npt=int(L[n0].split()[1]); A=[]; i=n0+1
    while len(A)<npt: A.append([float(x) for x in L[i].split()]); i+=1
    X=np.array(A).reshape(2,nx,ny,3)
    D={}; marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(L) if ln.startswith(('SCALARS','VECTORS'))]
    cl=int([l for l in L if l.startswith('CELL_DATA')][0].split()[1])
    for ii in range(len(marks)):
        a,t,nm=marks[ii]; b=marks[ii+1][0] if ii+1<len(marks) else len(L)
        body=[float(x) for l in L[a+1:b] for x in l.split() if not l.lookup() if False] if False else [float(x) for l in L[a+1:b] for x in l.split() if not l.startswith('LOOKUP')]
        D[nm]=np.array(body[:cl]) if t=='SCALARS' else np.array(body[:3*cl]).reshape(-1,3)
    return X,D
r0=0.5
Fx=0.0; Fy=0.0; Fxp=0.0
th0=0.0
for m,(fn,off) in enumerate([(dir+'/flow3d_block_1.vtk',0.0),(dir+'/flow3d_block_2.vtk',math.pi)]):
    X,D=read(fn); vel=D['velocity'][:(nx-1)*(ny-1)].reshape(ny-1,nx-1,3); pre=D['pressure'][:(nx-1)*(ny-1)].reshape(ny-1,nx-1)
    r1=X[0,1,0,0]; dr=(r1-r0)/2.0
    dth=math.pi/(ny-1); Ar=R_in*dth*dz
    for j in range(ny-1):
        th=off+(j+0.5)*dth; u=vel[j,0,0]; v=vel[j,0,1]; p=pre[j,0]
        ct=math.cos(th); st=math.sin(th); ut=-u*st+v*ct; tau=mu*ut/dr
        Fxp+=-p*ct*Ar; Fx+=(-p*ct+tau*(-st))*Ar; Fy+=(-p*st+tau*ct)*Ar
Cd=2*Fx/(U*U*1.0); Cdp=2*Fxp/(U*U*1.0)
# wake centerline: block1 j=0 (theta~0, x>0)
X,D=read(dir+'/flow3d_block_1.vtk'); vel=D['velocity'][:(nx-1)*(ny-1)].reshape(ny-1,nx-1,3)
xl=[]; ul=[]
for i in range(nx-1):
    xl.append(0.5*(X[0,i,0,0]+X[0,i+1,0,0])); ul.append(vel[0,i,0])
xl=np.array(xl); ul=np.array(ul); Lr=None
for i in range(1,len(xl)):
    if ul[i-1]<0 and ul[i]>=0: Lr=xl[i]; break
print('Cd=%.4f (pressure %.4f)  Cl=%.4f  Lr=%s  u_min(wake)=%.5f'%(Cd,Cdp,2*Fy/(U*U), ('%.3f'%Lr) if Lr else 'none', ul.min()))
