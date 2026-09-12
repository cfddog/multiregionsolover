#!/usr/bin/env python3
import sys, math
import numpy as np
dir=sys.argv[1]; mu=float(sys.argv[2]); U=1.0; R_in=0.5; ny=121; nx=121; dz=1.0
lines=open(dir+'/flow3d_block_1.vtk').read().splitlines()
n0=[n for n,l in enumerate(lines) if l.startswith('POINTS')][0]; pts=int(lines[n0].split()[1])
A=[]; i=n0+1
while len(A)<pts: A.append([float(x) for x in lines[i].split()]); i+=1
X=np.array(A).reshape(2,nx,ny,3)
D={}; marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(lines) if ln.startswith(('SCALARS','VECTORS'))]
cl=int([l for l in lines if l.startswith('CELL_DATA')][0].split()[1])
for ii in range(len(marks)):
    a,typ,nm=marks[ii]; b=marks[ii+1][0] if ii+1<len(marks) else len(lines)
    body=[float(x) for l in lines[a+1:b] for x in l.split() if not l.startswith('LOOKUP')]
    D[nm]=np.array(body[:cl]) if typ=='SCALARS' else np.array(body[:3*cl]).reshape(-1,3)
nc=(nx-1)*(ny-1)
vel=D['velocity'][:nc].reshape(ny-1,nx-1,3)   # [j,i,comp] (i fastest)
pre=D['pressure'][:nc].reshape(ny-1,nx-1)
# first radial cell i=0 (adjacent to wall), wall radius R_in=0.5
r0=X[0,0,0,0]*0+0.5; r1=X[0,1,0,0]
dr=(r1-r0)/2.0
rc=(r0+r1)/2.0
r_cell=X[0,0,:,0]  # not used
# cell centres theta
dth=2*math.pi/(ny-1)
th=(np.arange(ny-1)+0.5)*dth
Ar=0.5*dth*dz   # cylinder surface segment area
Fx_p=0.0; Fx_v=0.0; Fy_p=0.0; Fy_v=0.0
for j in range(ny-1):
    u=vel[j,0,0]; v=vel[j,0,1]; p=pre[j,0]
    ct=math.cos(th[j]); st=math.sin(th[j])
    ut=-u*st+v*ct                 # tangential velocity
    tau=mu*ut/dr                  # wall shear (tangent dir t=(-st,ct))
    # force on body = -p n dA (n outward from body = (ct,st)) + tau*t dA
    Fx_p += -p*ct*Ar; Fy_p += -p*st*Ar
    Fx_v += tau*(-st)*Ar; Fy_v += tau*(ct)*Ar
Fx=Fx_p+Fx_v
Cd=2*Fx/(1.0*U*U*1.0)
Cd_p=2*Fx_p/(1.0*U*U*1.0); Cd_v=2*Fx_v/(1.0*U*U*1.0)
# recirculation length along y=0 downstream: u=0 crossing
# find j nearest theta=0 (y=0,x>0): th=0 -> j=0 ; use x axis: cells along x at j=0
xline=[]; uline=[]
for i in range(nx-1):
    xc=0.5*(X[0,i,0,0]+X[0,i+1,0,0]); yc=0.0
    xline.append(xc); uline.append(vel[0,i,0])
xline=np.array(xline); uline=np.array(uline)
# first sign change from - to + after cylinder
Lr=None
for i in range(1,len(xline)):
    if uline[i-1]<0 and uline[i]>=0:
        Lr=xline[i]; break
print('Cd(total)=%.4f  (pressure %.4f + viscous %.4f)   Cl=%.4f'%(Cd,Cd_p,Cd_v,2*(Fy_p+Fy_v)/(U*U*1.0)))
print('recirculation length Lr=%.4f  (Lr/D=%.3f)'%(Lr if Lr else -1, (Lr if Lr else -1)/1.0))
print('u_min on wake centerline = %.5f'%uline.min())
