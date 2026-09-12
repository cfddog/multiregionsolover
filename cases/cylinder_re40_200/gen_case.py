#!/usr/bin/env python3
"""Cylinder polar O-grid (full 360 deg) for ReD=40/200; AC solver."""
import struct, array, math
D=1.0; R_in=D/2; R_out=25.0*D
nx, ny, nz = 121, 121, 2; beta=3.0
xs=[];ys=[];zs=[]
def r_of(i):
    xi=i/(nx-1); return R_in+(R_out-R_in)*(1.0-math.tanh(beta*(1.0-xi))/math.tanh(beta))
for k in range(nz):
    for j in range(ny):
        th=2*math.pi*j/ny
        for i in range(nx):
            r=r_of(i); xs.append(r*math.cos(th)); ys.append(r*math.sin(th)); zs.append(float(k))
def wr(f,data,fmt):
    raw=data.tobytes() if isinstance(data,array.array) else struct.pack(fmt,*data)
    f.write(struct.pack('I',len(raw))); f.write(raw); f.write(struct.pack('I',len(raw)))
with open('Mesh3d.x','wb') as f:
    wr(f,[1],'i'); wr(f,[nx,ny,nz],'iii')
    wr(f,array.array('d',xs),''); wr(f,array.array('d',ys),''); wr(f,array.array('d',zs),'')
ny4=ny//4
faces=[(1,nx,1,ny,1,1,3),(1,nx,1,ny,nz,nz,3),(1,1,1,ny,1,nz,2),
       (nx,nx,1,ny4+1,1,nz,22),(nx,nx,ny4+1,3*ny4+1,1,nz,21),(nx,nx,3*ny4+1,ny,1,nz,22)]
lines=['! cylinder polar O-grid','     1','%6d%6d%6d'%(nx,ny,nz),'cyl','%d'%8]
for f in faces: lines.append(''.join('%8d'%t for t in f))
lines.append(''.join('%8d'%t for t in (1,nx,1,1,1,nz,-1))); lines.append(''.join('%8d'%t for t in (1,nx,ny,ny,1,nz,1)))
lines.append(''.join('%8d'%t for t in (1,nx,ny,ny,1,nz,-1))); lines.append(''.join('%8d'%t for t in (1,nx,1,1,1,nz,1)))
open('bc3d.inp','w').write('\n'.join(lines)+'\n')
open('material.in','w').write('1\n2\n0.0 0.0 0.0\n')
print('mesh %dx%dx%d R_in=%.2f R_out=%.1f'%(nx,ny,nz,R_in,R_out))
