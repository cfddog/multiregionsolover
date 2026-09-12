#!/usr/bin/env python3
"""Cylinder polar O-grid as TWO patched blocks (theta 0..pi and pi..2pi),
so the seam is a normal block-block interface (no self-periodic). AC solver."""
import struct, array, math
D=1.0; R_in=D/2; R_out=25.0*D
nx, ny, nz = 121, 61, 2; beta=3.0
def r_of(i):
    xi=i/(nx-1); return R_in+(R_out-R_in)*(1.0-math.tanh(beta*(1.0-xi))/math.tanh(beta))
def wr(f,data,fmt):
    raw=data.tobytes() if isinstance(data,array.array) else struct.pack(fmt,*data)
    f.write(struct.pack('I',len(raw))); f.write(raw); f.write(struct.pack('I',len(raw)))
xs=[];ys=[];zs=[]
for m in range(2):
    for k in range(nz):
        for j in range(ny):
            th=m*math.pi+math.pi*j/(ny-1)
            for i in range(nx):
                r=r_of(i); xs.append(r*math.cos(th)); ys.append(r*math.sin(th)); zs.append(float(k))
with open('Mesh3d.x','wb') as f:
    wr(f,[2],'i'); wr(f,[nx,ny,nz,nx,ny,nz],'6i')
    # block1 then block2, each 3 coord records (separate layout)
    n=nx*ny*nz
    wr(f,array.array('d',xs[:n]),''); wr(f,array.array('d',ys[:n]),''); wr(f,array.array('d',zs[:n]),'')
    wr(f,array.array('d',xs[n:]),''); wr(f,array.array('d',ys[n:]),''); wr(f,array.array('d',zs[n:]),'')
def F(*t): return ''.join('%8d'%x for x in t)
h=ny//2
lines=['! cylinder two-block polar O-grid','     2',
       '%6d%6d%6d'%(nx,ny,nz),'blk-1','7',
       F(1,nx,1,ny,1,1,3), F(1,nx,1,ny,nz,nz,3), F(1,1,1,ny,1,nz,2),
       F(nx,nx,1,h+1,1,nz,22), F(nx,nx,h+1,ny,1,nz,21),
       F(1,nx,1,1,1,nz,-1), F(1,nx,ny,ny,1,nz,2),
       F(1,nx,ny,ny,1,nz,-1), F(1,nx,1,1,1,nz,2),
       '%6d%6d%6d'%(nx,ny,nz),'blk-2','7',
       F(1,nx,1,ny,1,1,3), F(1,nx,1,ny,nz,nz,3), F(1,1,1,ny,1,nz,2),
       F(nx,nx,1,h+1,1,nz,21), F(nx,nx,h+1,ny,1,nz,22),
       F(1,nx,1,1,1,nz,-1), F(1,nx,ny,ny,1,nz,1),
       F(1,nx,ny,ny,1,nz,-1), F(1,nx,1,1,1,nz,1)]
open('bc3d.inp','w').write('\n'.join(lines)+'\n')
open('material.in','w').write('2\n2 2\n0.0 0.0 0.0\n0.0 0.0 0.0\n')
print('2-block cylinder %dx%dx%d x2, R_in=%.2f R_out=%.1f'%(nx,ny,nz,R_in,R_out))
