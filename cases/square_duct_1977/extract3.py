import struct, os
d = '/home/sundong/Fortran_Project/OpenCFD-EC-1.16a/cases/square_duct_1977'
f = open(os.path.join(d, 'Mesh3d.x'), 'rb')
def rec():
    n = struct.unpack('I', f.read(4))[0]
    data = f.read(n); f.read(4)
    return data
nb = struct.unpack('i', rec())[0]
r = rec()
dims = struct.unpack('%di' % (len(r)//4), r)
B = []
for b in range(nb):
    ni, nj, nk = dims[3*b:3*b+3]
    a = struct.unpack('<%dd' % (3*ni*nj*nk), rec())
    B.append((ni, nj, nk, a))
sel = B[12:15]
for (ni, nj, nk, a) in sel:
    npts = ni*nj*nk
    xs, ys, zs = a[0:npts], a[npts:2*npts], a[2*npts:3*npts]
    # coordinate range along each index direction (edges of the block)
    def rng(vals, step, cnt):
        idx = lambda i: [vals[i + t*step] for t in range(cnt)]
        return min(idx(0)), max(idx(0)), min(idx(cnt-1)), max(idx(cnt-1))
    print('blk %dx%dx%d: x[%.4f,%.4f] y[%.4f,%.4f] z[%.4f,%.4f]'
          % (ni,nj,nk,min(xs),max(xs),min(ys),max(ys),min(zs),max(zs)))
    print('   i=1  corner xyz=(%.4f,%.4f,%.4f)   i=n xyz=(%.4f,%.4f,%.4f)'
          % (xs[0],ys[0],zs[0], xs[ni-1],ys[ni-1],zs[ni-1]))
    print('   j=1  xyz=(%.4f,%.4f,%.4f)   j=n xyz=(%.4f,%.4f,%.4f)'
          % (xs[0],ys[0],zs[0], xs[(nj-1)*ni],ys[(nj-1)*ni],zs[(nj-1)*ni]))
with open(os.path.join(d, 'Mesh3d.dat'), 'w') as g:
    g.write('%d\n' % len(sel))
    g.write(' '.join('%d %d %d' % (b[0],b[1],b[2]) for b in sel) + '\n')
    for (ni, nj, nk, a) in sel:
        npts = ni*nj*nk
        for c in range(3):
            g.write(' '.join('%.12e' % v for v in a[c*npts:(c+1)*npts]) + '\n')
print('wrote ASCII Mesh3d.dat (3 blocks)')
