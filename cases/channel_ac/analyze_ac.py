import numpy as np
lines=open('flow3d_block_1.vtk').read().splitlines()
def find(pre):
  for n,ln in enumerate(lines):
    if ln.startswith(pre): return n
  raise SystemExit('nf')
p0=find('POINTS'); pts=int(lines[p0].split()[1])
XYZ=np.array([[float(x) for x in l.split()] for l in lines[p0+1:p0+1+pts]])
cd=find('CELL_DATA'); cl=int(lines[cd].split()[1])
marks=[(n,ln.split()[0],ln.split()[1]) for n,ln in enumerate(lines) if ln.startswith(('SCALARS','VECTORS'))]
D={}
for ii in range(len(marks)):
  n0,typ,nm=marks[ii]; n1=marks[ii+1][0] if ii+1<len(marks) else len(lines)
  body=[float(x) for l in lines[n0+1:n1] for x in l.split() if not l.startswith('LOOKUP')]
  if typ=='SCALARS': D[nm]=np.array(body[:cl])
  else: D[nm]=np.array(body[:3*cl]).reshape(-1,3)
XYZr=XYZ.reshape(2,41,161,3)
Yc=0.5*(XYZr[0,1:,:,1]+XYZr[0,:-1,:,1])[:,0]   # 40
Xc=0.5*(XYZr[0,0,:-1,0]+XYZr[0,0,1:,0])          # 160
u=D['velocity'][:,0].reshape(40,160)   # [j,i]
pcell=D['pressure'].reshape(40,160)
upv=6.0*Yc*(1.0-Yc)
print(' i   xc     umax     mean    RMSvsPois')
for i in [20,40,60,80,100,120,140,158]:
  prof=u[:,i]
  print('%3d %6.3f  %7.4f  %7.4f   %8.5f'%(i,Xc[i],prof.max(),prof.mean(),np.sqrt(np.mean((prof-upv)**2))))
for i in (120,150,158):
  prof=u[:,i]
  print('\nstation i=%d x=%.3f  u_max=%.4f mean=%.4f'%(i,Xc[i],prof.max(),prof.mean()))
  print(' y_j    u_AC    u_Pois')
  for a,b,c in zip(Yc,prof,upv): print('%7.4f %7.4f %7.4f'%(a,b,c))
print('\ndp/dx between x=4.0 and 7.9:', (pcell[:,150].mean()-pcell[:,80].mean())/(Xc[150]-Xc[80]))
