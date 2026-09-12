#!/usr/bin/env python3
"""Generate the 1D LTNE transpiration-cooling sweep: 3 eps x 3 G x 3 q = 27 cases.
Template: cases/porous_ltne_1d/run_ac (AC solver, nx=401)."""
import os, shutil, math
repo='/home/sundong/Fortran_Project/OpenCFD-EC-1.16a'
Tmpl=os.path.join(repo,'cases/porous_ltne_1d/run_ac')
outroot=os.path.join(repo,'cases/porous_ltne_1d/sweep')
dp=1.0e-3; mu=1.846e-5; Pr=0.707; kf=0.0262
hc=31.4; Tinf=300.0
def hv_wakao(eps,G):
    Re=G*dp/mu; Nu=2.0+1.1*Re**0.6*Pr**(1.0/3.0)
    h=Nu*kf/dp; av=6.0*(1.0-eps)/dp
    return av*h, Re, Nu
rows=[]
for eps in (0.3,0.6,0.9):
    for G in (1.0,2.0,3.0):
        hv,Re,Nu=hv_wakao(eps,G)
        for q in (1.0e6,2.0e6,3.0e6):
            tag='eps%02d_G%d_q%g'%(int(round(eps*10)),int(G),q/1e6)
            d=os.path.join(outroot,tag); os.makedirs(d,exist_ok=True)
            for f in ('Mesh3d.x','bc3d.inp','material.in'):
                shutil.copy(os.path.join(Tmpl,f),os.path.join(d,f))
            # control.ec: set LS_U_in = G (rho=1)
            c=open(os.path.join(Tmpl,'control.ec')).read()
            c=c.replace('LS_U_in=2.0d0','LS_U_in=%.1fd0'%G)
            open(os.path.join(d,'control.ec'),'w').write(c)
            open(os.path.join(d,'porous.inp'),'w').write('1\n1  %.2f  1.0e-3  %.6e\n'%(eps,hv))
            open(os.path.join(d,'solid_bc.inp'),'w').write('1\n1\n2\n1  0.0  0.0  %.1f  %.1f\n4 -1.0  %.6e   0.0    0.0\n'%(hc,Tinf,q))
            rows.append((tag,eps,G,q,hv,Re,Nu))
import json
json.dump(rows,open(os.path.join(outroot,'cases.json'),'w'),indent=0)
print('created',len(rows),'cases under',outroot)
for r in rows[:3]: print(r)
