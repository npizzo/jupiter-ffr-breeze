from pathlib import Path
import json,csv,h5py,numpy as np
root=Path('/work/pi_nicholas_pizzo_uri_edu/jupiter-ffr-breeze');result={}
for label,omega,job in [('case3',1.5,'2089750'),('case4',3.,'2089751')]:
 screen=root/'data/chan-eddies-v3-screen-2089476'/label;parent=root/('data/chan-eddies-v3-'+label+'-'+job+'/parent')
 rows=[]
 for d in [screen,parent]:rows += [{k:float(v) for k,v in r.items()} for r in csv.DictReader((d/'diagnostics.csv').open())]
 rows.sort(key=lambda r:r['time']);a,b=rows[0],rows[-1]
 heat=sum(.5*(x['external_power']+y['external_power'])*(y['time']-x['time']) for x,y in zip(rows,rows[1:]));defect=b['total_energy']-a['total_energy']-heat
 states=[]
 for p in sorted((screen/'fields3d').glob('*.jld2')):
  with h5py.File(p,'r') as f:
   rho=f['density'][()].astype(float);u=f['u_face'][()].astype(float);v=f['v_face'][()].astype(float);w=f['w_face'][()].astype(float);t=float(f['t'][()])
  rx=.5*(rho+np.roll(rho,1,axis=2));ry=.5*(rho+np.roll(rho,1,axis=1));rz=.5*(rho[1:]+rho[:-1]);mx=rx*u;my=ry*v;mw=np.zeros_like(w);mw[1:-1]=rz*w[1:-1]
  iv=.25*(my+np.roll(my,1,axis=2)+np.roll(my,-1,axis=1)+np.roll(np.roll(my,1,axis=2),-1,axis=1))
  iu=.25*(mx+np.roll(mx,-1,axis=2)+np.roll(mx,1,axis=1)+np.roll(np.roll(mx,-1,axis=2),1,axis=1))
  mwc=.5*(mw[:-1]+mw[1:]);iw=.5*(mwc+np.roll(mwc,1,axis=2))
  mxc=.5*(mx+np.roll(mx,-1,axis=2));ix=.5*(mxc[:-1]+mxc[1:])
  fz=2*omega*np.sin(np.deg2rad(67.5));fy=2*omega*np.cos(np.deg2rad(67.5));dv=16/rho.size
  horiz=float(np.sum(fz*(u*iv-v*iu))*dv);tilt=float((-np.sum(fy*u*iw)+np.sum(fy*w[1:-1]*ix))*dv)
  ke=float(.5*(np.sum(rx*u*u)+np.sum(ry*v*v)+np.sum(rz*w[1:-1]**2))*dv)
  states.append(dict(time=t,kinetic_energy=ke,coriolis_power_horizontal=horiz,coriolis_power_tilted=tilt,coriolis_power_over_heating=(horiz+tilt)/.5))
 result[label]=dict(first_time=a['time'],last_time=b['time'],integrated_external_heat=heat,total_energy_change=b['total_energy']-a['total_energy'],energy_defect=defect,defect_over_initial_kinetic_energy=defect/states[0]['kinetic_energy'],defect_over_integrated_heat=defect/heat,last_power_defect=.5*b['energy_residual_over_input'],states=states)
result['limitations']='Coriolis-work estimate uses saved Float32 face velocities/density and the constant-Cartesian staggered interpolation. It excludes time-discretization errors and is not a full energy-operator decomposition.'
(root/'data/chan-eddies-v3-screen-2089476/ENERGY_EXPLANATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
