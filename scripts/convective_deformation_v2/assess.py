from pathlib import Path
import sys,json,tomllib
import numpy as np,h5py
root=Path.cwd();reports=[]
for folder in sys.argv[1:]:
 p=root/folder/'parent';m=tomllib.loads((p/'metadata.toml').read_text());d=np.genfromtxt(p/'diagnostics.csv',delimiter=',',names=True);rows=[]
 for f in sorted(p.glob('full_*.jld2')):
  with h5py.File(f) as h:
   t=float(h['t'][()]);z=h['z'][:];mask=z<.55;rho=h['density'][:];w=.5*(h['w_face'][1:]+h['w_face'][:-1]);T=h['temperature'][:];w=w.astype(float);T=T.astype(float)
   spec=np.sum(abs(np.fft.fft2(w[mask],axes=(-2,-1)))**2,axis=0);n=w.shape[-1];rad=np.hypot(np.fft.fftfreq(n)[:,None]*n,np.fft.fftfreq(n)[None,:]*n);spec[0,0]=0.;den=float(spec.sum());short=float(spec[rad>=n/8].sum()/den) if den else 0.
   shell=np.floor(rad+.5).astype(int);power=np.bincount(shell.ravel(),weights=spec.ravel());peak=int(np.argmax(power[1:])+1);peakcells=float(n/peak)
   heatflux=np.mean(rho*2.5*w*(T-T.mean(axis=(-2,-1),keepdims=True)),axis=(-2,-1));n2=12*np.diff(np.log(h['theta'][:].astype(float)),axis=0)/(1/len(z));cap=(z[:-1]>.7)&(z[:-1]<.9)
   rows.append({'time':t,'convective_w_rms':float(np.sqrt(np.mean(w[mask]**2))),'w_energy_wavelength_le_8dx':short,'peak_w_wavelength_cells':peakcells,'max_lower_convective_enthalpy_flux_over_input':float(np.max(heatflux[mask])/m['initialization']['enthalpy_flux']),'min_interior_cap_mean_N2':float(np.min(n2.mean(axis=(-2,-1))[cap]))})
 final=rows[-1];N=m['initialization']['N2'];U=m['initialization']['U'];late=rows[-2:]
 gates={'completed_10_turnovers':m['status']=='COMPLETE' and m.get('final_turnovers',0)>=9.999,'active_vertical_flow':all(q['convective_w_rms']>.1*U for q in late),'resolved_horizontal_w':all(q['w_energy_wavelength_le_8dx']<.05 and q['peak_w_wavelength_cells']>=8 for q in late),'interior_cap_retained':all(q['min_interior_cap_mean_N2']>.8*N for q in late),'affordable_100_turnover_projection':m['wall_seconds_actual']*10<3600}
 r={'source':folder,'status':m['status'],'scope':'Pilot screening only; passing is not resolution convergence or a proof of convection-driven LSV formation. No automatic extension.', 'gates':gates,'all_screening_passed':all(gates.values()),'snapshots':rows,'initial_N2':N,'initial_U':U,'wall_seconds':m['wall_seconds_actual'],'note':'Cap check is for 0.7<z/H<0.9 below the imposed cooling layer; whole-column numerical guards remain active. Positive convective heat transport is reported separately and must be interpreted with the time history.'}
 reports.append(r)
out=root/'runs/convective-deformation-v2/pilot-assessment.json';out.write_text(json.dumps(reports,indent=2)+'\n');print(json.dumps(reports,indent=2))
