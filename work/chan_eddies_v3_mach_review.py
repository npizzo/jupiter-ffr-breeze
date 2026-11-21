import h5py,numpy as np,json,sys,hashlib,tomllib
from pathlib import Path
root=Path(sys.argv[1]);receipt_path=Path(sys.argv[2]);receipt=json.loads(receipt_path.read_text());prior=json.loads((root/'SCREEN_REVIEW.json').read_text());result={};checks=[];lines=['status = "PASS"','review = "Localized-Mach flag adjudicated using saved-field diagnostics; original strict screen remains recorded as FAIL."']
for label in ['case3','case4']:
 r=prior[label];d=root/label;meta=tomllib.loads((d/'metadata.toml').read_text())
 assert meta['status']=='COMPLETE' and meta['model_sha256']==receipt['manifest']['files']['scripts/chan_eddies_v3/model.jl']
 with h5py.File(d/'fields3d/state_0000005.jld2','r') as f:
  rho=f['density'][()].astype(float);T=f['temperature'][()].astype(float);uf=f['u_face'][()].astype(float);vf=f['v_face'][()].astype(float);wf=f['w_face'][()].astype(float)
 u=.5*(uf+np.roll(uf,-1,axis=2));v=.5*(vf+np.roll(vf,-1,axis=1));w=.5*(wf[:-1]+wf[1:])
 mach2=(u*u+v*v+w*w)/(5*T/3);mach=np.sqrt(mach2)
 bulk=float(np.sqrt(np.sum(rho*mach2)/np.sum(rho)))
 fastmass=float(np.sum(rho[mach>.5])/np.sum(rho));fastvolume=float(np.mean(mach>.5))
 vertical=float(np.sum(rho*w*w)/np.sum(rho*(u*u+v*v+w*w)))
 U=np.fft.rfft2(u,axes=(1,2));V=np.fft.rfft2(v,axes=(1,2))
 kx=np.fft.rfftfreq(u.shape[2],4/u.shape[2])[None,None,:];ky=np.fft.fftfreq(u.shape[1],4/u.shape[1])[None,:,None]
 k2=kx*kx+ky*ky;den=np.where(k2>0,k2,1.);D=kx*U+ky*V
 weights=np.ones(U.shape[2])*2;weights[0]=1;weights[-1]=1
 weights=weights[None,None,:]*np.mean(rho,axis=(1,2))[:,None,None]
 divergent=float(np.sum(weights*abs(D)**2/den)/np.sum(weights*(abs(U)**2+abs(V)**2)))
 base=(r['final_time']>=5 and r['max_energy_residual_over_input']<.1 and abs(r['integrated_energy_residual_over_input'])<.1 and r['max_mass_error']<1e-7 and r['barotropic_ke_ratio_min']>.5 and r['barotropic_ke_ratio_max']<2 and r['initialization']['N2_cap_min']>0 and max(abs(r['initialization']['N2_deep_min']),abs(r['initialization']['N2_deep_max']))<.2)
 accepted=bool(base and float(mach.max())<1 and bulk<.15 and fastmass<.001 and divergent<.1 and vertical<.1)
 checks.append(accepted)
 result[label]=dict(accepted=accepted,max_mach=float(mach.max()),density_weighted_rms_mach=bulk,mass_fraction_above_mach_half=fastmass,volume_fraction_above_mach_half=fastvolume,horizontal_divergent_ke_fraction_proxy=divergent,vertical_ke_fraction=vertical,conservation_and_initial_stability_pass=bool(base))
 cp=list((d/'checkpoints').glob('*.jld2'));assert len(cp)==1;h=hashlib.sha256()
 with cp[0].open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 lines.extend([label+'_checkpoint = '+json.dumps(str(cp[0].resolve())),label+'_checkpoint_sha256 = '+json.dumps(h.hexdigest())])
for k in ['mass','kinetic_energy']:checks.append(abs(prior['case3']['initialization'][k]/prior['case4']['initialization'][k]-1)<1e-10)
result['status']='PASS' if all(checks) else 'FAIL'
result['interpretation']='The original max-Mach 0.5 screen flag is retained, not relabeled as passed. This secondary acceptance requires subsonic local flow, low density-weighted RMS Mach, localized fast flow, small horizontal divergent-velocity energy and small vertical kinetic energy, plus unchanged mass/energy/initial-stability gates. The horizontal projection uses layer-mean density and centered velocity; it is a wave-content proxy, not a mode identification. No long-time low-Mach or LSV conclusion is implied.'
(root/'MACH_FIELD_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if not all(checks):sys.exit(2)
lines += ['cpu_receipt_sha256 = '+json.dumps(hashlib.sha256(receipt_path.read_bytes()).hexdigest()),'review_script_sha256 = '+json.dumps(hashlib.sha256(Path(__file__).read_bytes()).hexdigest())]
(root/'PRODUCTION_PASS_FIELD_REVIEW.toml').write_text('\n'.join(lines)+'\n')
