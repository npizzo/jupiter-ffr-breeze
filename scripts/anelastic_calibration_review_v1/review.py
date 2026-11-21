"""Read-only review/render of completed calibration. Execute on Unity only."""
import os, sys, json, csv, hashlib, subprocess
from pathlib import Path
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

root, out = map(Path, sys.argv[1:3]); out.mkdir(parents=True, exist_ok=True)
F=9.375e-7; area=4.; length=2.; depth=.8; nu=2.5e-5; f=.5
rows={}
for sub in ('parent_initial','parent'):
    with open(root/sub/'diagnostics.csv') as h:
        for r in csv.DictReader(h):
            r={k:float(v) for k,v in r.items()}; rows[r['time']]=r
diag=[rows[t] for t in sorted(rows)]
d={k:np.array([r[k] for r in diag]) for k in diag[0]}
frames={}
for sub in ('parent_initial','parent'):
    for p in sorted((root/sub).glob('frame_*.jld2')):
        with h5py.File(p) as h: frames[float(h['t'][()])]=p
times=np.array(sorted(frames)); assert len(times)==401 and times[0]==0 and times[-1]==2000
assert np.allclose(np.diff(times),5)
vmax=0.; cmin=np.inf; cmax=-np.inf; records=[]; spectra=[]; capneg=[]
for t in times:
    with h5py.File(frames[t]) as h:
        zeta=np.array(h['vorticity_slices'],dtype=float); c=h['tracer_slices'][()]
        vmax=max(vmax,float(np.max(np.abs(zeta/f)))); cmin=min(cmin,float(c.min())); cmax=max(cmax,float(c.max()))
        n=zeta.shape[-1]; dx=length/n; k=np.fft.fftfreq(n)*n; kx,ky=np.meshgrid(k,k); kr=np.hypot(kx,ky)
        u=h['ubar'][()].astype(float); v=h['vbar'][()].astype(float)
        uh=np.fft.fft2(u)/n**2; vh=np.fft.fft2(v)/n**2; energy=.5*(abs(uh)**2+abs(vh)**2)
        zb=(np.roll(v,-1,1)-np.roll(v,1,1)-np.roll(u,-1,0)+np.roll(u,1,0))/(2*dx)
        z=h['z'][()]; th=h['theta_profile'][()]; zmid=.5*(z[1:]+z[:-1]); N2=.8*np.diff(np.log(th))/np.diff(z)
        capneg.append(float(np.mean(N2[zmid>.85]<0)))
        rec={'t':t,'barotropic_ke_saved_centered':float(energy.sum()),'gravest_ke':float(energy[np.isclose(kr,1)].sum()),'low_k_le2_ke':float(energy[(kr>0)&(kr<=2)].sum()),'barotropic_skew':float(np.mean(zb**3)/max(np.mean(zb**2)**1.5,1e-99)),'cap_negative_gradient_fraction_above_085':capneg[-1]}
        for j in range(3):
            a=zeta[j]; pos=float(np.maximum(a,0).sum()*dx**2); neg=float(np.minimum(a,0).sum()*dx**2)
            rec.update({f'gamma_plus_{j}':pos,f'gamma_minus_{j}':neg,f'gamma_imbalance_{j}':abs(pos+neg)/max(pos-neg,1e-99),f'skew_{j}':float(np.mean(a**3)/max(np.mean(a*a)**1.5,1e-99))})
        records.append(rec)
        spectra.append(np.bincount(np.floor(kr+.5).astype(int).ravel(),weights=energy.ravel()))
        if t==times[-1]: final={key:h[key][()] for key in ('z','theta_profile','theta_reference_profile','density_profile','video_heights','w_rms_profile','convective_theta_flux_profile')}
with open(out/'mode_and_circulation.csv','w') as h:
    w=csv.DictWriter(h,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
with open(out/'diagnostics.csv','w') as h:
    w=csv.DictWriter(h,fieldnames=list(diag[0])); w.writeheader(); w.writerows(diag)
windows=[]
for start,end in [(0,500),(500,1000),(1000,1500),(1500,2000)]:
    a=np.flatnonzero((d['time']>=start)&(d['time']<=end)); i,j=a[0],a[-1]
    windows.append({'start':start,'end':end,'heat_storage_fraction':float((d['theta_content'][j]-d['theta_content'][i])/(area*F*(end-start))),'urms_mean':float(d['urms_conv'][a].mean()),'wrms_mean':float(d['wrms_conv'][a].mean()),'Ro_spectral_volume_mean':float(d['Ro_spectral_volume'][a].mean())})
summary={'simulation_job':'2058241','review_job':os.getenv('SLURM_JOB_ID'),'stage':'calibration only; no LSV detection or inverse-transfer claim','time_units':'nondimensional, not Earth days','frame_count':len(times),'time_range':[0.,2000.],'geometry':'2 x 2 x 1, horizontally doubly periodic; impermeable stress-free top and bottom','resolution':[128,128,64],'f':f,'bottom_theta_flux':F,'final':diag[-1],'windows':windows,'max_mach':float(d['max_mach'].max()),'max_budget_residual_after_t5_over_input':float(abs(d['budget_residual_over_flux'][d['time']>5]).max()),'max_circulation_relative_imbalance':max(r[f'gamma_imbalance_{j}'] for r in records for j in range(3)),'final_mode_diagnostics':records[-1],'vertical_Re_final':float(d['wrms_conv'][-1]*depth/nu),'video_heights':final['video_heights'].tolist(),'scale':{'vorticity':'zeta/f symmetric asinh, fixed all frames/depths, full extrema','vorticity_limit':vmax,'vorticity_linear_width':.02,'tracer':'absolute concentration, asinh centered at 0.1; tick labels are physical concentrations; fixed full extrema','tracer_linear_width':.01,'tracer_center':.1,'tracer_range':[cmin,cmax]},'limitations':['Thermal storage and cap erosion must be reviewed before production.','Saved centered Float32 vorticity permits circulation checks only to serialization precision.','Barotropic spectra use saved density-weighted depth-mean centered velocity; they are not an energy-transfer budget.','Window means are calibration summaries, not independent statistical ensembles.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2))
plt.rcParams.update({'font.size':9})
fig,ax=plt.subplots(2,3,figsize=(12,7),layout='constrained')
ax[0,0].plot(d['time'],d['urms_conv'],label='total velocity rms'); ax[0,0].plot(d['time'],d['wrms_conv'],label='vertical rms'); ax[0,0].legend(); ax[0,0].set_title('Convection develops')
ax[0,1].plot(d['time'],d['theta_source']/(area*F),label='net source / bottom input'); ax[0,1].plot(d['time'],d['convective_theta_flux_mid']/F,label='mid-layer convective flux / input'); ax[0,1].axhline(0,c='gray',lw=.6); ax[0,1].legend(); ax[0,1].set_title('Thermal balance remains transient')
ax[0,2].plot(final['theta_profile'],final['z'],label='final horizontal mean'); ax[0,2].plot(final['theta_reference_profile'],final['z'],'--',label='reference'); ax[0,2].axhline(.8,c='gray'); ax[0,2].legend(); ax[0,2].set_title('Potential-temperature profile'); ax[0,2].set_ylabel('z')
for key,label in [('barotropic_ke_saved_centered','total barotropic'),('gravest_ke','gravest |k|=1'),('low_k_le2_ke','0<|k|<=2')]: ax[1,0].plot(times,[r[key] for r in records],label=label)
ax[1,0].legend(); ax[1,0].set_title('Absolute depth-mean kinetic energy')
for t in (500,1000,1500,2000):
    e=spectra[list(times).index(t)]; ax[1,1].loglog(np.arange(1,len(e)),e[1:],label=f't={t}')
ax[1,1].legend(); ax[1,1].set_title('Barotropic KE by integer radial shell'); ax[1,1].set_xlabel('k / (2π/L)')
ax[1,2].plot(times,[r['barotropic_skew'] for r in records],label='depth-mean vorticity');
for j,height in enumerate(final['video_heights']): ax[1,2].plot(times,[r[f'skew_{j}'] for r in records],alpha=.65,label=f'z={height:.3f}')
ax[1,2].legend(); ax[1,2].set_title('Vorticity skewness, not circulation')
for a in (ax[0,0],ax[0,1],ax[1,0],ax[1,2]): a.set_xlabel('nondimensional time')
fig.suptitle('Anelastic weak-cap calibration — t=0–2000, not an LSV production run'); fig.savefig(out/'diagnostics.png',dpi=130); plt.close(fig)
ffmpeg=imageio_ffmpeg.get_ffmpeg_exe(); matplotlib.rcParams['animation.ffmpeg_path']=ffmpeg
fig,axs=plt.subplots(2,3,figsize=(10.8,7.2)); fig.subplots_adjust(left=.065,right=.88,bottom=.085,top=.87,wspace=.13,hspace=.2)
vn=AsinhNorm(linear_width=.02,vmin=-vmax,vmax=vmax); cn=AsinhNorm(linear_width=.01,vmin=cmin-.1,vmax=cmax-.1)
ims=[]
for j,height in enumerate(final['video_heights']):
    a=axs[0,j]; im=a.imshow(np.zeros((128,128)),origin='lower',extent=(0,2,0,2),cmap='RdBu_r',norm=vn,interpolation='nearest'); ims.append(im); a.set_title(f'z={height:.4f} above bottom')
for j in range(3): ims.append(axs[1,j].imshow(np.zeros((128,128)),origin='lower',extent=(0,2,0,2),cmap='BrBG',norm=cn,interpolation='nearest'))
for a in axs.flat: a.set_xlabel('x'); a.set_aspect('equal')
axs[0,0].set_ylabel('y — vorticity'); axs[1,0].set_ylabel('y — passive tracer')
cax=fig.add_axes([.9,.535,.016,.315]); cb=fig.colorbar(ims[0],cax=cax); cb.set_label('ζ/f (positive = cyclonic)')
cax=fig.add_axes([.9,.105,.016,.315]); cb=fig.colorbar(ims[3],cax=cax); ticks=np.linspace(cmin,cmax,5); cb.set_ticks(ticks-.1); cb.set_ticklabels([f'{v:.3f}' for v in ticks]); cb.set_label('tracer concentration')
title=fig.suptitle('',fontsize=13); fig.text(.5,.015,'Fixed full-range asinh scales across all depths and times • dry anelastic • f=+0.5 • no β or shear',ha='center',fontsize=9)
writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-crf','27','-preset','fast','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(out/'evolution.mp4'),dpi=100):
    for idx,t in enumerate(times):
        with h5py.File(frames[t]) as h:
            zz=h['vorticity_slices'][()]/f; cc=h['tracer_slices'][()]-.1
        for j in range(3): ims[j].set_data(zz[j]); ims[j+3].set_data(cc[j])
        title.set_text(f'Weak-cap calibration • t={t:.0f} (nondimensional) • 128×128×64')
        writer.grab_frame()
        if t in (0,1000,2000): fig.savefig(out/f'frame-{int(t)}.jpg',dpi=100)
        if idx%50==0: print(f'render {idx+1}/{len(times)}',flush=True)
plt.close(fig)
check=subprocess.run([ffmpeg,'-v','error','-i',str(out/'evolution.mp4'),'-f','null','-'],capture_output=True,text=True)
assert check.returncode==0 and not check.stderr, check.stderr
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
provenance={'input_root_on_unity':str(root),'simulation_job':'2058241','review_job':os.getenv('SLURM_JOB_ID'),'render_script_sha256':sha(Path(__file__)),'first_saved_frame_sha256':sha(frames[times[0]]),'last_saved_frame_sha256':sha(frames[times[-1]]),'movie_sha256':sha(out/'evolution.mp4'),'movie_frames':len(times),'fps':12,'duration_seconds':len(times)/12,'full_movie_decode':'PASS','source_frame_times':'exact saved times 0:5:2000; restart t=5 deduplicated','array_orientation':'HDF5 depth,y,x; displayed origin lower; no interpolation','diagnostic_note':'Cap minimum gradient and net heating are physical transient diagnostics; potential-temperature residual is a separate discrete-time budget check.'}
(out/'provenance.json').write_text(json.dumps(provenance,indent=2)); print(json.dumps(summary),flush=True)
