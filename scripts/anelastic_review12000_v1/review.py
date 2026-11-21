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

root, continuation, continuation2, out = map(Path, sys.argv[1:5]); out.mkdir(parents=True, exist_ok=True)
segments=[root/"parent_initial",root/"parent",continuation/"parent",continuation2/"parent"]
F=9.375e-7; area=4.; length=2.; depth=.8; nu=2.5e-5; f=.5
rows={}; offsets=[]
for segment in segments:
    offset=rows[max(rows)]['turnovers_depth'] if rows else 0.
    offsets.append(offset)
    with open(segment/'diagnostics.csv') as h:
        incoming=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(h)]
    for r in incoming:
        r['turnovers_depth']+=offset
        if r['time'] not in rows: rows[r['time']]=r
assert all(np.isclose(np.diff(sorted(rows)),1))
diag=[rows[t] for t in sorted(rows)]
d={k:np.array([r[k] for r in diag]) for k in diag[0]}
frames={}; duplicate_checks=[]
for segment in segments:
    for p in sorted(segment.glob('frame_*.jld2')):
        with h5py.File(p) as h:
            t=float(h['t'][()])
            if t in frames:
                with h5py.File(frames[t]) as old:
                    diffs={}
                    for key in ('vorticity_slices','tracer_slices','ubar','vbar','theta_profile'):
                        a=old[key][()];b=h[key][()]
                        assert np.allclose(a,b,rtol=1e-6,atol=1e-12),(t,key)
                        diffs[key]=float(np.max(abs(a-b)))
                    duplicate_checks.append({'t':t,'max_abs_differences':diffs})
            else: frames[t]=p
times=np.array(sorted(frames)); assert len(times)==2401 and times[0]==0 and times[-1]==12000
assert np.allclose(np.diff(times),5)
print('Validated complete saved timeline and restart joins',flush=True)
vmax=0.; cmin=np.inf; cmax=-np.inf; records=[]; spectra=[]; capneg=[]; snapshots={}; parseval_errors=[]
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
        lowmask=(kr>0)&(kr<=2)
        zlow=np.fft.ifft2(np.fft.fft2(zb)*lowmask).real
        expected=rows[t]['barotropic_ke'];parseval_errors.append(abs(energy.sum()-expected)/max(expected,1e-99))
        assert parseval_errors[-1]<1e-6,(t,parseval_errors[-1])
        assert np.all(np.isfinite(zeta)) and np.all(np.isfinite(c))
        if t in (6000,9000,12000): snapshots[t]=(zb.copy(),zlow.copy(),z.copy(),th.copy(),N2.copy(),zmid.copy())
        rec={'t':t,'barotropic_ke_saved_centered':float(energy.sum()),'gravest_ke':float(energy[np.isclose(kr,1)].sum()),'low_k_le2_ke':float(energy[(kr>0)&(kr<=2)].sum()),'barotropic_skew':float(np.mean(zb**3)/max(np.mean(zb**2)**1.5,1e-99)),'cap_negative_gradient_fraction_above_085':capneg[-1],'low_k_zeta_skew':float(np.mean(zlow**3)/max(np.mean(zlow*zlow)**1.5,1e-99)),'gravest_ke_fraction':float(energy[np.isclose(kr,1)].sum()/max(energy.sum(),1e-99)),'low_k_le2_ke_fraction':float(energy[lowmask].sum()/max(energy.sum(),1e-99))}
        for j in range(3):
            a=zeta[j]; pos=float(np.maximum(a,0).sum()*dx**2); neg=float(np.minimum(a,0).sum()*dx**2)
            rec.update({f'gamma_plus_{j}':pos,f'gamma_minus_{j}':neg,f'gamma_imbalance_{j}':abs(pos+neg)/max(pos-neg,1e-99),f'skew_{j}':float(np.mean(a**3)/max(np.mean(a*a)**1.5,1e-99))})
        for j in range(3):
            alow=np.fft.ifft2(np.fft.fft2(zeta[j])*lowmask).real
            rec[f'low_k_depth_coherence_{j}']=float(np.mean(alow*zlow)/max(np.sqrt(np.mean(alow**2)*np.mean(zlow**2)),1e-99))
        records.append(rec)
        spectra.append(np.bincount(np.floor(kr+.5).astype(int).ravel(),weights=energy.ravel()))
        if t==times[-1]: final={key:h[key][()] for key in ('z','theta_profile','theta_reference_profile','density_profile','video_heights','w_rms_profile','convective_theta_flux_profile')}
with open(out/'mode_and_circulation.csv','w') as h:
    w=csv.DictWriter(h,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
with open(out/'diagnostics.csv','w') as h:
    w=csv.DictWriter(h,fieldnames=list(diag[0])); w.writeheader(); w.writerows(diag)
windows=[]
for start,end in [(k,k+500) for k in range(0,12000,500)]:
    a=np.flatnonzero((d['time']>=start)&(d['time']<=end)); i,j=a[0],a[-1]
    windows.append({'start':start,'end':end,'heat_storage_fraction':float((d['theta_content'][j]-d['theta_content'][i])/(area*F*(end-start))),'urms_mean':float(d['urms_conv'][a].mean()),'wrms_mean':float(d['wrms_conv'][a].mean()),'Ro_spectral_volume_mean':float(d['Ro_spectral_volume'][a].mean()),'mid_convective_flux_over_input_mean':float(d['convective_theta_flux_mid'][a].mean()/F),'cap_min_gradient_mean':float(d['cap_min_dtheta_dz'][a].mean())})
summary={'simulation_jobs':['2058241','2071578','2077087'],'review_job':os.getenv('SLURM_JOB_ID'),'stage':'calibration only; no LSV detection or inverse-transfer claim','time_units':'nondimensional, not Earth days','frame_count':len(times),'time_range':[0.,12000.],'geometry':'2 x 2 x 1, horizontally doubly periodic; impermeable stress-free top and bottom','resolution':[128,128,64],'f':f,'bottom_theta_flux':F,'final':diag[-1],'windows':windows,'max_mach':float(d['max_mach'].max()),'max_budget_residual_after_t5_over_input':float(abs(d['budget_residual_over_flux'][d['time']>5]).max()),'max_circulation_relative_imbalance':max(r[f'gamma_imbalance_{j}'] for r in records for j in range(3)),'final_mode_diagnostics':records[-1],'vertical_Re_final':float(d['wrms_conv'][-1]*depth/nu),'final_tracer_range':[float(c.min()),float(c.max())],'video_heights':final['video_heights'].tolist(),'scale':{'vorticity':'zeta/f symmetric asinh, fixed all frames/depths, full extrema','vorticity_limit':vmax,'vorticity_linear_width':.02,'tracer':'absolute concentration, asinh centered at 0.1; tick labels are physical concentrations; fixed full extrema','tracer_linear_width':.01,'tracer_center':.1,'tracer_range':[cmin,cmax]},'limitations':['Thermal storage and cap erosion must be reviewed before production.','Saved centered Float32 vorticity permits circulation checks only to serialization precision.','Barotropic spectra use saved density-weighted depth-mean centered velocity; they are not an energy-transfer budget.','Window means are calibration summaries, not independent statistical ensembles.']}
summary['mode_windows']=[]
for start,end in [(1500,2000)]+[(k,k+1000) for k in range(2000,12000,1000)]:
    rr=[r for r in records if start<=r['t']<=end]
    keys=['barotropic_ke_saved_centered','gravest_ke','gravest_ke_fraction','low_k_le2_ke','low_k_le2_ke_fraction','barotropic_skew','low_k_zeta_skew','cap_negative_gradient_fraction_above_085']+[f'low_k_depth_coherence_{j}' for j in range(3)]+[f'skew_{j}' for j in range(3)]
    summary['mode_windows'].append({'start':start,'end':end,'means':{k:float(np.mean([r[k] for r in rr])) for k in keys},'std':{k:float(np.std([r[k] for r in rr])) for k in keys}})
summary['milestones']={str(t):next(r for r in records if r['t']==t) for t in (2000,4000,6000,9000,12000)}
summary['validation']={'max_spectral_vs_saved_barotropic_ke_relative_error':max(parseval_errors),'duplicate_frame_checks':duplicate_checks,'turnover_offsets_per_segment':offsets,'all_2401_times_present':True}
(out/'summary.json').write_text(json.dumps(summary,indent=2))
fig,ax=plt.subplots(2,3,figsize=(11,7),layout='constrained')
limits=[max(np.max(abs(snapshots[t][i]/f)) for t in snapshots) for i in (0,1)]
for j,t in enumerate((6000,9000,12000)):
    for i in (0,1):
        im=ax[i,j].imshow(snapshots[t][i]/f,origin='lower',extent=(0,2,0,2),cmap='RdBu_r',norm=AsinhNorm(linear_width=.005,vmin=-limits[i],vmax=limits[i]),interpolation='nearest');ax[i,j].set_title(f't={t}');ax[i,j].set_xlabel('x')
        if j==2:fig.colorbar(im,ax=ax[i,:].tolist(),label='ζ / f')
ax[0,0].set_ylabel('Full depth-mean vorticity · y');ax[1,0].set_ylabel('Only 0 < |k| ≤ 2 · y')
fig.suptitle('Growth of large scales: unfiltered and explicitly low-pass depth-mean flow')
fig.savefig(out/'large_scales.png',dpi=130);plt.close(fig)
plt.rcParams.update({'font.size':9})
fig,ax=plt.subplots(2,3,figsize=(12,7),layout='constrained')
ax[0,0].plot(d['time'],d['urms_conv'],label='total velocity rms'); ax[0,0].plot(d['time'],d['wrms_conv'],label='vertical rms'); ax[0,0].legend(); ax[0,0].set_title('Convection develops')
ax[0,1].plot(d['time'],d['theta_source']/(area*F),label='net source / bottom input'); ax[0,1].plot(d['time'],d['convective_theta_flux_mid']/F,label='mid-layer convective flux / input'); ax[0,1].axhline(0,c='gray',lw=.6); ax[0,1].legend(); ax[0,1].set_title('Thermal balance remains transient')
ax[0,2].plot(final['theta_profile'],final['z'],label='final horizontal mean'); ax[0,2].plot(final['theta_reference_profile'],final['z'],'--',label='reference'); ax[0,2].axhline(.8,c='gray'); ax[0,2].legend(); ax[0,2].set_title('Potential-temperature profile'); ax[0,2].set_ylabel('z')
for key,label in [('barotropic_ke_saved_centered','total barotropic'),('gravest_ke','gravest |k|=1'),('low_k_le2_ke','0<|k|<=2')]: ax[1,0].plot(times,[r[key] for r in records],label=label)
ax[1,0].legend(); ax[1,0].set_title('Absolute depth-mean kinetic energy')
for t in (2000,6000,9000,12000):
    e=spectra[list(times).index(t)]; ax[1,1].loglog(np.arange(1,len(e)),e[1:],label=f't={t}')
ax[1,1].legend(); ax[1,1].set_title('Barotropic KE by integer radial shell'); ax[1,1].set_xlabel('k / (2π/L)')
ax[1,2].plot(times,[r['barotropic_skew'] for r in records],label='depth-mean vorticity');
for j,height in enumerate(final['video_heights']): ax[1,2].plot(times,[r[f'skew_{j}'] for r in records],alpha=.65,label=f'z={height:.3f}')
ax[1,2].legend(); ax[1,2].set_title('Vorticity skewness, not circulation')
for a in (ax[0,0],ax[0,1],ax[1,0],ax[1,2]): a.set_xlabel('nondimensional time')
fig.suptitle('Anelastic weak-cap calibration — t=0–12000, unchanged weak-cap configuration'); fig.savefig(out/'diagnostics.png',dpi=130); plt.close(fig)
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
with writer.saving(fig,str(out/'evolution-master.mp4'),dpi=100):
    for idx,t in enumerate(times):
        with h5py.File(frames[t]) as h:
            zz=h['vorticity_slices'][()]/f; cc=h['tracer_slices'][()]-.1
        for j in range(3): ims[j].set_data(zz[j]); ims[j+3].set_data(cc[j])
        title.set_text(f'Weak-cap continuation • t={t:.0f} (nondimensional) • 128×128×64')
        writer.grab_frame()
        if t in (0,6000,9000,12000): fig.savefig(out/f'frame-{int(t)}.jpg',dpi=100)
        if idx%50==0: print(f'render {idx+1}/{len(times)}',flush=True)
plt.close(fig)
# Two-pass web delivery preserves resolution, cadence and every frame, within browser upload limits.
master=out/'evolution-master.mp4'; movie=out/'evolution.mp4'
base=[ffmpeg,'-y','-v','error','-i',str(master),'-map','0:v:0','-c:v','libx264','-b:v','900k','-preset','slow','-threads','4','-pix_fmt','yuv420p','-an','-passlogfile',str(out/'encode-pass')]
subprocess.run(base+['-pass','1','-f','null','/dev/null'],check=True)
subprocess.run(base+['-pass','2','-movflags','+faststart',str(movie)],check=True)
assert movie.stat().st_size<25_000_000, movie.stat().st_size
for sec in (100,199.9):
    subprocess.run([ffmpeg,'-y','-v','error','-ss',str(sec),'-i',str(movie),'-frames:v','1',str(out/f'encoded-{sec}.jpg')],check=True)
check=subprocess.run([ffmpeg,'-v','error','-i',str(out/'evolution.mp4'),'-f','null','-'],capture_output=True,text=True)
assert check.returncode==0 and not check.stderr, check.stderr
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
provenance={'input_roots_on_unity':[str(root),str(continuation),str(continuation2)],'simulation_jobs':['2058241','2071578','2077087'],'review_job':os.getenv('SLURM_JOB_ID'),'render_script_sha256':sha(Path(__file__)),'first_saved_frame_sha256':sha(frames[times[0]]),'last_saved_frame_sha256':sha(frames[times[-1]]),'movie_sha256':sha(out/'evolution.mp4'),'movie_frames':len(times),'fps':12,'duration_seconds':len(times)/12,'full_movie_decode':'PASS','source_frame_times':'exact saved times 0:5:12000; any overlapping saved frames checked; actual restart equality recorded by simulation pickup receipts','array_orientation':'HDF5 depth,y,x; displayed origin lower; no interpolation','diagnostic_note':'Cap minimum gradient and net heating are physical transient diagnostics; potential-temperature residual is a separate discrete-time budget check.'}
provenance['duplicate_frame_checks']=duplicate_checks
provenance['all_frames_finite']=True
provenance['web_encoding']={'codec':'H264','resolution':[1080,720],'target_video_bitrate':900000,'passes':2,'all_saved_frames_retained':True,'bytes':movie.stat().st_size}
(out/'provenance.json').write_text(json.dumps(provenance,indent=2)); print('REVIEW_AND_MOVIE_PASS '+str(out),flush=True)
