"""Unity-only resolution/equilibration review and synchronized three-height movies."""
import csv,json,sys,tomllib,subprocess,hashlib
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import imageio_ffmpeg
out=Path(sys.argv[1]);roots=[Path(x) for x in sys.argv[2:]];out.mkdir(parents=True,exist_ok=False)
metas=[tomllib.loads((p/'metadata.toml').read_text()) for p in roots]
assert len(roots)==3 and all(m['status'] in ('COMPLETE','WALL_LIMIT') for m in metas)
limits=np.zeros(3); fscale=1.28; files=[]; metrics=[]
for root,m in zip(roots,metas):
 frames=sorted(root.glob('frame_*.jld2'));assert len(frames)>=2
 tt=[]
 for path in frames:
  with h5py.File(path) as d:
   tt.append(float(d['t'][()])); a=d['vorticity_slices'][:]/fscale
   assert a.shape[0]==3 and np.isfinite(a).all()
   limits=np.maximum(limits,np.max(abs(a),axis=(1,2)))
 assert np.all(np.diff(tt)>0)
 files.append(frames)
 for path in sorted((root/'fields3d').glob('*.jld2')):
  with h5py.File(path) as d:
   w=.5*(d['w_face'][:-1].astype(float)+d['w_face'][1:].astype(float));t=float(d['t'][()]);rho=d['density'][:].astype(float);u=d['u_face'][:].astype(float);v=d['v_face'][:].astype(float)
  nz,ny,nx=w.shape;w-=w.mean(axis=(1,2),keepdims=True)
  E=abs(np.fft.rfft2(w,axes=(1,2),norm='forward'))**2;fac=np.ones(nx//2+1);fac[1:-1]=2;E*=fac[None,None,:]
  kx=2*np.pi*np.fft.rfftfreq(nx,4/nx)[None,:];ky=2*np.pi*np.fft.fftfreq(ny,4/ny)[:,None];K=np.hypot(kx,ky)
  e=E.mean(axis=0);e[0,0]=0;kmean=np.sum(e*K)/max(e.sum(),1e-300)
  metrics.append(dict(case=m['name'],t=t,w_rms=float(np.sqrt(np.mean(w*w))),gradient_length=float(1/max(kmean,1e-300)),wavelength=float(2*np.pi/max(kmean,1e-300)),fraction_w_energy_shorter_than_6dx=float(e[K>2*np.pi/(6*4/nx)].sum()/max(e.sum(),1e-300))))
(out/'resolution.json').write_text(json.dumps(metrics,indent=2)+'\n')
summary=[];assets=[]
ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
limits=np.maximum(limits,1e-10)
for root,m,frames in zip(roots,metas,files):
 fig,axs=plt.subplots(2,3,figsize=(15,7.4),layout='constrained',dpi=100)
 images=[]
 with h5py.File(frames[0]) as d:heights=d['video_heights'][:]
 for j in range(3):
  for row in range(2):
   norm=AsinhNorm(linear_width=.1,vmin=-limits[j],vmax=limits[j]) if row==0 else AsinhNorm(linear_width=.01,vmin=0,vmax=1)
   img=axs[row,j].imshow(np.zeros((m['ny'],m['nx'])),origin='lower',extent=(0,4,0,4),cmap='RdBu_r' if row==0 else 'viridis',norm=norm,aspect='equal')
   axs[row,j].set(xlabel='x/H',ylabel='y/H',title=f"{'Vorticity' if row==0 else 'Tracer'} | z/H={heights[j]:.4f}")
   fig.colorbar(img,ax=axs[row,j],label='ζ / 1.28 (common benchmark scale)' if row==0 else 'C',shrink=.85)
   images.append((row,j,img))
 title=fig.suptitle('')
 movie=out/(m['name']+'.mp4');poster=out/(m['name']+'.png')
 writer=imageio_ffmpeg.write_frames(str(movie),size=(1500,740),fps=12,codec='libx264',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',quality=7,macro_block_size=2,output_params=['-crf','24','-preset','fast','-movflags','+faststart','-threads','1']);writer.send(None)
 times=[]
 for path in frames:
  with h5py.File(path) as d:t=float(d['t'][()]);vort=d['vorticity_slices'][:]/fscale;dye=d['tracer_slices'][:]
  assert np.isfinite(dye).all();times.append(t)
  for row,j,img in images:img.set_data((vort if row==0 else dye)[j])
  title.set_text(f"Preliminary Cai benchmark: {m['name']} | t={t:.1f} | t F={t*m['flux']:.2f}\n128 × 128 × 64; Ω={m['omega']}; {m['status']} — not a GRS reproduction")
  fig.canvas.draw();writer.send(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
 writer.close();fig.savefig(poster,dpi=100);plt.close(fig)
 subprocess.run([ffmpeg,'-v','error','-i',str(movie),'-f','null','-'],check=True)
 diag=list(csv.DictReader((root/'diagnostics.csv').open()));late=[r for r in diag if float(r['time'])>=.8*m['final_time']]
 summary.append(dict(case=m['name'],status=m['status'],final_time=m['final_time'],nominal_thermal_times=m['final_time']*m['flux'],late_velocity_rms=float(np.mean([float(r['rms_velocity']) for r in late])),late_top_flux_over_bottom=float(np.mean([float(r['top_heat_flux']) for r in late]))/m['flux'],note='Late means alone do not establish equilibrium. Low resolution; no automatic LSV conclusion.'))
 info=dict(case=m['name'],source=str(root.resolve()),frames=len(frames),times=[times[0],times[-1]],video_heights=heights.tolist(),vorticity_common_reference=1.28,vorticity_limits=limits.tolist(),vorticity_asinh_width=.1,tracer_limits=[0,1],tracer_asinh_width=.01,decode='PASS',metadata=m)
 (out/(m['name']+'.json')).write_text(json.dumps(info,indent=2)+'\n')
 assets += [movie,poster,out/(m['name']+'.json')]
nonrot=next(s for s in summary if s['case']=='B0_nonrotating')
for s,m in zip(summary,metas):s['late_reference_Roc_provisional']=nonrot['late_velocity_rms']/(2*m['omega']) if m['omega'] else None
(out/'summary.json').write_text(json.dumps(dict(cases=summary,interpretation='Provisional Rossby uses actual nonrotating velocity; only meaningful after reference equilibrium. Inspect resolution.json before interpreting absence or sign of LSVs.'),indent=2)+'\n')
assets += [out/'summary.json',out/'resolution.json']
(out/'publication.json').write_text(json.dumps(dict(status='RENDERED_AND_DECODE_VERIFIED',assets=[dict(name=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in assets]),indent=2)+'\n')
print('CAI_BENCHMARK_REVIEW_COMPLETE')
