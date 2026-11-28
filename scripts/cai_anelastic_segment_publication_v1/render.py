from pathlib import Path
import csv, json, hashlib, subprocess, shutil, tomllib
import numpy as np
import h5py, imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm, FuncNorm
from matplotlib.animation import FFMpegWriter

root=Path.cwd()
src=Path('/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127/anelastic-B7-2342242/parent')
out=root/'rendered/cai-anelastic-b7-first-segment'
out.mkdir(parents=True,exist_ok=False)
times=[]; qs=[]; cs=[]
for p in sorted(src.glob('frame_*.jld2')):
    with h5py.File(p) as h:
        times.append(float(h['t'][()])); qs.append(h['vorticity_slices'][:]/float(h['f'][()]));cs.append(h['tracer_slices'][:])
        heights=h['video_heights'][:]
Q=np.asarray(qs);C=np.asarray(cs);times=np.asarray(times)
assert list(times)==[0,10,20,30,40] and Q.shape==C.shape and Q.shape[1]==3
assert np.isfinite(Q).all() and np.isfinite(C).all()
with (src/'diagnostics.csv').open() as f: diagnostics=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
outcome=tomllib.loads((src.parent/'OUTCOME.toml').read_text())
summary={'status':'PARTIAL_WALL_LIMIT','job':2342242,'target_time':5000,'checkpoint_time':outcome['final_time'],
         'movie_times':times.tolist(),'grid':[512,512,101],'omega':2.56,'beta':0,'start':'rest with mean-free temperature perturbations, seed314159',
         'final_saved_diagnostics':diagnostics[-1], 'source_directory':str(src),
         'notes':['Normal wall-limit stop with checkpoint; target t=5000 not reached.',
                  'Only five saved times, separated by10 nondimensional units. Playback repeats each saved image without interpolating fields.',
                  'This is an early thermal/convective transient, not a mature LSV or evidence of sign selection.',
                  'The energy proxy is not a conserved anelastic energy; its residual is not interpreted as numerical energy loss.',
                  'At t40 peak Mach~0.58 and max theta/reference departure~0.28 caution against assuming small-perturbation anelastic fidelity.',
                  'Spatial vorticity textures do not establish inverse transfer.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
shutil.copy2(src/'diagnostics.csv',out/'diagnostics.csv')
shutil.copy2(src/'metadata.toml',out/'metadata.toml')
qmax=np.maximum(np.max(abs(Q),axis=(0,2,3)),1e-8)
cmin=float(C.min());cmax=float(C.max());means=C.mean(axis=(-2,-1));stds=C.std(axis=(-2,-1))
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':10})
fig,axs=plt.subplots(2,3,figsize=(16,9.2),dpi=100)
fig.subplots_adjust(left=.05,right=.95,bottom=.07,top=.86,wspace=.42,hspace=.36)
ims=[];labels=[]
for r in range(2):
    for j in range(3):
        norm=AsinhNorm(linear_width=max(qmax[j]/50,1e-4),vmin=-qmax[j],vmax=qmax[j]) if r==0 else FuncNorm((lambda x:np.arcsinh((x-.1)/.01), lambda x:.1+.01*np.sinh(x)),vmin=cmin,vmax=cmax)
        im=axs[r,j].imshow(np.zeros(Q.shape[-2:]),origin='lower',extent=[0,4,0,4],norm=norm,cmap='RdBu_r' if r==0 else 'BrBG',interpolation='nearest')
        ims.append(im);labels.append(axs[r,j].set_title(''))
        axs[r,j].set(xlabel='x / convective depth',ylabel='y / convective depth')
        fig.colorbar(im,ax=axs[r,j],fraction=.046,pad=.03,label='ζ/f₀ (positive = cyclonic)' if r==0 else 'Absolute tracer concentration')
title=fig.suptitle('',y=.97,fontsize=14)
def update(i):
    title.set_text(f'Cai B7 anelastic counterpart · Ω=2.56 · 512² ×101 · t={times[i]:g}\nPartial early transient; saved interval Δt=10; checkpoint t=46.52, target5000\nFixed scales across time; vorticity limits differ by depth; no field interpolation')
    for j in range(3):
        ims[j].set_data(Q[i,j]);labels[j].set_text(f'z/H={heights[j]:.4f}')
        ims[j+3].set_data(C[i,j]);labels[j+3].set_text(f'z/H={heights[j]:.4f} · mean C={means[i,j]:.5f}\nRMS variation={stds[i,j]:.4f}')
update(len(times)-1);fig.savefig(out/'snapshot_t40.jpg',dpi=100,pil_kwargs={'quality':90})
fn=out/'vorticity_tracer_0_40.mp4'
writer=FFMpegWriter(fps=10,codec='libx264',bitrate=1800,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(fn),dpi=100):
    for i in range(len(times)):
        update(i)
        for _ in range(15):writer.grab_frame()
plt.close(fig)
nf,duration=imageio_ffmpeg.count_frames_and_secs(str(fn));assert nf==75
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(fn),'-f','null','-'],check=True)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'media.json').write_text(json.dumps({'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames':75,'duration':duration,'unique_saved_times':times.tolist(),
    'heights':heights.tolist(),'vorticity_limit_per_depth':qmax.tolist(),'vorticity_asinh_width_per_depth':np.maximum(qmax/50,1e-4).tolist(),
    'tracer_absolute_limits':[cmin,cmax],'tracer_mapping':'Asinh centered at0.1 with width0.01 over actual full range; fixed across all times/depths, no clipping or plane-mean subtraction.', 'assets':assets},indent=2)+'\n')
print('ALL_MEDIA_VERIFIED',out,flush=True)
