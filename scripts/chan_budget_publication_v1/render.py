from pathlib import Path
import json,csv,hashlib,subprocess,shutil,tomllib
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter

root=Path.cwd();out=root/'rendered/chan-budget-resolution';out.mkdir(parents=True,exist_ok=False)
summary=json.loads((root/'runs/chan-budget-resolution-v3/analysis.json').read_text())
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
shutil.copy2(root/'data/chan-budget-resolution-v3-cpu-2221330/prepared/TRANSFER_PASS.toml',out/'transfer_audit.toml')
cases=[('baseline','2221768_0','Control · 256² × 80'),('refined','2222246_1','Refined · 384² × 120'),('smaller_dt','2222246_2','Smaller timestep · 256² × 80')]
times=np.arange(200,210.01,.5);Q=[];C=[];Z=[];D=[]
for key,job,label in cases:
 d=root/f'data/chan-budget-resolution-v3-run-{job}/branch'
 assert tomllib.loads((d/'metadata.toml').read_text())['status']=='COMPLETE'
 frames={int(p.stem[6:]):p for p in d.glob('frame_*.jld2')};qs=[];cs=[]
 for t in times:
  with h5py.File(frames[round(t*1000)]) as h:
   assert abs(float(h['t'][()])-t)<1e-7
   qs.append(h['vorticity_slices'][:]/float(h['f'][()]));cs.append(h['tracer_slices'][:]);heights=h['video_heights'][:]
 Q.append(np.stack(qs));C.append(np.stack(cs));Z.append(heights)
 assert np.isfinite(Q[-1]).all() and np.isfinite(C[-1]).all()
 with (d/'diagnostics.csv').open() as f:rows=list(csv.DictReader(f))
 D.append({k:np.array([float(r[k]) for r in rows]) for k in rows[0]})
 shutil.copy2(d/'diagnostics.csv',out/f'{key}_diagnostics.csv')
qmax=max(float(np.max(abs(q))) for q in Q);cmin=min(float(c.min()) for c in C);cmax=max(float(c.max()) for c in C)
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':9,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for (key,job,label),d,style in zip(cases,D,['-','--',':']):
 t=d['time']-200
 axs[0].plot(t,d['residual_step'],style,label=label,lw=2)
 integer=np.abs(t-np.rint(t))<1e-7
 axs[1].plot(t[integer],d['heat_sparse'][integer]-d['heat_step'][integer],style,label=label,lw=2)
axs[0].set(xlabel='Elapsed time from checkpoint',ylabel='Incremental energy residual',title='Refinement does not reduce the deficit')
axs[1].set(xlabel='Elapsed time from checkpoint',ylabel='Once-per-unit heat minus every-step heat',title='Heat sampling error is much smaller')
for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.savefig(out/'budget.svg');plt.close(fig)
for kind in ['vorticity','tracer']:
 fig,axs=plt.subplots(3,3,figsize=(12,12),dpi=100)
 fig.subplots_adjust(left=.065,right=.90,bottom=.055,top=.90,wspace=.30,hspace=.37)
 if kind=='vorticity':norm=AsinhNorm(linear_width=.02,vmin=-qmax,vmax=qmax);cmap='RdBu_r';clabel='ζ/f₀ (positive cyclonic)'
 else:norm=FuncNorm((lambda x:np.arcsinh((x-.1)/.005),lambda x:.1+.005*np.sinh(x)),vmin=cmin,vmax=cmax);cmap='BrBG';clabel='Absolute tracer concentration C'
 ims=[];titles=[]
 for i in range(3):
  for j in range(3):
   im=axs[i,j].imshow(np.zeros((2,2)),origin='lower',extent=[0,4,0,4],interpolation='nearest',cmap=cmap,norm=norm)
   ims.append(im);titles.append(axs[i,j].set_title('',fontsize=9));axs[i,j].set(xlabel='x/H',ylabel='y/H')
 cax=fig.add_axes([.925,.19,.016,.57]);fig.colorbar(ims[0],cax=cax,label=clabel)
 title=fig.suptitle('',y=.975,fontsize=14)
 def update(n):
  for i in range(3):
   for j in range(3):
    a=Q[i][n,j] if kind=='vorticity' else C[i][n,j]
    ims[3*i+j].set_data(a)
    label=f'{cases[i][2]}\nz/H={Z[i][j]:.5f}'
    if kind=='tracer':label+=f' · mean={a.mean():.5f}, σ={a.std():.2g}'
    titles[3*i+j].set_text(label)
  title.set_text(f'Short energy-budget test · {kind}\nt = {times[n]:.1f} · Ω=3, β=0, no jet · fixed shared color scale')
 update(len(times)-1);fig.savefig(out/f'{kind}.jpg',dpi=100,pil_kwargs={'quality':88})
 fn=out/f'{kind}_200_210.mp4';writer=FFMpegWriter(fps=3,codec='libx264',bitrate=1100,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
 with writer.saving(fig,str(fn),dpi=100):
  for n in range(len(times)):update(n);writer.grab_frame()
 plt.close(fig)
 nf,duration=imageio_ffmpeg.count_frames_and_secs(str(fn));assert nf==21 and abs(duration-7)<.1
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(fn),'-f','null','-'],check=True)
 print(kind,'FULL_DECODE_PASS',fn.stat().st_size,flush=True)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
media={'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames_per_movie':21,'fps':3,'times':[200,210],
 'vorticity_scale':{'limits':[-qmax,qmax],'asinh_width':.02},'tracer_scale':{'limits':[cmin,cmax],'center':.1,'asinh_width':.005},
 'row_order':[k for k,j,l in cases],'heights':[z.tolist() for z in Z],'assets':assets}
(out/'media.json').write_text(json.dumps(media,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)
