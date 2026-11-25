from pathlib import Path
import hashlib,json,subprocess,tomllib
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter

root=Path.cwd();out=root/'rendered/chan-three-case2000-v1'
out.mkdir(parents=True,exist_ok=False)
cases=[('Anelastic','2214011_1','chan-clean-pair-v1-run-2214011_1'),
       ('Compressible · potential temperature','2214011_0','chan-clean-pair-v1-run-2214011_0'),
       ('Compressible · total energy','2237362','chan-total-energy-long-v1-run-2237362')]
times=np.arange(0,2001,5);nframes=len(times)
Q=np.empty((3,nframes,3,256,256),np.float32);C=np.empty_like(Q)
metadata=[];heights=None
for c,(label,job,dirname) in enumerate(cases):
 src=root/'data'/dirname/'parent';meta=tomllib.loads((src/'metadata.toml').read_text())
 assert meta['status']=='COMPLETE' and meta['final_time']==2000
 for k,v in {'nx':256,'ny':256,'nz':80,'omega':3,'beta':0,'jet_speed':0,'bottom_flux':.03125,'cooling_time':20,'seed':314159}.items():assert meta[k]==v,(label,k)
 metadata.append({'label':label,'job':job,'configuration':meta})
 for i,t in enumerate(times):
  with h5py.File(src/f'frame_{int(t*1000):09d}.jld2') as h:
   assert abs(float(h['t'][()])-t)<1e-7
   hs=h['video_heights'][:]
   if heights is None:heights=hs
   assert np.array_equal(hs,heights)
   Q[c,i]=h['vorticity_slices'][:]/meta['f'];C[c,i]=h['tracer_slices'][:]
 print('CACHED',label,flush=True)
assert np.isfinite(Q).all() and np.isfinite(C).all()
qmax=float(np.max(abs(Q)));cmin=float(C.min());cmax=float(C.max())
late_rms=np.std(C[:,300:],axis=(-2,-1),dtype=np.float64)
width=float(np.clip(np.median(late_rms)/3,1e-8,.005))
means=C.mean(axis=(-2,-1),dtype=np.float64);rms=C.std(axis=(-2,-1),dtype=np.float64)
provenance={'cases':metadata,'time_range':[0,2000],'frames':nframes,'fps':12,'interval':5,
 'row_order':[x[0] for x in cases],'heights':heights.tolist(),
 'vorticity':{'variable':'zeta/f0','asinh_width':.02,'limits':[-qmax,qmax],'shared_across':'all cases, depths and times'},
 'tracer':{'variable':'absolute concentration C','center':.1,'asinh_width':width,'limits':[cmin,cmax],'shared_across':'all cases, depths and times','mean_subtraction':False,'per_frame_variance_normalization':False},
 'final_tracer_means':means[:,-1].tolist(),'final_tracer_rms':rms[:,-1].tolist(),
 'notes':['Each movie is3rows by3depth columns. Rows are anelastic, compressible potential temperature, compressible total energy.',
 'Same nondimensional time at every panel; later output from older cases is deliberately excluded from this matched-time comparison.',
 'Thermal states evolve differently; neither matching initial parameters nor visual appearance alone isolates a causal mechanism.',
 'Anelastic and potential-temperature cases use the original solver environment; total energy uses PR1165 pin8273be89d9357ca2ac572347c45db8c647e3c45c and its isolated environment.',
 'All fields and diagnostics are rendered on Unity. Saved tracer fields areFloat32; gradients below stored precision cannot be recovered.']}
(out/'metadata.json').write_text(json.dumps(provenance,indent=2)+'\n')
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':10})
for kind in ['vorticity','tracer']:
 isq=kind=='vorticity'
 fig,axs=plt.subplots(3,3,figsize=(15.6,14),dpi=100)
 fig.subplots_adjust(left=.10,right=.87,bottom=.045,top=.91,wspace=.12,hspace=.25)
 norm=AsinhNorm(linear_width=.02,vmin=-qmax,vmax=qmax) if isq else FuncNorm((lambda x:np.arcsinh((x-.1)/width),lambda x:.1+width*np.sinh(x)),vmin=cmin,vmax=cmax)
 images=[];titles=[]
 for c,(label,_,_) in enumerate(cases):
  for j in range(3):
   ax=axs[c,j];images.append(ax.imshow(np.zeros((256,256)),origin='lower',extent=[0,4,0,4],interpolation='nearest',cmap='RdBu_r' if isq else 'BrBG',norm=norm))
   titles.append(ax.set_title('',fontsize=10));ax.set(xlabel='x/H',ylabel='y/H')
  pos=axs[c,0].get_position();fig.text(.027,.5*(pos.y0+pos.y1),label,rotation=90,ha='center',va='center',fontsize=13,fontweight='bold')
 cax=fig.add_axes([.915,.12,.018,.70]);cb=fig.colorbar(images[0],cax=cax,label='ζ / f₀' if isq else 'Absolute tracer concentration C')
 if not isq:
  ticks=.1+width*np.sinh(np.linspace(np.arcsinh((cmin-.1)/width),np.arcsinh((cmax-.1)/width),7))
  cb.set_ticks(ticks,labels=[f'{x:.6g}' for x in ticks])
 title=fig.suptitle('',y=.98,fontsize=16)
 def update(i):
  title.set_text(('Vorticity' if isq else 'Absolute tracer concentration')+f' · three matched cases · t={times[i]}\nΩ=3, β=0, no jet · shared fixed color scale · '+('positive ζ/f₀ is cyclonic' if isq else 'no mean subtraction or per-frame rescaling'))
  for c in range(3):
   for j in range(3):
    k=3*c+j;images[k].set_data(Q[c,i,j] if isq else C[c,i,j])
    titles[k].set_text(f'z/H={heights[j]:.5f}'+('' if isq else f'\nmean C={means[c,i,j]:.7f}; σ={rms[c,i,j]:.2e}'))
 update(nframes-1);fig.savefig(out/f'{kind}_2000.jpg',dpi=100,pil_kwargs={'quality':90})
 fn=out/f'{kind}_0_2000.mp4';writer=FFMpegWriter(fps=12,codec='libx264',bitrate=2200,extra_args=['-threads','2','-maxrate','2800k','-bufsize','5600k','-pix_fmt','yuv420p','-movflags','+faststart'])
 with writer.saving(fig,str(fn),dpi=100):
  for i in range(nframes):
   update(i);writer.grab_frame()
   if i%100==0:print(kind,'FRAME',i,flush=True)
 plt.close(fig)
 nf,duration=imageio_ffmpeg.count_frames_and_secs(str(fn));assert nf==nframes
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(fn),'-f','null','-'],check=True)
 print(kind,'FULL_DECODE_PASS',fn.stat().st_size,flush=True)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'media.json').write_text(json.dumps({'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames_per_movie':nframes,'fps':12,'time_range':[0,2000],'assets':assets},indent=2)+'\n')
print('ALL_MEDIA_VERIFIED',flush=True)
