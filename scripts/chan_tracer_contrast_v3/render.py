from pathlib import Path
import ast,json,hashlib,subprocess,sys
import h5py,numpy as np,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
root=Path.cwd();case=sys.argv[1]
tree=ast.parse((root/'scripts/chan_nov26_completed_v2/render.py').read_text());cases=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='cases' for t in n.targets));cfg=cases[case]
prior=root/'rendered'/('chan-nov26-'+case);m=json.loads((prior/'media.json').read_text());out=root/'rendered'/('chan-contrast-v3-'+case);out.mkdir(exist_ok=False)
times=np.arange(0,cfg['end']+1,5);C=[]
for t in times:
 part=next(root/'data'/p/'parent' for stop,p in cfg['sources'] if t<=stop)
 with h5py.File(part/f'frame_{t*1000:09d}.jld2') as h:C.append(h['tracer_slices'][:]);heights=h['video_heights'][:]
C=np.array(C);means=C.mean(axis=(-2,-1),dtype=float);rms=C.std(axis=(-2,-1),dtype=float);A=C-means[:,:,None,None];lim=float(abs(A).max());width=m['tracer_contrast']['asinh_width'];inner=min(lim/4,10**round(np.log10(4*width)));ticks=[-lim,-inner,0,inner,lim]
fig,axs=plt.subplots(1,3,figsize=(16,5.3),dpi=100);fig.subplots_adjust(left=.05,right=.93,bottom=.13,top=.75,wspace=.6);ims=[]
for j,ax in enumerate(axs):
 im=ax.imshow(A[0,j],origin='lower',extent=[0,4,0,4],interpolation='nearest',cmap='BrBG',norm=AsinhNorm(linear_width=width,vmin=-lim,vmax=lim));ims.append(im);ax.set(xlabel='x/H',ylabel='y/H');cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.03,label='C − plane mean');cb.set_ticks(ticks,labels=[f'{v:.2g}' for v in ticks]);cb.ax.minorticks_off()
title=fig.suptitle('',fontsize=13)
movie=out/f'tracer_contrast_0_{cfg["end"]}.mp4';writer=imageio_ffmpeg.write_frames(str(movie),(1600,530),fps=20,codec='libx264',bitrate='1800k',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=2,output_params=['-threads','2','-movflags','+faststart']);writer.send(None)
for i,t in enumerate(times):
 title.set_text(f'{cfg["title"]} · Ω=3, β=0, no jet · t={t}\n256² × 80 · C minus instantaneous plane mean · fixed scale, no RMS normalization')
 for j,im in enumerate(ims):im.set_data(A[i,j]);axs[j].set_title(f'z/H={heights[j]:.5f}\nmean C={means[i,j]:.8f}; RMS variation={rms[i,j]:.2e}')
 fig.canvas.draw();writer.send(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
writer.close();fig.savefig(out/f'tracer_contrast_{cfg["end"]}.jpg',dpi=100);plt.close(fig)
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==len(times)
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
r={'status':'MEDIA_VERIFIED','frames':frames,'duration':duration,'full_decode_passed':True,'contrast_limits':[-lim,lim],'asinh_width':width,'ticks':ticks,'definition':'C minus instantaneous horizontal-plane mean. Fixed full-data range and fixed asinh width; no RMS normalization. Actual mean/RMS printed. Five explicit colorbar ticks to avoid overlap.','assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}}
(out/'contrast_media.json').write_text(json.dumps(r,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)
