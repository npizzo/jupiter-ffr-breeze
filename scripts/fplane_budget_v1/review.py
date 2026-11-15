"""Audit exact native balance, sparse reconstruction, and matched-grid sensitivity on Unity."""
from pathlib import Path
import sys,json,tomllib,subprocess,hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio_ffmpeg
from gate import check_case
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=False,parents=True)
result={};common={}
def load_native(p):
 with h5py.File(p) as h:return {k:np.asarray(h[k]).T for k in ['Z','origin_Z','cumulative_parts','cumulative_residual','third_moment_storage','third_moment_parts','third_moment_residual']}
def rms(x,mask):return float(np.sqrt(np.mean(x[mask]**2)))
for name in ['mature512','mature1024']:
 d=root/name;meta=tomllib.loads((d/'metadata.toml').read_text());assert meta['status']=='COMPLETE'
 native_check=check_case(d);start=90.5*86400;end=90.75*86400
 files=sorted((d/'reduced').glob('state_*.npz'));files=[p for p in files if start<=int(p.stem.split('_')[1])<=end];assert len(files)==37
 rows=[json.loads(p.with_suffix('.json').read_text()) for p in files]
 with np.load(files[-1]) as last:z=last['Z'].copy();y=last['y'].copy()
 mask=np.broadcast_to(np.abs(y)[None,:]<=.4*meta['ly'],z.shape)
 n0=load_native(d/'native_vorticity'/f'budget_{round(start):012d}.jld2');n1=load_native(d/'native_vorticity'/f'budget_{round(end):012d}.jld2')
 exactparts=(n1['cumulative_parts']-n0['cumulative_parts'])/(end-start);storage=(n1['Z']-n0['Z'])/(end-start)
 exact_remainder=storage-exactparts.sum(2)
 temporal={}
 for stride in [1,2,6]:
  selected=files[::stride];pred=np.zeros_like(z);previous=None
  for p in selected:
   with np.load(p) as v:
    q=sum(v[k] for k in ['vertical_source','horizontal_source','barotropic_advection','mean_divergence_source']);t=float(v['time_days'])*86400
   if previous is not None:pred+=(q+previous[1])/2*(t-previous[0])
   previous=(t,q)
  pred/=(end-start)
  mismatch=exactparts[:,:,0]+exactparts[:,:,1]-pred
  temporal[str(stride*10)]={'sample_minutes':stride*10,'reconstruction_rms':rms(pred,mask),'native_advection_plus_coriolis_rms':rms(exactparts[:,:,:2].sum(2),mask),'reconstruction_mismatch_rms':rms(mismatch,mask),'reconstruction_mismatch_over_native_adv_cor':rms(mismatch,mask)/max(rms(exactparts[:,:,:2].sum(2),mask),1e-30)}
  if stride==1:finest=pred.copy()
  else:temporal[str(stride*10)]['difference_from_10min_over_10min_rms']=rms(pred-finest,mask)/max(rms(finest,mask),1e-30)
 case={'native_check':native_check,'window_days':[90.5,90.75],'grid':[meta['nx'],meta['ny'],meta['nz']],'native_storage_rms':rms(storage,mask),'native_closure_rms':rms(exact_remainder,mask),'native_closure_over_storage':rms(exact_remainder,mask)/max(rms(storage,mask),1e-30),'term_rms':dict(zip(['advection','coriolis','explicit_diffusion','external_forcing'],[rms(exactparts[:,:,i],mask) for i in range(4)])),'temporal':temporal,'third_moment_storage':(n1['third_moment_storage']-n0['third_moment_storage']).tolist(),'third_moment_native_parts':(n1['third_moment_parts']-n0['third_moment_parts']).tolist(),'final_moments':rows[-1]['moments']['interior'],'remap':tomllib.loads((d/'remap_validation.toml').read_text()),'wall_seconds':meta['wall_seconds']}
 result[name]=case
 def avg2(a):return (a[::2,::2]+a[1::2,::2]+a[::2,1::2]+a[1::2,1::2])/4
 common[name]={'Z':z if meta['nx']==512 else avg2(z),'storage':storage if meta['nx']==512 else avg2(storage),'advection_coriolis':exactparts[:,:,:2].sum(2) if meta['nx']==512 else avg2(exactparts[:,:,:2].sum(2))}
y=(np.arange(512)+.5)*6000000/512-3000000;mask=np.broadcast_to(abs(y)[None,:]<=2400000,(512,512));differences={}
for key in common['mature512']:
 a=common['mature512'][key];b=common['mature1024'][key];differences[key]={'coarse_rms':rms(a,mask),'restricted_fine_rms':rms(b,mask),'difference_rms_over_coarse':rms(b-a,mask)/max(rms(a,mask),1e-30)}
fig,ax=plt.subplots(1,3,figsize=(15,4.8),constrained_layout=True)
for n in result:
 c=result[n];mins=[10,20,60];ax[0].plot(mins,[c['temporal'][str(x)]['reconstruction_mismatch_over_native_adv_cor'] for x in mins],marker='o',label=n)
ax[0].set(xlabel='Centered-source sampling interval [min]',ylabel='Mismatch / native advection + Coriolis RMS',title='Reconstruction and output cadence');ax[0].legend()
ax[1].bar(list(result),[result[n]['native_closure_over_storage'] for n in result]);ax[1].set(yscale='log',ylabel='Native residual / storage RMS',title='Every-stage native budget closure')
ax[2].bar(list(differences),[v['difference_rms_over_coarse'] for v in differences.values()]);ax[2].set(ylabel='Fine–coarse difference / coarse RMS',title='Common 512²-grid sensitivity');ax[2].tick_params(axis='x',rotation=25)
fig.savefig(out/'native_budget_audit.png',dpi=150);plt.close(fig)
movie=out/'mature512_over_mature1024.mp4'
subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),str(root/'mature512'),str(root/'mature1024'),'--labels','512² native-budget restart','1024² native-budget restart','--title','Mature heated f-plane | dense budget audit | days 90–90.75','--tracer','--fps','4','--crf','28','--output',str(movie)],check=True)
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True);assert movie.stat().st_size<24000000
report={'status':'COMPLETE','compute_location':'Unity','cases':result,'common_grid_differences':differences,'limits':['This is an 18-hour remapped restart, not long-time convergence or independent statistical realizations.','Only the last six hours enter the dense-window comparison; twelve hours may not remove all remap adjustment.','The native explicit advection term includes the actual WENO discretization; its difference from a centered reconstruction is not automatically physical dissipation.','Every-stage native integrals remove sampling ambiguity from the native balance; the V/H continuum decomposition may still need finer temporal and spatial sampling.','Vertical resolution, boundary and cap sensitivity are not solved by horizontal refinement.'],'recommendation':'Review native closure, cadence and common-grid differences before choosing longer continuations or a one-minute source-output probe.'}
(out/'analysis.json').write_text(json.dumps(report,indent=2).replace('/work/pi_nicholas_pizzo_uri_edu/jupiter-ffr-breeze/',''))
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','full_decode_passed':True,'frames':frames,'duration_seconds':duration,'compute_location':'Unity','assets':assets},indent=2))
print('DENSE_NATIVE_BUDGET_REVIEW_COMPLETE',flush=True)
