from pathlib import Path
import sys,json,csv,tomllib,hashlib,subprocess,shutil
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter

root=Path.cwd();scratch=Path('/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127')
case=sys.argv[1];assert case in ('marcus','anelastic')
if case=='marcus':
    sources=[scratch/'marcus-CV1-fourday-2345592/parent'];name='marcus-cv1-four-days';end=345600.;fps=12;repeat=1
    title='Marcus CV1-like reconstruction · initialized anticyclone + shear'
else:
    chain=json.loads((root/'runs/cai-anelastic-continue-v1/chain-state.json').read_text())
    sources=[scratch/'anelastic-B7-2342242/parent']+[Path(x['output_directory'])/'parent' for x in chain['jobs'] if 'outcome' in x]
    name='cai-anelastic-b7-through240';end=240.;fps=10;repeat=4
    title='Cai B7 anelastic counterpart · Ω=2.56 · continued from rest'
out=root/'rendered'/name;out.mkdir(parents=True,exist_ok=False)
frames={}
for src in sources:
    for p in src.glob('frame_*.jld2'):
        with h5py.File(p) as h:t=float(h['t'][()])
        if t<=end+1e-8:frames[round(t,7)]=p
times=np.array(sorted(frames));assert times[0]==0 and times[-1]==end
if case=='marcus':assert len(times)==161
else:assert np.array_equal(times,np.arange(0,241,10))
Q=[];C=[];stats=[]
for t in times:
    with h5py.File(frames[t]) as h:
        f=float(h['f'][()]);q=h['vorticity_slices'][:]/f;c=h['tracer_slices'][:];heights=h['video_heights'][:]
        nx=q.shape[-1];ny=q.shape[-2];lx=float(h['lx'][()]);ly=float(h['ly'][()])
    if case=='marcus':
        x=-lx/2+(np.arange(nx)+.5)*lx/nx;y=-ly/2+(np.arange(ny)+.5)*ly/ny
        record={'time_seconds':float(t),'earth_days':float(t/86400)}
        for j,z in enumerate(heights):
            mask=(x[None,:]/(1.5*4.6e6))**2+((y[:,None]-z/np.tan(np.deg2rad(-23)))/4.6e6)**2<=1
            record[f'raw_fixed_core_circulation_z{j}']=float((q[j]*f)[mask].sum()*lx/nx*ly/ny)
            record[f'excess_over_background_circulation_z{j}']=float((q[j]*f-9.1e-6)[mask].sum()*lx/nx*ly/ny)
        stats.append(record)
        ix=np.where(abs(x)<=18e6)[0];iy=np.where(abs(y)<=12e6)[0]
        q=q[:,iy][:,:,ix];c=c[:,iy][:,:,ix];extent=[(x[ix[0]]-.5*lx/nx)/1000,(x[ix[-1]]+.5*lx/nx)/1000,(y[iy[0]]-.5*ly/ny)/1000,(y[iy[-1]]+.5*ly/ny)/1000]
    else:extent=[0,lx,0,ly]
    Q.append(q);C.append(c)
Q=np.asarray(Q);C=np.asarray(C);assert np.isfinite(Q).all() and np.isfinite(C).all()
native=[]
for n,src in enumerate(sources):
    shutil.copy2(src/'diagnostics.csv',out/f'diagnostics_segment{n:02d}.csv')
    with (src/'diagnostics.csv').open() as fp:native.extend({k:float(v) for k,v in r.items()} for r in csv.DictReader(fp) if float(r['time'])<=end+1e-8)
meta=tomllib.loads((sources[-1]/'metadata.toml').read_text())
summary={'case':case,'status':'COMPLETE' if case=='marcus' else 'PARTIAL_CONTINUING','movie_time_end':end,'source_segments':[s.parent.name for s in sources],
         'grid':[384,384,256] if case=='marcus' else [512,512,101], 'saved_times':times.tolist(),'native_final':native[-1],
         'max_recorded_mach':max(r['max_mach'] for r in native),'notes':[]}
if case=='marcus':
    summary.update(target_earth_days=4,elapsed_gpu_job='08:27:00',initial=stats[0],final=stats[-1],
        centre_raw_circulation_retention=stats[-1]['raw_fixed_core_circulation_z1']/stats[0]['raw_fixed_core_circulation_z1'],
        centre_excess_circulation_retention=stats[-1]['excess_over_background_circulation_z1']/stats[0]['excess_over_background_circulation_z1'],
        max_mass_flux_constraint=max(r['mass_flux_divergence_relative'] for r in native),max_relative_theta_departure=max(r['max_relative_theta_departure'] for r in native))
    summary['notes']=['Four Earth days test short-term persistence of an initialized anticyclone, not spontaneous formation or sign selection.',
      'Reconstructed reference atmosphere and numerical/boundary departures are documented in the source README; not an exact-author-input reproduction.',
      'Movie shows a fixed central crop x±18000km,y±12000km from the105000x70000km domain. No tracking or time-varying crop.',
      'Southern hemisphere: positive raw zeta is anticyclonic; negative zeta/f is anticyclonic.',
      'Circulation uses a fixed ellipse, includes background shear, and can change through deformation/advection. Excess circulation subtracts the prescribed9.1e-6/s background, not an independently evolved shear field.',
      'Global kinetic energy is dominated by the background flow; its constancy is not evidence of vortex-energy conservation.']
    with (out/'core_circulation.csv').open('w') as fp:w=csv.DictWriter(fp,fieldnames=list(stats[0]));w.writeheader();w.writerows(stats)
else:
    summary.update(checkpoint_time=248.44410780816244,target_time=5000)
    summary['notes']=['Time is nondimensional. Partial evolving convection; no mature vortex-sign conclusion.',
      'The continuation restores physical prognostic state, but v1 restart diagnostics sampled derived temperature before refresh. Its segment profile-drift and energy-proxy baselines are invalid; do not interpret those columns as budgets. The v2 continuation corrects this.',
      'Mach and fixed-reference potential-temperature departures are measured after integration state updates. A low mass-flux residual does not establish physical validity of the anelastic approximation.',
      'Each frame is repeated four times for readable playback; no interpolation of fields.']
shutil.copy2(sources[-1]/'metadata.toml',out/'metadata.toml')
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
qmax=np.maximum(abs(Q).max(axis=(0,2,3)),1e-8);means=C.mean(axis=(-2,-1),dtype=float);stds=C.std(axis=(-2,-1),dtype=float)
center=.5 if case=='marcus' else .1
width=.01 if case=='marcus' else float(max(1e-7,np.min(np.median(stds[-5:],axis=0))/3))
cmin=float(C.min());cmax=float(C.max())
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':10})
fig,axs=plt.subplots(2,3,figsize=(16,9),dpi=100);fig.subplots_adjust(left=.055,right=.955,bottom=.065,top=.86,wspace=.43,hspace=.35)
ims=[];labels=[]
for r in range(2):
    for j in range(3):
        norm=AsinhNorm(linear_width=max(qmax[j]/30,1e-4),vmin=-qmax[j],vmax=qmax[j]) if r==0 else FuncNorm((lambda x:np.arcsinh((x-center)/width),lambda x:center+width*np.sinh(x)),vmin=cmin,vmax=cmax)
        im=axs[r,j].imshow(np.zeros(Q.shape[-2:]),origin='lower',extent=extent,interpolation='nearest',norm=norm,cmap='RdBu_r' if r==0 else 'BrBG')
        ims.append(im);labels.append(axs[r,j].set_title(''));axs[r,j].set(xlabel='x [km]' if case=='marcus' else 'x/H',ylabel='y [km]' if case=='marcus' else 'y/H')
        cb=fig.colorbar(im,ax=axs[r,j],fraction=.045,pad=.03,label='ζ/f₀; negative = anticyclonic' if r==0 else 'Absolute tracer C')
        if r==1:
            tick=center+width*np.sinh(np.linspace(np.arcsinh((cmin-center)/width),np.arcsinh((cmax-center)/width),5));cb.set_ticks(tick,labels=[f'{x:.6g}' for x in tick])
heading=fig.suptitle('',fontsize=13,y=.98)
def update(i):
    tm=f'{times[i]/86400:.3f} Earth days' if case=='marcus' else f't={times[i]:g} (nondimensional)'
    heading.set_text(title+' · '+tm+'\n'+('Fixed central crop; initialized-vortex persistence, no convection' if case=='marcus' else 'Partial continuation; checkpoint t=248.44, target5000')+'\nFixed color mappings; vorticity ranges differ by depth; absolute tracer, no framewise rescaling')
    for j in range(3):
        height=f'z={heights[j]/1000:.2f} km relative to1bar' if case=='marcus' else f'z/H={heights[j]:.4f}'
        ims[j].set_data(Q[i,j]);labels[j].set_text(height)
        ims[j+3].set_data(C[i,j]);labels[j+3].set_text(height+f'\nmean C={means[i,j]:.6f}; RMS={stds[i,j]:.2g}')
update(len(times)-1);fig.savefig(out/'final.jpg',dpi=100,pil_kwargs={'quality':90})
fn=out/'vorticity_tracer.mp4';writer=FFMpegWriter(fps=fps,codec='libx264',bitrate=2000,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(fn),dpi=100):
    for i in range(len(times)):
        update(i)
        for _ in range(repeat):writer.grab_frame()
plt.close(fig)
nf,duration=imageio_ffmpeg.count_frames_and_secs(str(fn));assert nf==len(times)*repeat
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(fn),'-f','null','-'],check=True)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'media.json').write_text(json.dumps({'status':'MEDIA_VERIFIED','full_decode_pass':True,'frames':nf,'duration':duration,'heights':heights.tolist(),'fixed_extent':extent,'vorticity_limits_per_depth':qmax.tolist(),'tracer':{'center':center,'asinh_width':width,'min':cmin,'max':cmax},'assets':assets},indent=2)+'\n')
print('MEDIA_VERIFIED',name,flush=True)
