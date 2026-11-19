"""Unity-only publication preparation; never infer a null result from failed fidelity checks."""
from pathlib import Path
import sys,json,tomllib,hashlib,subprocess,csv
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

out=Path(sys.argv[1]);roots=[Path(p) for p in sys.argv[2:]]
assert len(roots)==2
out.mkdir(parents=True,exist_ok=False)
plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
report={'scope':'Isolated warm plume, stable versus neutral interface; not LSV formation or convergence.',
        'cases':{},'native_budget':'Density-weighted layer mean curl; exact SSPRK stage integration including vertical transport. Physical stretching/tilting are separate reconstructions.',
        'fidelity_criteria':'Flag <8 cells in characteristic w/N; >10% w power at wavelength<4dx; >0.1% tracer in top 5 km or outer horizontal strips. These are screening criteria, not convergence proof.',
        'cautions':['Finite amplitude 2 K thermal from rest; finite heat pulse, no continuing bottom heat.',
                    'No sponge; gravity waves can encounter the rigid lid before tracer does.',
                    'Stable and neutral hydrostatic density profiles differ above the interface.',
                    'A pressure projection is nonlocal; absence of tracer at boundaries is insufficient to prove boundary independence.']}
data={};vmax=1e-9;cmax=1.
for root in roots:
    m=tomllib.loads((root/'metadata.toml').read_text());name=m['name']
    frames=sorted((root/'frames').glob('frame_*.jld2'));assert frames
    arrays=[]
    for path in frames:
        with h5py.File(path) as h:
            a={k:np.asarray(h[k]) for k in ['t','z','f','x','y','vorticity','tracer','w','divergence']}
        assert a['vorticity'].shape[0]==3 and np.isfinite(a['vorticity']).all()
        vmax=max(vmax,float(np.max(np.abs(a['vorticity']/a['f']))));cmax=max(cmax,float(np.max(a['tracer'])))
        arrays.append(a)
    data[name]=(root,m,arrays)
for name,(root,m,arrays) in data.items():
    with (root/'diagnostics.csv').open() as h:diag=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(h)]
    case={'source':str(root),'metadata':m,'diagnostics':diag,'full_field_screen':[],'movies':[],
          'conditional_interface_statistics':[],
          'conditional_note':'Eulerian w>RMS(w) and tracer C>0.001. Stretching reconstructed as -(f+zeta)*horizontal divergence; not substituted for native stage-integrated terms.'}
    for a in arrays:
        for k,zz in enumerate(a['z']):
            w=a['w'][k];mask=(w>np.sqrt(np.mean(w*w)))&(a['tracer'][k]>.001)
            if not mask.any():continue
            f=float(a['f']);zeta=a['vorticity'][k];div=a['divergence'][k]
            case['conditional_interface_statistics'].append({'time':float(a['t']),'height_m':float(zz),'area_fraction':float(mask.mean()),
                'mean_zeta_over_f':float(np.mean(zeta[mask])/f),'mean_divergence_over_f':float(np.mean(div[mask])/f),
                'mean_reconstructed_stretching_over_f2':float(np.mean(-(f+zeta[mask])*div[mask])/f**2)})
    # Stored fields are z,y,x in h5py; scalar 3D arrays are Float32 publication diagnostics.
    for path in sorted((root/'fields').glob('field_*.jld2')):
        with h5py.File(path) as h:
            z=h['z'][:];wface=h['w_face'][:].astype(float);C=h['tracer'][:].astype(float);t=float(h['t'][()]);x=h['x'][:];y=h['y'][:]
        w=.5*(wface[:-1]+wface[1:]);dx=float(x[1]-x[0]);dz=float(z[1]-z[0])
        kx=np.fft.fftfreq(len(x),dx);ky=np.fft.fftfreq(len(y),dx)
        near=np.hypot(ky[:,None],kx[None,:])>1/(4*dx)
        fractions=[]
        for zz in (75000,79500,80500,82000):
            k=int(np.argmin(abs(z-zz)));p=abs(np.fft.fft2(w[k]-w[k].mean()))**2
            fractions.append({'height_m':float(z[k]),'w_near_grid_fraction':float(p[near].sum()/p.sum()) if p.sum()>1e-25 else None})
        hit=z[np.max(C,axis=(1,2))>.01]
        cap=(z>=78000)&(z<=82000);wr=float(np.max(abs(w[cap])));N=np.sqrt(m['cap_n2'])
        case['full_field_screen'].append({'time':t,'height_of_tracer_C_gt_0_01_m':float(hit.max()) if hit.size else None,
            'tracer_penetration_cells':float(max(0,hit.max()-80000)/dz) if hit.size else 0,
            'characteristic_w_over_N_cells':wr/N/dz if N>0 else None,
            'top_5km_w_max_relative_to_global':float(np.max(abs(w[z>105000]))/max(np.max(abs(w)),1e-30)),
            'spectra':fractions})
    budgets=sorted((root/'budget').glob('budget_*.jld2'));assert budgets
    with h5py.File(budgets[-1]) as h:
        residual=h['residual'][:];zeta=h['z'][:];origin=h['origin'][:];parts=h['cumulative'][:]
        case['budget']={'checks':int(h['checks'][()]),'max_step_absolute_residual':float(h['max_step_absolute_residual'][()]),
                        'cumulative_max_absolute_residual':float(abs(residual).max()),'change_max_abs':float(abs(zeta-origin).max()),
                        'cumulative_terms_max_abs_by_layer_and_term':np.max(abs(parts),axis=(-1,-2)).tolist(),
                        'layers_m':[[70000,80000],[80000,84000],[84000,110000]],
                        'term_order':['horizontal_advection','vertical_advection','coriolis','viscosity']}
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    artists=[];extent=[float(arrays[0]['x'][0])/1000,float(arrays[0]['x'][-1])/1000,float(arrays[0]['y'][0])/1000,float(arrays[0]['y'][-1])/1000]
    for row in range(2):
        for col in range(3):
            norm=AsinhNorm(linear_width=.01,vmin=-vmax,vmax=vmax) if row==0 else AsinhNorm(linear_width=.01,vmin=0,vmax=cmax)
            im=axes[row,col].imshow(np.zeros_like(arrays[0]['tracer'][col]),origin='lower',extent=extent,cmap='RdBu_r' if row==0 else 'viridis',norm=norm,interpolation='nearest')
            axes[row,col].set(xlabel='x [km]',ylabel='y [km]',title=f"z = {float(arrays[0]['z'][col])/1000:.2f} km")
            fig.colorbar(im,ax=axes[row,col],label='ζ / f · asinh' if row==0 else 'Tracer C · asinh');artists.append(im)
    title=fig.suptitle('')
    video=out/f'{name}_plume.mp4'
    writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-crf','20','-pix_fmt','yuv420p','-movflags','+faststart'])
    with writer.saving(fig,str(video),100):
        for a in arrays:
            for col in range(3):artists[col].set_data(a['vorticity'][col]/a['f']);artists[3+col].set_data(a['tracer'][col])
            title.set_text(f"{name.capitalize()} upper layer · isolated +2 K plume · t = {float(a['t'])/60:.1f} min · {m['status']}")
            writer.grab_frame()
    fig.savefig(out/f'{name}_plume.png',dpi=120);plt.close(fig)
    # Decode every frame on Unity; a successful encode alone is insufficient.
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-i',str(video),'-f','null','-'],check=True)
    case['movies']=[{'file':video.name,'frame_count':len(arrays),'fps':12,'decode':'PASS','sha256':hashlib.sha256(video.read_bytes()).hexdigest()}]
    case['movie_scales']={'vorticity':'zeta/f, fixed across both cases, depths and time','vorticity_asinh_width':.01,'vorticity_limits':[-vmax,vmax],'tracer_asinh_width':.01,'tracer_limits':[0,cmax]}
    report['cases'][name]=case
(out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print('PLUME_REVIEW_AND_MOVIES_COMPLETE',out,flush=True)
