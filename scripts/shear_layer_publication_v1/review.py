"""Run only on Unity. Description of pilot contrasts, not a causal budget closure."""
import csv,json,sys,tomllib,subprocess
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio_ffmpeg
root,out=map(Path,sys.argv[1:]);out.mkdir(exist_ok=True,parents=True)
rows=list(csv.DictReader((out/'layer_statistics.csv').open()))
energy=[]
for cap in ('stable','neutral'):
 for shear in ('cyclonic','none','anticyclonic'):
  name=f'{cap}_{shear}';m=tomllib.loads((root/name/'metadata.toml').read_text())
  assert m['status']=='COMPLETE'
  for p in sorted((out/'layers'/name).glob('frame_*.jld2')):
   with h5py.File(p) as d:
    for i,layer in enumerate(('lower','interface','upper','full')):
     u,v=d['u'][i],d['v'][i];ny,nx=u.shape
     u=u-u.mean(axis=1,keepdims=True);v=v-v.mean(axis=1,keepdims=True)
     U=np.fft.fft2(np.concatenate((u,u[::-1])),norm='forward')
     V=np.fft.fft2(np.concatenate((v,-v[::-1])),norm='forward')
     E=.5*(abs(U)**2+abs(V)**2)
     ky=2*np.pi*np.fft.fftfreq(2*ny,m['ly']/ny)[:,None]
     kx=2*np.pi*np.fft.fftfreq(nx,m['lx']/nx)[None,:]
     K=np.hypot(kx,ky);large=(K>0)&(K<2*np.pi/(2*m['turbulence_wavelength']))&(kx!=0)
     assert np.isclose(E.sum(),.5*np.mean(u*u+v*v),rtol=1e-6,atol=1e-12)
     energy.append(dict(case=name,time_days=float(d['t'][()])/86400,layer=layer,
         nonzonal_ke_m2_s2=float(E.sum()),large_nonzonal_ke_m2_s2=float(E[large].sum()),
         resolved_high_k_energy_fraction=float(E[K>.6*min(np.pi*nx/m['lx'],np.pi*ny/m['ly'])].sum()/max(E.sum(),1e-30))))
with (out/'layer_energy.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(energy[0]));w.writeheader();w.writerows(energy)
fig,axs=plt.subplots(2,3,figsize=(15,8),layout='constrained',sharex=True)
for col,shear in enumerate(('cyclonic','none','anticyclonic')):
 for cap,style in [('stable','-'),('neutral','--')]:
  for layer,color in [('lower','tab:orange'),('upper','tab:blue'),('full','black')]:
   rr=[r for r in rows if r['case']==f'{cap}_{shear}' and r['layer']==layer]
   axs[0,col].plot([float(r['time_days']) for r in rr],[float(r['signed_skew_prime']) for r in rr],style,color=color,label=f'{cap}: {layer}')
   ee=[r for r in energy if r['case']==f'{cap}_{shear}' and r['layer']==layer]
   axs[1,col].plot([r['time_days'] for r in ee],[r['large_nonzonal_ke_m2_s2'] for r in ee],style,color=color)
 axs[0,col].axhline(0,color='.5',lw=.5);axs[0,col].set_title(f'{shear} initial shear')
 axs[1,col].set_xlabel('Earth days')
axs[0,0].set_ylabel('Rotation-aligned skewness of nonzonal vorticity')
axs[1,0].set_ylabel('Nonzonal KE at wavelengths >800 km [m²/s²]')
axs[0,2].legend(fontsize=8,ncol=2)
fig.suptitle('Shear × upper stratification: ten-day, single-seed screening\nNegative skewness is anticyclonic; skewness alone does not establish an LSV')
fig.savefig(out/'layer_contrasts.png',dpi=140);plt.close(fig)
movies=[]
for shear in ('cyclonic','none','anticyclonic'):
 movie=out/f'{shear}_stable_vs_neutral.mp4'
 subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),
    str(root/f'stable_{shear}'),str(root/f'neutral_{shear}'),'--labels','Stable upper layer','Neutral upper target',
    '--title',f'{shear.capitalize()} initial shear | dry f-plane pilot','--tracer','--output',str(movie)],check=True)
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-i',str(movie),'-f','null','-'],check=True)
 movies.append(movie.name)
(out/'review.json').write_text(json.dumps(dict(status='COMPLETE',movies=movies,
 caveats=['Single seed and ten days: not robust sign selection or a mature LSV claim.',
 'Initial shear evolves freely; compare actual zonal flow before interpreting null effects.',
 'Neutral-target control changes upper hydrostatic density and thermodynamic response as well as stability.',
 'Centered stretching/tilting are not a closed physical budget.',
 'No inverse-cascade inference from low-mode energy or morphology alone.'],
 statistics='Layer-wise nonzonal vorticity skewness; separately tracked total mean vorticity; absolute nonzonal KE above twice seed wavelength'),indent=2)+'\n')
print('SHEAR_LAYER_REVIEW_AND_MOVIES_COMPLETE',flush=True)
