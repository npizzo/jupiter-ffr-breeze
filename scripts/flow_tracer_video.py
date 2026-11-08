"""One synchronized movie: three depth columns, vorticity above passive concentration."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from vorticity_depth_videos import collect

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('case',type=Path);p.add_argument('--output',type=Path,required=True)
p.add_argument('--label',required=True);p.add_argument('--preview',action='store_true')
p.add_argument('--cloud',action='store_true',help='Add a third row of liquid+ice mass fraction')
p.add_argument('--fps',type=int,default=12);p.add_argument('--crf',type=int,default=28)
a=p.parse_args();m,files,times,vlimits=collect(a.case,a.preview)
if a.fps<1:raise ValueError('fps must be positive')
rows=[('vorticity_slices','Vertical vorticity ζ/f₀','RdBu_r'),('tracer_slices','Passive tracer concentration','viridis')]
if a.cloud:rows.append(('condensate_slices','Liquid + ice mass fraction','magma'))
tracer_range=[0.0,1.0];cloud_max=1e-12
for path in files:
 with h5py.File(path) as d:
  for key,_,_ in rows:
   field=d[key][:]
   if field.shape!=(3,m['ny'],m['nx']) or not np.isfinite(field).all():raise ValueError(f'Invalid {key}: {path}')
  c=d['tracer_slices'][:];tracer_range=[min(tracer_range[0],float(c.min())),max(tracer_range[1],float(c.max()))]
  if a.cloud:cloud_max=max(cloud_max,float(d['condensate_slices'][:].max()))
fig,axes=plt.subplots(len(rows),3,figsize=(15,4*len(rows)),layout='constrained',squeeze=False)
artists=[];scales=[]
for row,(key,label,cmap) in enumerate(rows):
 row_artists=[];row_scales=[]
 for col,z in enumerate(m['video_heights_actual_m']):
  if key=='vorticity_slices':limits=[-max(float(vlimits[col]),1e-12),max(float(vlimits[col]),1e-12)]
  elif key=='tracer_slices':limits=tracer_range
  else:limits=[0,cloud_max]
  with h5py.File(files[0]) as d:field=d[key][col]/m['f'] if key=='vorticity_slices' else d[key][col]
  ax=axes[row,col];im=ax.imshow(field,origin='lower',cmap=cmap,vmin=limits[0],vmax=limits[1],interpolation='nearest',extent=(0,m['lx']/1000,-m['ly']/2000,m['ly']/2000))
  ax.set(xlabel='x [km]',ylabel='y [km]',title=f'{label}\nz = {z/1000:.2f} km above bottom')
  fig.colorbar(im,ax=ax,shrink=.85);row_artists.append(im);row_scales.append(limits)
 artists.append(row_artists);scales.append(row_scales)
title=fig.suptitle('');a.output.parent.mkdir(parents=True,exist_ok=True)
writer=FFMpegWriter(fps=a.fps,codec='libx264',extra_args=['-crf',str(a.crf),'-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(a.output),dpi=100):
 for path,t in zip(files,times):
  with h5py.File(path) as d:
   for row,(key,_,_) in enumerate(rows):
    field=d[key][:]/m['f'] if key=='vorticity_slices' else d[key][:]
    for col in range(3):artists[row][col].set_data(field[col])
  clock=f'Earth day {t/86400:.3f}' if times[-1]>=86400 else f'{t/60:.1f} minutes'
  title.set_text(f'{"Partial run | " if a.preview else ""}{a.label} | {clock}')
  writer.grab_frame()
fig.savefig(a.output.with_suffix('.png'),dpi=100);plt.close(fig)
a.output.with_suffix('.json').write_text(json.dumps({'case':m['name'],'label':a.label,'layout':f'{len(rows)} rows, three depth columns','row_fields':[x[0] for x in rows],'height_m':m['video_heights_actual_m'],'times_seconds':times.tolist(),'fps':a.fps,'fixed_color_limits':scales,'temporal_interpolation':False,'preview':a.preview,'note':'Concentration is a passive mass fraction, not cloud brightness. Vorticity scales differ by depth; concentration scale is shared. No field is clipped.'},indent=2)+'\n')
print('Rendered synchronized flow/concentration:',a.output)
