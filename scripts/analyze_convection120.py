"""Compare five-day convection, explicit day-2/day-5 windows and resolution."""
import argparse,csv,gzip,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze_controlled import read_case,metrics
from analyze_neutral import snapshot,plot_sections
from analyze_turbulence import velocity_map,spectrum
from analyze_convection import save_figure

KEYS=('kinetic_energy_J_m2','w_rms_center_m_s','vorticity_rms_over_f','buoyancy_power_center_W_m2')
def window(c,start,end):
 t=c['t'];keep=(t>=start*3600)&(t<=end*3600);ts=t[keep]
 assert ts[0]==start*3600 and ts[-1]==end*3600
 return {k:float(np.trapezoid(c['d'][k][keep],ts)/(ts[-1]-ts[0])) for k in KEYS}
def changes(a,b):
 # Near-zero signed buoyancy work in controls makes percentage ratios misleading.
 return {k:100*(a[k]/b[k]-1) for k in KEYS[:3]}
def main():
 p=argparse.ArgumentParser();p.add_argument('campaign',type=Path);p.add_argument('--previous',type=Path,required=True);p.add_argument('--output',type=Path,default=Path('gallery/convection120'));args=p.parse_args()
 root,dest=args.campaign,args.output;dest.mkdir(parents=True,exist_ok=True)
 summary={'duration_hours':120,'windows_hours':{'day2':[24,48],'day5':[96,120]},'cases':{},'pairs':{}}
 groups=[]
 for nz,color in [(64,'#0072B2'),(128,'#D55E00')]:
  cases=[]
  for mode,style in [('heated','-'),('cooled','--')]:
   name=f'convection120_f030_{mode}_h256_z{nz}';folder=root/name
   c=read_case(folder,f'{nz} levels: '+('heated + cooled' if mode=='heated' else 'cooling only'),color,style)
   assert c['m']['stop_seconds']==432000 and c['m']['nx']==256 and c['m']['ny']==256
   assert c['m']['beta']==0 and c['m']['jet_speed']==0 and c['m']['initial_eddy_rossby']==.1
   r={k:v for k,v in metrics(c).items() if not k.startswith('late_mean_')}
   r['day2']=window(c,24,48);r['day5']=window(c,96,120);r['day5_vs_day2_percent']=changes(r['day5'],r['day2'])
   c['initial_map']=velocity_map(folder/'frame_00000.jld2');c['final_map']=velocity_map(folder/'frame_00120.jld2')
   k,e=spectrum(c['final_map']);r['final_slice_low_mode_1_2_energy_fraction']=float(e[1:3].sum()/e[1:].sum())
   r['final_slice_energy_weighted_length_m']=float(c['m']['lx']/(2*np.pi)*np.sum(e[1:]/k[1:])/e[1:].sum())
   summary['cases'][name]=r;cases.append(c)
   with gzip.open(dest/(name+'_diagnostics.csv.gz'),'wt',newline='') as f:
    writer=csv.writer(f);writer.writerow(c['d'].dtype.names);writer.writerows([list(row) for row in c['d']])
  a,b=cases
  assert all(np.array_equal(a['initial_map'][k],b['initial_map'][k]) for k in a['initial_map'])
  assert np.array_equal(a['initial_theta'],b['initial_theta'])
  match=('nx','ny','nz','lx','ly','height','cap_height','cap_n2','convective_n2','f','beta','viscosity','thermal_diffusivity','cap_diffusivity','diffusivity_taper_width','upper_cooling_flux','perturbation','seed','turbulence_seed','turbulence_rossby','turbulence_wavelength','turbulence_vertical_ratio','script_sha256','initial_turbulence_script_sha256')
  assert all(a['m'][k]==b['m'][k] for k in match)
  summary['pairs'][f'z{nz}']={w+'_heated_vs_control_percent':changes(summary['cases'][a['name']][w],summary['cases'][b['name']][w]) for w in ('day2','day5')}
  summary['pairs'][f'z{nz}']['day5_buoyancy_work_increment_W_m2']=summary['cases'][a['name']]['day5']['buoyancy_power_center_W_m2']-summary['cases'][b['name']]['day5']['buoyancy_power_center_W_m2']
  groups.append(cases)
  frames=[snapshot(root/c['name']/'frame_00120.jld2') for c in cases]
  plot_sections(cases,frames,dest,f'z{nz}_final_sections',[c['label'] for c in cases],f'30× heating: 120 hours, {nz} levels')
  hours=[24,48,72,120];frames=[snapshot(root/a['name']/f'frame_{h:05d}.jld2') for h in hours]
  assert all(frame['t']==h*3600 for frame,h in zip(frames,hours))
  plot_sections([a]*4,frames,dest,f'z{nz}_heated_snapshots',[f'{h} h' for h in hours],f'30× heating: {nz} levels')
  fig,axes=plt.subplots(1,4,figsize=(14,4),layout='constrained')
  limit=max(float(np.max(abs(frame['zeta']))/a['m']['f']) for frame in frames)
  for ax,frame,h in zip(axes,frames,hours):
   im=ax.imshow(frame['zeta']/a['m']['f'],origin='lower',cmap='RdBu_r',vmin=-limit,vmax=limit,extent=(0,1200,0,1200))
   ax.set(title=f'{h} h',xlabel='x [km]',ylabel='y [km]')
  fig.colorbar(im,ax=list(axes),shrink=.8,label='ζ/f near 20 km (fixed scale)');fig.suptitle(f'Heated run, {nz} levels: horizontal vorticity')
  save_figure(fig,dest,f'z{nz}_vorticity_snapshots')
 summary['vertical_refinement_day5_percent']=changes(summary['cases'][groups[1][0]['name']]['day5'],summary['cases'][groups[0][0]['name']]['day5'])
 old=json.loads(args.previous.read_text())['cases']['strong48_f030_heated']
 mapping={'kinetic_energy_J_m2':'late_mean_ke_J_m2','w_rms_center_m_s':'late_mean_w_rms_m_s','vorticity_rms_over_f':'late_mean_vorticity_rms_over_f','buoyancy_power_center_W_m2':'late_mean_buoyancy_power_W_m2'}
 summary['horizontal_refinement_day2_percent']=changes(summary['cases'][groups[0][0]['name']]['day2'],{k:old[v] for k,v in mapping.items()})
 summary['note']='Same seed, grid-dependent sampling; two grids do not establish convergence. Slice spectra do not measure spectral transfer.'
 (dest/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
 for ax,key,title in zip(axes.flat,KEYS,['Column KE [J/m²]','Vertical RMS speed [m/s]','Volume RMS ζ/f','Buoyancy work [W/m²]']):
  for cases in groups:
   for c in cases:ax.plot(c['t']/3600,c['d'][key],c['style'],color=c['color'],label=c['label'])
  ax.set(xlabel='Time [Earth hours]',ylabel=title);ax.grid(alpha=.2)
  if key in KEYS[:3]:ax.set_yscale('log')
 axes[0,0].legend(fontsize=8);fig.suptitle('Five-day 30× convection: heated and cooling-only pairs')
 save_figure(fig,dest,'evolution')
 fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
 for cases in groups:
  for c in cases:
   k,e=spectrum(c['final_map']);axes[0].loglog(k[1:],e[1:],c['style'],color=c['color'],label=c['label'])
   axes[1].plot(c['fields']['heat_flux_profile'],c['z']/1000,c['style'],color=c['color'],label=c['label'])
 axes[0].set(xlabel='Horizontal shell index',ylabel='Slice horizontal KE [m²/s²]');axes[0].legend(fontsize=7)
 axes[1].set(xlabel='Resolved upward heat flux [W/m²]',ylabel='Height [km]');axes[1].axhline(80,color='gray',lw=.5)
 fig.suptitle('120 hours: near-20-km spectra and vertical heat transport');save_figure(fig,dest,'spectra_heat_flux')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
