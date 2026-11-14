from pathlib import Path
import json,sys
import numpy as np
from scipy.special import expit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(sys.argv[1]);r=json.loads((out/'diagnostics.json').read_text());rows=r['snapshots'];t=[x['time_days'] for x in rows]
fig,ax=plt.subplots(2,3,figsize=(16,9),constrained_layout=True)
for col,kind in enumerate(['full','interior','nonzonal_interior']):
 with np.load(out/(kind+'_scatter.npz')) as d:
  i=np.linspace(0,len(d['Z'])-1,min(30000,len(d['Z']))).astype(int)
  im=ax[0,col].scatter(d['Psi'][i]/1e8,d['Z'][i]/abs(rows[0]['f']),c=d['day'][i],s=.7,rasterized=True,cmap='viridis')
  ax[0,col].set(xlabel='Ψ / 10⁸ m² s⁻¹',ylabel='Z / |f|',title=kind.replace('_',' '))
  fig.colorbar(im,ax=ax[0,col],label='Earth days')
  fit=r['pooled_fits'][kind]
  if fit['status']=='FIT':
   px=np.linspace(np.min(d['Psi']),np.max(d['Psi']),300);pq=(px-fit['Psi_mean'])/fit['Psi_standard_deviation']
   qz=fit['Z_mean']+fit['Z_standard_deviation']*np.polynomial.polynomial.polyval(pq,fit['quadratic_standardized_coefficients'])
   lz=fit['logistic_sigma2']+(fit['logistic_sigma1']-fit['logistic_sigma2'])*expit(fit['logistic_standardized_intercept']+fit['logistic_standardized_slope']*pq)
   ax[0,col].plot(px/1e8,qz/abs(rows[0]['f']),color='black',lw=1,label='quadratic')
   ax[0,col].plot(px/1e8,lz/abs(rows[0]['f']),color='darkred',lw=1,label='logistic')
   ax[0,col].legend(loc='lower right',fontsize=8)
   ax[0,col].text(.03,.97,f"Quadratic R²={fit['quadratic_R2']:.3f}\nLogistic R²={fit['logistic_R2']:.3f}",transform=ax[0,col].transAxes,va='top',bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
for scale in ['0','50','100','200']:
 ax[1,0].plot(t,[x['moments']['interior']['full_Z'][scale]['raw_skewness'] for x in rows],label=f'{scale} km FWHM')
ax[1,0].set(xlabel='Earth days',ylabel='Density-depth-mean Z skewness',title='Interior; 50 km unresolved on 256²');ax[1,0].legend()
ax[1,1].plot(t,[x['moments']['interior']['FuSun_reference']['rotation_signed_skewness_volume'] for x in rows],label='sign(f) Sζ, 3D volume')
ax[1,1].plot(t,[x['moments']['interior']['FuSun_reference']['2p5_Ro_volume'] for x in rows],label='2.5 Ro, reference only')
ax[1,1].set(xlabel='Earth days',ylabel='3D vorticity statistic',title='Fu–Sun comparison: different physical regime');ax[1,1].legend()
for k in ['V','H','barotropic_advection','mean_divergence']:
 ax[1,2].plot(t,[x['source']['interior']['third_moment_source_s_inv4'][k] for x in rows],label=k)
ax[1,2].set(xlabel='Earth days',ylabel='3〈Z² source〉 [s⁻⁴]',title='Centered source reconstruction; remainder retained');ax[1,2].legend(fontsize=8)
fig.suptitle(r['metadata']['name']+' | hypotheses tested, not assumed');fig.savefig(out/'diagnostic_summary.png',dpi=150);plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
for scale in ['0','50','100','200']:
 for col,kind in enumerate(['full_Z','nonzonal_Z']):
  p=rows[-1]['moments']['interior'][kind][scale]['PDF_zeta_over_abs_f'];x=np.array(p['bin_edges']);pdf=np.array(p['density'])
  ax[col].plot((x[1:]+x[:-1])/2,pdf,label=f'{scale} km FWHM')
for a,title in zip(ax,['Full interior Z','Z with zonal mean removed']):
 a.set(xlabel='Z / |f|',ylabel='Probability density',title=title);a.legend()
fig.suptitle(r['metadata']['name']+f" | final PDF at day {t[-1]:g}");fig.savefig(out/'vorticity_PDF.png',dpi=150);plt.close(fig)
