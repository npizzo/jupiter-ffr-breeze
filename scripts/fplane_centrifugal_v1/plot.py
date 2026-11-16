from pathlib import Path
import json,sys,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(sys.argv[1]);r=json.loads((out/'analysis.json').read_text());h=np.load(out/'height_statistics.npz')
plt.rcParams.update({'font.size':11})
z=h['z']/1000;bins=h['bins'];bc=(bins[1:]+bins[:-1])/2;pdf=h['pdf'];m=h['moments'];face=h['native_face']
names=['Bottom 0–10 km','Lower 10–60 km','Upper convection 60–80 km','Cap 80–120 km']
masks=[z<10,(z>=10)&(z<60),(z>=60)&(z<80),z>=80]
fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
for name,sel in zip(names,masks):
    pp=pdf[sel].mean(0);axs[0,0].semilogy(bc,np.where(pp>0,pp,np.nan),label=name)
axs[0,0].axvline(-1,c='k',ls='--',lw=1);axs[0,0].axvline(1,c='.6',ls=':',lw=1)
axs[0,0].set(xlim=(-2,2),xlabel=r'$\zeta/f$',ylabel='Probability density',title='Vorticity distributions by height band')
axs[0,0].legend(fontsize=9)
axs[0,1].plot(h['raw_skewness'],z,label='Raw skewness')
axs[0,1].plot(h['central_skewness'],z,'--',label='Central skewness')
axs[0,1].axvline(0,c='.6',lw=.7);axs[0,1].set(xlabel='Pooled skewness',ylabel='Height (km)',title='Height dependence');axs[0,1].legend(fontsize=9)
axs[1,0].plot(m[:,0,5,:].mean(0)*100,z,label=r'Centered: $\zeta/f<-1$')
axs[1,0].plot(m[:,0,6,:].mean(0)*100,z,label=r'Centered: $\zeta/f>1$')
axs[1,0].plot(face[:,2,:].mean(0)*100,z,'--',label=r'Native corners: $\zeta/f<-1$')
axs[1,0].plot(face[:,3,:].mean(0)*100,z,'--',label=r'Native corners: $\zeta/f>1$')
axs[1,0].set(xlabel='Area × snapshot fraction (%)',ylabel='Height (km)',title='Threshold exceedances and operator sensitivity');axs[1,0].legend(fontsize=8)
n2=np.array([x['mean_N2_profile'] for x in r['extremes_by_day']])
axs[1,1].plot(n2.mean(0),z);axs[1,1].axvline(0,c='.6',lw=.7)
axs[1,1].set(xlabel=r'Mean-profile $N^2$ (s$^{-2}$)',ylabel='Height (km)',title='Convecting region and maintained stable cap')
for ax in [axs[0,1],axs[1,0],axs[1,1]]:ax.axhline(80,c='.5',ls=':',lw=1);ax.set_ylim(0,120)
fig.suptitle('Mature heated f-plane LSV · days 90–100\n11 daily 3D states; central 80% of channel width; cell-centered curl',fontsize=15)
fig.savefig(out/'height_distributions.png',dpi=150);plt.close(fig)

selected=json.loads((out/'selected_vortex_profiles.json').read_text())[:3]
fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained')
for j,c in enumerate(selected):
    p=c['profile'];rr=np.array(p['r_m'])/1000
    axs[0,j].plot(rr,p['azimuthal_zeta_over_f'],label='Mean sampled curl')
    axs[0,j].plot(rr,p['circulation_zeta_over_f'],'--',label='Circulation derivative')
    axs[0,j].axhline(-1,c='k',ls=':');axs[0,j].set(ylabel=r'Ring-mean $\zeta/f$',xlabel='Radius (km)',
        title=f"Day {c['day']:g}, z={c['z_m']/1000:.2f} km\nCore ζ/f={c['core_zeta_over_f']:.2f}; coherence={c['coherence']:.2f}")
    axs[1,j].plot(rr,p['Phi_over_f2'],label='Using sampled curl')
    axs[1,j].plot(rr,p['Phi_circulation_over_f2'],'--',label='Using circulation derivative')
    axs[1,j].axhline(0,c='k',lw=.8);axs[1,j].set(xlabel='Radius (km)',ylabel=r'$\Phi/f^2$',
        title='Conditional screen: '+('passes' if c['axisymmetric_screen_pass'] else 'does not pass'))
    for ax in axs[:,j]:ax.axvspan(0,4*r['dx_m']/1000,color='.9',zorder=-1);ax.legend(fontsize=8)
fig.suptitle('Strong anticyclonic candidates: an axisymmetric stability screen\nGray region: fewer than four radial grid cells; profiles are diagnostic, not stability eigenmodes',fontsize=14)
fig.savefig(out/'vortex_profiles.png',dpi=150);plt.close(fig)

energy=json.loads((out/'energy.json').read_text());track=json.loads((out/'tracking.json').read_text());tt=track['timeline']
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
days=[e['day'] for e in energy]
axs[0,0].plot(days,[e['deviation_ke_m2_s2'] for e in energy],'o-');axs[0,0].set(xlabel='Earth day',ylabel=r'm$^2$ s$^{-2}$',title='Depth-dependent kinetic energy (daily states)')
axs[0,1].plot(days,[e['barotropic_to_deviations_m2_s3'] for e in energy],'o-',label='Depth-mean shear → deviations')
axs[0,1].plot(days,[e['buoyancy_production_proxy_m2_s3'] for e in energy],'o-',label='Thermal covariance proxy')
axs[0,1].axhline(0,c='.6',lw=.7);axs[0,1].set(xlabel='Earth day',ylabel=r'm$^2$ s$^{-3}$',title='Instantaneous energy-production diagnostics');axs[0,1].legend(fontsize=8)
axs[1,0].plot([t['day'] for t in tt],[t['q_min'] for t in tt],label='Most anticyclonic')
axs[1,0].plot([t['day'] for t in tt],[t['q_max'] for t in tt],label='Most cyclonic')
axs[1,0].axhline(-1,c='k',ls=':');axs[1,0].axhline(1,c='.6',ls=':');axs[1,0].set(xlabel='Earth day',ylabel=r'$\zeta/f$',title='Three-hour extrema at 19.69 km');axs[1,0].legend(fontsize=8)
longest=sorted([v for v in track['tracks'].values() if v[0]['sign']<0],key=len,reverse=True)[:5]
for seq in longest:axs[1,1].plot([s['day'] for s in seq],[s['q_over_f'] for s in seq],'.-',label=f"Track {seq[0]['id']}")
axs[1,1].axhline(-1,c='k',ls=':');axs[1,1].set(xlabel='Earth day',ylabel=r'Core $\zeta/f$',title='Five longest tentative anticyclone tracks');axs[1,1].legend(fontsize=8)
fig.suptitle('Temporal evidence and its limits\nDaily 3D sampling cannot resolve growth on the 0.96-hour inverse-f time; energy curves are not a closed budget',fontsize=14)
fig.savefig(out/'temporal_evidence.png',dpi=150);plt.close(fig)

with np.load(out/'lsv_maps.npz') as m:Z=m['Z_over_f'];x=m['x']/1000;y=m['y']/1000
lsv=json.loads((out/'lsv_profiles.json').read_text());fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
lim=max(abs(Z.min()),abs(Z.max()));im=axs[0].imshow(Z.T,origin='lower',extent=[0,x[-1]+r['dx_m']/2000,y[0]-r['dy_m']/2000,y[-1]+r['dy_m']/2000],cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='equal')
for c in lsv:axs[0].plot(c['x_m']/1000,c['y_m']/1000,'ko' if c['sign']>0 else 'ks',ms=6,fillstyle='none')
axs[0].set(xlabel='x (km)',ylabel='y (km)',title='Day-100 density-weighted depth-mean vorticity');fig.colorbar(im,ax=axs[0],label='Z/f',shrink=.8)
for c in lsv:axs[1].plot(np.array(c['profile']['r_m'])/1000,c['profile']['Phi_over_f2'],label='Cyclonic maximum' if c['sign']>0 else 'Anticyclonic minimum')
axs[1].axhline(0,c='k',lw=.8);axs[1].set(xlabel='Radius (km)',ylabel=r'$\Phi/f^2$',title='Local profiles around the depth-mean extrema');axs[1].legend()
fig.savefig(out/'lsv_context.png',dpi=150);plt.close(fig)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.glob('*.png')}
(out/'figure_manifest.json').write_text(json.dumps({'status':'PASS','compute_location':'Unity','figures':assets},indent=2)+'\n')
print('PLOTS_COMPLETE',assets,flush=True)
