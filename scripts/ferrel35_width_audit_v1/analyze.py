"""Read-only matched-width audit. Execute only on Unity; no simulation changes."""
from pathlib import Path
import hashlib, json, tomllib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path('rendered/ferrel35-width-audit-v1')
OUT.mkdir(parents=True, exist_ok=False)
WINDOWS = [(50, 75), (75, 100), (50, 100)]
ROOTS = {'narrow': Path('data/cascade35-v1/cascade35_heated'),
         'wide': Path('data/ferrel35-width-v4/wide')}
FIELDS = ['zonal_u', 'zonal_v', 'zonal_w', 'zonal_theta', 'eddy_uv', 'eddy_uw', 'mass_streamfunction']
QA = ['mass_divergence_max', 'mass_flux_reconstruction_max', 'top_streamfunction_max']
NATIVE_NAMES = ['advection', 'coriolis', 'diffusion', 'external_forcing']
PHYSICAL = ['nx','nz','lx','height','f','beta','latitude','heat_flux','gravity','theta0',
            'molar_mass','heat_capacity','base_pressure','cap_height','cap_n2','cap_core_bottom',
            'cap_core_top','cap_min_fraction','sponge_time','cooling_time','viscosity',
            'thermal_diffusivity','tracer_diffusivity','cap_diffusivity','mixing_model',
            'balance_heat','cooling_model','output_seconds','initial_flow_family',
            'turbulence_wavelength','turbulence_rossby','turbulence_seed']

def derivative(a, h, axis):
    # Matches the production observer: centered interior, first-order one-sided edge.
    return np.gradient(a, h, axis=axis, edge_order=1)

def smooth(a, dz, dy):
    # Display only. All main metrics use unsmoothed data.
    for axis, width in [(0, 5000/dz), (1, 200000/dy)]:
        rad=int(np.ceil(4*width)); x=np.arange(-rad,rad+1)
        kernel=np.exp(-.5*(x/width)**2);kernel/=kernel.sum()
        a=np.apply_along_axis(lambda v:np.convolve(np.pad(v,rad,mode='symmetric'),kernel,mode='valid'),axis,a)
    return a

def rms(a, mask, rho):
    weights=np.broadcast_to(rho[:,None],mask.shape)[mask]
    return float(np.sqrt(np.sum(weights*a[mask]**2)/weights.sum()))

def cosine(a,b,mask,rho):
    weights=np.broadcast_to(rho[:,None],mask.shape)[mask];weights/=weights.sum()
    return float(np.sum(weights*a[mask]*b[mask])/max(rms(a,mask,rho)*rms(b,mask,rho),1e-30))

def psi_center(a):
    return (a[:-1,:-1]+a[1:,:-1]+a[:-1,1:]+a[1:,1:])/4

def load_case(label, root):
    files={}; metas=[]
    for seg in range(5,11):
        p=root/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text())
        assert m['status']=='COMPLETE'
        metas.append({'segment':seg,'metadata_sha256':hashlib.sha256((p/'metadata.toml').read_bytes()).hexdigest()})
        for f in p.glob('frame_*.jld2'):
            i=int(f.stem.split('_')[1])
            if 400<=i<=800:
                assert i not in files;files[i]=f
    assert sorted(files)==list(range(400,801)), 'Missing three-hour snapshots'
    with h5py.File(files[400]) as d:
        rho=d['density'][:];y=d['y'][:];z=d['z'][:]
    dy=m['ly']/m['ny'];dz=m['height']/m['nz']
    assert np.allclose(np.diff(y),dy) and np.allclose(np.diff(z),dz)
    averages=[{} for _ in WINDOWS]; states={}; qa={k:0. for k in QA}; digest=hashlib.sha256()
    for i in range(400,801):
        with h5py.File(files[i]) as d:
            assert float(d['t'][()])==i*10800
            assert np.array_equal(d['density'][:],rho)
            a={k:d[k][:] for k in FIELDS}
            for k in QA:qa[k]=max(qa[k],float(d[k][()]))
            u,v,w=a['zonal_u'],a['zonal_v'],a['zonal_w']
            a['mean_y']=-derivative(rho[:,None]*u*v,dy,1)/rho[:,None]
            a['mean_z']=-derivative(rho[:,None]*u*w,dz,0)/rho[:,None]
            a['eddy_y']=-derivative(rho[:,None]*a['eddy_uv'],dy,1)/rho[:,None]
            a['eddy_z']=-derivative(rho[:,None]*a['eddy_uw'],dz,0)/rho[:,None]
            a['fv']=(m['f']+m['beta']*y)[None,:]*v
            a['v_square']=v*v;a['w_square']=w*w
            a['eddy_ke_conversion']=rho[:,None]*(a['eddy_uv']*derivative(u,dy,1)+a['eddy_uw']*derivative(u,dz,0))
            if label=='wide':
                a['native_snapshot']=d['native_rhou_rhs'][:]/rho[None,:,None]
                a['reconstruction']=a['native_snapshot'][0]-a['mean_y']-a['mean_z']-a['eddy_y']-a['eddy_z']
                # Independently reproduce the observer's saved decomposition.
                for key in ['mean_y','mean_z','eddy_y','eddy_z']:
                    native_key={'mean_y':'centered_mean_advection_y','mean_z':'centered_mean_advection_z',
                                'eddy_y':'centered_eddy_advection_y','eddy_z':'centered_eddy_advection_z'}[key]
                    assert np.allclose(a[key],d[native_key][:]/rho[:,None],rtol=2e-9,atol=2e-14),key
            for k,vv in a.items():
                assert np.isfinite(vv).all(),(i,k)
                if k in FIELDS:digest.update(np.asarray(vv,dtype='<f8').tobytes())
            if i in [400,600,800]:states[i]=u.copy()
            for j,(lo,hi) in enumerate(WINDOWS):
                if lo*8<=i<=hi*8:
                    weight=(.5 if i in [lo*8,hi*8] else 1)/((hi-lo)*8)
                    for k,vv in a.items():averages[j][k]=averages[j].get(k,0)+weight*vv
        if i%100==0:print(label,'read day',i/8,flush=True)
    exact=[]
    if label=='wide':
        for lo,hi in WINDOWS:
            total=np.zeros((4,len(z),len(y)));observed=np.zeros((len(z),len(y)));coverage=[]
            for seg in range(6,11):
                start=max(lo,(seg-1)*10);end=min(hi,seg*10)
                if end<=start:continue
                # The endpoint contains the segment's origin state, even when its duplicate frame was omitted.
                with h5py.File(files[end*8]) as d:
                    origin=float(d['native_budget_origin_time'][()])/86400
                    assert origin==(seg-1)*10
                    assert float(d['native_budget_end_time'][()])/86400==end
                    end_parts=d['native_integrated_parts'][:];end_state=d['native_rhou'][:]
                    if start==origin:begin_parts=np.zeros_like(end_parts);begin_state=d['native_rhou_origin'][:]
                if start!=origin:
                    with h5py.File(files[start*8]) as d:
                        assert float(d['native_budget_origin_time'][()])/86400==origin
                        begin_parts=d['native_integrated_parts'][:];begin_state=d['native_rhou'][:]
                total+=end_parts-begin_parts;observed+=end_state-begin_state;coverage.append([start,end])
            assert sum(b-a for a,b in coverage)==hi-lo
            duration=(hi-lo)*86400
            parts=total/duration/rho[None,:,None];tendency=observed/duration/rho[:,None]
            direct=(states[hi*8]-states[lo*8])/duration
            assert np.allclose(tendency,direct,rtol=1e-7,atol=1e-14),'native/centered state mismatch'
            residual=tendency-parts.sum(axis=0)
            assert np.max(abs(residual))<1e-10*max(np.max(abs(parts)),1e-15)+1e-14
            exact.append({'parts':parts,'tendency':tendency,'residual':residual,'coverage':coverage})
    return dict(label=label,meta=m,files=files,rho=rho,y=y,z=z,dy=dy,dz=dz,means=averages,
                states=states,exact=exact,qa=qa,provenance=metas,diagnostic_sha256=digest.hexdigest())

cases={k:load_case(k,p) for k,p in ROOTS.items()}
n,w=cases['narrow'],cases['wide']
for k in PHYSICAL:assert n['meta'].get(k)==w['meta'].get(k),(k,n['meta'].get(k),w['meta'].get(k))
assert n['dy']==w['dy'] and np.array_equal(n['rho'],w['rho']) and np.array_equal(n['z'],w['z'])
assert np.array_equal(n['y'][abs(n['y'])<=2400000],w['y'][abs(w['y'])<=2400000])

report={'status':'COMPLETE','compute_location':'Unity','windows_earth_days':WINDOWS,
        'common_region':'|y| <= 2400 km, 10 <= z <= 75 km; density-weighted unsmoothed norms',
        'wall_region':'Within 600 km of either wall, 10 <= z <= 75 km; each side reported separately',
        'comparison':'Matched physical parameters and spacing; different realized random channel modes; one seed per width',
        'equation':'du_bar/dt = Coriolis + mean_y + mean_z + eddy_y + eddy_z + reconstruction + sampling + diffusion + external forcing + closure residual (wide only)',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':{}}

for label,c in cases.items():
    m,rho,y,z=c['meta'],c['rho'],c['y'],c['z'];height=(z[:,None]>=10000)&(z[:,None]<=75000)
    masks={'interior':height&(abs(y[None,:])<=2400000),
           'south_wall':height&(y[None,:]<-m['ly']/2+600000),
           'north_wall':height&(y[None,:]>m['ly']/2-600000)}
    c['masks']=masks
    rec={'width_km':m['ly']/1000,'grid':[m['nx'],m['ny'],m['nz']], 'qa':c['qa'],
         'metadata':c['provenance'],'diagnostics_read_sha256':c['diagnostic_sha256'],'windows':[]}
    for j,(lo,hi) in enumerate(WINDOWS):
        a=c['means'][j];tendency=(c['states'][hi*8]-c['states'][lo*8])/((hi-lo)*86400)
        terms={key:a[key] for key in ['fv','mean_y','mean_z','eddy_y','eddy_z']}
        terms['mean_total']=a['mean_y']+a['mean_z'];terms['eddy_total']=a['eddy_y']+a['eddy_z'];terms['tendency']=tendency
        if label=='wide':
            ex=c['exact'][j]
            terms.update(reconstruction=a['reconstruction'],sampling_advection=ex['parts'][0]-a['native_snapshot'][0],
                         coriolis=ex['parts'][1],diffusion=ex['parts'][2],external_forcing=ex['parts'][3],closure=ex['residual'],
                         sampled_fv_difference=ex['parts'][1]-a['fv'])
            assert np.allclose(terms['mean_total']+terms['eddy_total']+terms['reconstruction']+terms['sampling_advection'],ex['parts'][0],rtol=1e-8,atol=1e-14)
        else:
            terms['unresolved_residual']=tendency-a['fv']-terms['mean_total']-terms['eddy_total']
        a['terms']=terms
        rs={}
        for region,mask in masks.items():
            rr={key:rms(a[key],mask,rho) for key in ['zonal_u','zonal_v','zonal_w']}
            rr['psi_rms_kg_m_s']=rms(psi_center(a['mass_streamfunction']),mask,rho)
            rr['v_persistence']=rr['zonal_v']/max(np.sqrt(np.average(a['v_square'][mask],weights=np.broadcast_to(rho[:,None],mask.shape)[mask])),1e-30)
            rr['w_persistence']=rr['zonal_w']/max(np.sqrt(np.average(a['w_square'][mask],weights=np.broadcast_to(rho[:,None],mask.shape)[mask])),1e-30)
            rr['acceleration_rms_m_s2']={k:rms(v,mask,rho) for k,v in terms.items()}
            rr['eddy_vs_minus_fv_cosine']=cosine(terms['eddy_total'],-a['fv'],mask,rho)
            rr['eddy_vs_minus_fv_cosine_smoothed']=cosine(smooth(terms['eddy_total'],c['dz'],c['dy']),smooth(-a['fv'],c['dz'],c['dy']),mask,rho)
            rr['eddy_projection_on_minus_fv']=rr['eddy_vs_minus_fv_cosine']*rms(terms['eddy_total'],mask,rho)/max(rms(a['fv'],mask,rho),1e-30)
            # Regional shear conversion, not regional force work: flux through its edges is omitted.
            rr['resolved_eddy_to_mean_shear_conversion_W_m3']=float(a['eddy_ke_conversion'][mask].mean())
            rs[region]=rr
        rec['windows'].append({'days':[lo,hi],'regions':rs})
    mask=masks['interior'];a,b=c['means'][:2]
    rec['half_window_cosines']={key:cosine(a[key],b[key],mask,rho) for key in ['zonal_u','zonal_v','zonal_w']}
    rec['half_window_cosines']['psi']=cosine(psi_center(a['mass_streamfunction']),psi_center(b['mass_streamfunction']),mask,rho)
    rec['half_window_cosines']['psi_smoothed']=cosine(smooth(psi_center(a['mass_streamfunction']),c['dz'],c['dy']),smooth(psi_center(b['mass_streamfunction']),c['dz'],c['dy']),mask,rho)
    report['cases'][label]=rec

# Independent regression against the published narrow-channel audit.
old=json.loads(Path('rendered/ferrel35-day100-audit-v2/summary.json').read_text())
ref=next(r for r in old['windows'] if r['window_earth_days']=='50–100')
got=report['cases']['narrow']['windows'][2]['regions']['interior']
for key,target in [('zonal_u','mean_u_rms_m_s'),('zonal_v','mean_v_rms_m_s'),('zonal_w','mean_w_rms_m_s')]:
    assert np.isclose(got[key],ref[target],rtol=1e-11,atol=1e-14),(key,got[key],ref[target])
report['independent_checks']=['Published narrow day50–100 u/v/w RMS reproduced to rtol 1e-11',
    'Wide saved centered terms independently reconstructed at every snapshot',
    'Wide cumulative native integrals cover all window seconds without overlap or gaps',
    'Native momentum increments agree with independent endpoint zonal velocities',
    'Identical central y coordinates, density, grid spacing and listed physical parameters']
report['limitations']=[
    'Exact all-step native momentum integrals exist only for the wide run; narrow unresolved residual mixes omitted native terms, reconstruction and sampling.',
    'Wide mean/eddy decomposition is centered and sampled every three hours; its reconstruction and sampling errors are explicitly retained, not interpreted as eddy forcing.',
    'Wall association is descriptive, not unique causal attribution: width also changes allowed channel modes and total f range; one realization per width.',
    'Half-window persistence is not a confidence interval; the flow need not be statistically stationary.',
    'No residual/Lagrangian circulation or thermodynamic indirectness test; passive dye is not ammonia chemistry.',
    'The prescribed cap, bottom and thermal feedback remain unchanged; horizontal resolution is not established as converged.']

# Full-width views use common units, limits, and smoothing widths; statistics remain raw.
fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
for row,key in enumerate(['mass_streamfunction','zonal_u']):
    fields=[smooth(c['means'][2][key],c['dz'],c['dy']) for c in cases.values()]
    scale=1000 if row==0 else 1;lim=max(abs(a).max() for a in fields)/scale
    for col,(label,c) in enumerate(cases.items()):
        a=fields[col]/scale
        yy=(np.arange(c['meta']['ny']+1)*c['dy']-c['meta']['ly']/2)/1000 if row==0 else c['y']/1000
        zz=np.arange(c['meta']['nz']+1)*c['dz']/1000 if row==0 else c['z']/1000
        im=axs[row,col].pcolormesh(yy,zz,a,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
        if row==0:axs[row,col].contour(yy,zz,a,levels=np.linspace(-lim,lim,11),colors='k',linewidths=.35)
        for x in [-2400,2400]:axs[row,col].axvline(x,color='k',ls=':',lw=.7)
        axs[row,col].axhline(80,color='k',ls='--',lw=.7)
        axs[row,col].set(title=f'{label.capitalize()} · '+('mass streamfunction' if row==0 else 'zonal velocity'),xlabel='y north [km]',ylabel='Height [km]')
        fig.colorbar(im,ax=axs[row,col],label='Ψ / 1000 [kg m⁻¹ s⁻¹]' if row==0 else 'u [m/s]')
fig.suptitle('35°N · days 50–100 · same grid spacing and forcing\nDisplay smoothing: 200 km in y, 5 km in z. Dotted lines: common central region; dashed: cap base.')
fig.savefig(OUT/'circulation.png',dpi=140);plt.close(fig)

fig,axs=plt.subplots(2,3,figsize=(15,8),layout='constrained')
ps=[smooth(psi_center(c['means'][j]['mass_streamfunction']),c['dz'],c['dy']) for c in cases.values() for j in range(3)]
lim=max(abs(a[:,abs(c['y'])<=2400000]).max() for c in cases.values() for a in [smooth(psi_center(c['means'][j]['mass_streamfunction']),c['dz'],c['dy']) for j in range(3)])/1000
for row,(label,c) in enumerate(cases.items()):
    take=abs(c['y'])<=2400000
    for col,win in enumerate(WINDOWS):
        a=smooth(psi_center(c['means'][col]['mass_streamfunction']),c['dz'],c['dy'])/1000
        im=axs[row,col].pcolormesh(c['y'][take]/1000,c['z']/1000,a[:,take],cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
        axs[row,col].set(ylim=(10,75),xlabel='y north [km]',ylabel='Height [km]',title=f'{label.capitalize()} · days {win[0]}–{win[1]}')
        fig.colorbar(im,ax=axs[row,col],label='Ψ / 1000 [kg m⁻¹ s⁻¹]')
fig.suptitle('Does interior overturning persist? Same physical region and common scale\nDisplay smoothing: 200 km in y, 5 km in z; each panel is an independently averaged window.')
fig.savefig(OUT/'interior_persistence.png',dpi=140);plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
profile_report={}
for label,c in cases.items():
    a=c['means'][2];zh=(c['z']>=10000)&(c['z']<=75000);weights=c['rho'][zh];weights/=weights.sum()
    profs=[np.sqrt(np.sum(weights[:,None]*field[zh]**2,axis=0)) for field in [a['zonal_u'],psi_center(a['mass_streamfunction']),a['zonal_v'],a['terms']['eddy_total']]]
    side_records={}
    for side,sgn in [('south',-1),('north',1)]:
        sel=sgn*c['y']>0;dist=c['meta']['ly']/2-sgn*c['y'];order=np.argsort(dist[sel]);dd=dist[sel][order]/1000
        side_records[side]={}
        for ax,prof,name in zip(axs.flat,profs,['u','psi','v','eddy_acceleration']):
            val=prof[sel][order]
            # Fixed 200-km bins; no averaging across the two walls.
            bins=np.arange(0,c['meta']['ly']/2000+200,200);bx=[];by=[]
            for low,high in zip(bins[:-1],bins[1:]):
                use=(dd>=low)&(dd<high)
                if use.any():bx.append(float(dd[use].mean()));by.append(float(val[use].mean()))
            ax.plot(bx,by,label=f'{label} {side}',color='C0' if label=='narrow' else 'C1',ls='-' if side=='south' else '--')
            peak=int(np.argmax(by));side_records[side][name]={'peak_distance_km':bx[peak],'peak_value':by[peak],'distance_km':bx,'rms':by}
    profile_report[label]=side_records
for ax,title,unit in zip(axs.flat,['Zonal velocity','Mass streamfunction','Meridional velocity','Resolved eddy acceleration'],['u RMS [m/s]','Ψ RMS [kg m⁻¹ s⁻¹]','v RMS [m/s]','acceleration RMS [m/s²]']):
    ax.axvline(600,color='gray',ls=':',lw=1);ax.set(title=title,xlabel='Distance inward from each wall [km]',ylabel=unit);ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Days 50–100 · separate north/south wall profiles\nDensity-weighted RMS over heights 10–75 km, then 200-km distance bins; vertical line = 600-km wall strip.')
fig.savefig(OUT/'wall_distance.png',dpi=140);plt.close(fig)
report['wall_profiles']=profile_report

c=cases['wide'];terms=c['means'][2]['terms'];mask=c['masks']['interior'];rho=c['rho']
fig,axs=plt.subplots(3,3,figsize=(15,10),layout='constrained')
names=['coriolis','eddy_total','mean_total','reconstruction','sampling_advection','diffusion','external_forcing','tendency','closure']
take=abs(c['y'])<=2400000;zh=(c['z']>=10000)&(c['z']<=75000)
display=[smooth(terms[k],c['dz'],c['dy'])[np.ix_(zh,take)]*1e6 for k in names]
lim=max(abs(a).max() for a in display)
for ax,key,a in zip(axs.flat,names,display):
    im=ax.pcolormesh(c['y'][take]/1000,c['z'][zh]/1000,a,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
    ax.set(title=f'{key.replace("_"," ")}\nraw RMS {rms(terms[key],mask,rho):.2e} m/s²',xlabel='y north [km]',ylabel='Height [km]')
fig.colorbar(im,ax=list(axs.flat),label='Acceleration [10⁻⁶ m/s²]')
fig.suptitle('Wide channel, interior momentum balance · days 50–100\nCommon display scale; 200-km / 5-km smoothing. Exact native total; sampled centered mean/eddy split.')
fig.savefig(OUT/'momentum_balance.png',dpi=140);plt.close(fig)
(OUT/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in OUT.iterdir() if p.is_file()}
(OUT/'PASS.json').write_text(json.dumps({'status':'PASS','compute_location':'Unity','frames_read_per_case':401,'checks':report['independent_checks'],'assets':assets},indent=2)+'\n')
print(json.dumps({k:{'days50_100':v['windows'][2]['regions'],'persistence':v['half_window_cosines']} for k,v in report['cases'].items()},indent=2),flush=True)
