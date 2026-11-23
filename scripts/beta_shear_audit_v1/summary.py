"""Postprocess on Unity only; publish aggregates, never model fields."""
import csv, json, pathlib, tomllib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=pathlib.Path('data')
out=R/'beta-shear-deep-audit';out.mkdir(exist_ok=True)
def rows(p):
    return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(p.open())]
def yz(case,t):
    with h5py.File(R/'beta-shear-linear-v1-2138365'/case/f'yz_{100*t:07d}.jld2') as h:
        return {k:h[k][()] for k in ['u','v','w','rho','T','dy','dz']}
result={'scope':'Diagnostic audit of rejected shear pilot. No production release.',
        'jobs':{'native_operators':2138128,'paired_refinement':2138365,'forensic_replay':2138558}}
energy=[]
for p in sorted((R/'beta-shear-native-audit-2138128').glob('*.toml')):
    a=tomllib.loads(p.read_text());name=p.stem.rsplit('_',1)[0]
    root=R/'beta-shear-v1-pilot-2134379'/'parent' if name=='original_3d' else R/'beta-shear-diag-v2-2136871'/name
    d=rows(root/'diagnostics.csv');times=np.array([r['time'] for r in d]);i=int(np.argmin(abs(times-a['time'])))
    lo=max(0,i-1);hi=min(len(d)-1,i+1)
    fd=(d[hi]['total_energy']-d[lo]['total_energy'])/(times[hi]-times[lo])
    a.update(case=name,measured_energy_derivative=fd,comparison_window=[float(times[lo]),float(times[hi])],
             kinetic_advection_defect=a['momentum_advection_work']+a['advective_mass_kinetic'],
             thermodynamic_pressure_defect=a['thermodynamic_advection']+a['potential_energy_advection']+a['pressure_gravity_work'])
    a['relative_derivative_mismatch']=abs(a['predicted_energy_derivative']-fd)/max(abs(fd),1e-20)
    energy.append(a)
result['energy_audits']=energy
orig=rows(R/'beta-shear-v1-pilot-2134379/parent/diagnostics.csv')
result['original_end']={k:orig[-1][k] for k in ['time','min_temperature','max_mach','total_energy','kinetic_energy','cumulative_energy_defect_since_restart','energy_defect_over_restart_ke']}
result['original_initial_energy']=orig[0]['total_energy'];result['original_initial_ke']=orig[0]['kinetic_energy']
curves={};paired={}
for variant in ['coarse','fine','traditional']:
    samples=[]
    for t in range(31):
        b=yz(variant+'_base',t);a=yz(variant+'_probe',t);n=b['rho'].size
        delta={k:a[k]-b[k] for k in ['u','v','w']}
        vw=np.sqrt((np.sum(delta['v']**2)+np.sum(delta['w']**2))/n)
        allv=np.sqrt(sum(np.sum(v*v) for v in delta.values())/n)
        s={'time':t,'cross_plane_difference_rms':float(vw),'velocity_difference_rms':float(allv),
           'temperature_difference_rms':float(np.sqrt(np.mean((a['T']-b['T'])**2)))}
        if variant=='coarse':
            c=yz('coarse_double',t)
            s['doubling_relative_error']=float(np.sqrt(sum(np.sum((c[k]-b[k]-2*delta[k])**2) for k in delta))/max(2*allv*np.sqrt(n),1e-30))
        samples.append(s)
    fit={}
    for lo,hi in [(10,25),(15,30)]:
        x=np.arange(lo,hi+1);y=np.log([samples[t]['velocity_difference_rms'] for t in x]);slope,intercept=np.polyfit(x,y,1)
        residual=y-(slope*x+intercept);ss=np.sum((y-y.mean())**2)
        fit[f'{lo}-{hi}']={'amplitude_growth_rate':float(slope),'r_squared':float(1-np.sum(residual**2)/ss)}
    paired[variant]={'samples':samples,'fits':fit,'amplification_t30':samples[-1]['velocity_difference_rms']/samples[0]['velocity_difference_rms']}
    curves[variant]=np.array([s['velocity_difference_rms'] for s in samples])
result['paired_tests']=paired
# Independent directional-derivative check of the discrete energy functional.
# Native face momenta, arithmetic density interpolation, closed normal walls.
a=yz('coarse_base',0);rho=a['rho'];q=(rho*a['T'])**.6;nz,ny=rho.shape
def faces(r):return [r,np.concatenate([r[:,:1],.5*(r[:,:-1]+r[:,1:]),r[:,-1:]],1),np.concatenate([r[:1],.5*(r[:-1]+r[1:]),r[-1:]],0)]
rf=faces(rho);M=[r*a[k] for r,k in zip(rf,['u','v','w'])]
rng=np.random.default_rng(1234);dr=rho*.01*rng.normal(size=rho.shape);dq=q*.01*rng.normal(size=q.shape)
dm=[r*.01*rng.normal(size=r.shape) for r in rf];dm[1][:,[0,-1]]=0;dm[2][[0,-1],:]=0
z=(np.arange(nz)+.5)[:,None]*float(a['dz']);vol=float(a['dy'])*float(a['dz'])
def E(r,q,m):return vol*(1.5*np.sum(q**(5/3))+np.sum(12*z*r)+sum(np.sum(mi*mi/(2*ri)) for mi,ri in zip(m,faces(r))))
drf=faces(dr)
analytic=vol*(2.5*np.sum(q**(2/3)*dq)+np.sum(12*z*dr)+sum(np.sum(mi*di/ri-.5*mi*mi*dri/ri**2) for mi,di,ri,dri in zip(M,dm,rf,drf)))
checks=[]
for eps in [1e-2,1e-3,1e-4]:
    fd=(E(rho+eps*dr,q+eps*dq,[m+eps*d for m,d in zip(M,dm)])-E(rho-eps*dr,q-eps*dq,[m-eps*d for m,d in zip(M,dm)]))/(2*eps)
    checks.append({'epsilon':eps,'analytic':float(analytic),'finite_difference':float(fd),'relative_error':float(abs(fd-analytic)/abs(analytic))})
assert min(x['relative_error'] for x in checks)<1e-6
result['energy_functional_derivative_check']=checks
replay=R/'beta-shear-forensic-v1-2138558/replay'
if (replay/'step_health.csv').exists():
    health=rows(replay/'step_health.csv');bad=[r for r in health if not all(np.isfinite(v) for v in r.values()) or min(r[k] for k in ['min_rho','min_rhotheta','min_temperature','min_pressure'])<=0]
    result['forensic']={'health_last_rows':health[-8:],'first_bad_health':bad[0] if bad else None}
    if (replay/'FIRST_BAD_STATE.toml').exists():result['forensic']['capture']=tomllib.loads((replay/'FIRST_BAD_STATE.toml').read_text())
    if (replay/'FIRST_BAD_STATE.jld2').exists():
        fields={}
        with h5py.File(replay/'FIRST_BAD_STATE.jld2') as h:
            for k in ['rho','rhotheta','T','P','u','v','w','nu']:
                v=h[k][:];finite=np.isfinite(v);nonpositive=finite&(v<=0);badmask=~finite | (nonpositive if k in ['rho','rhotheta','T','P'] else False)
                idx=np.argwhere(badmask)
                fields[k]={'shape':list(v.shape),'nonfinite_count':int((~finite).sum()),'nonpositive_count':int(nonpositive.sum()),'finite_min':float(v[finite].min()) if finite.any() else None,'finite_max':float(v[finite].max()) if finite.any() else None,'first_bad_zyx_indices':idx[:6].tolist(),'bad_z_index_range':[int(idx[:,0].min()),int(idx[:,0].max())] if len(idx) else None}
        result['forensic']['fields']=fields
    if (replay/'diagnostics.csv').exists():
        rr=rows(replay/'diagnostics.csv');result['forensic']['diagnostic_last']=rr[-1]
        oo=min(orig,key=lambda r:abs(r['time']-150));rr0=rr[0]
        result['forensic']['restart_comparison']={k:{'original':oo[k],'replay':rr0[k]} for k in ['time','total_energy','min_temperature','max_mach','kinetic_energy']}
def clean(v):
    if isinstance(v,dict):return {k:clean(a) for k,a in v.items()}
    if isinstance(v,list):return [clean(a) for a in v]
    return None if isinstance(v,float) and not np.isfinite(v) else v
result=clean(result);(out/'audit-summary.json').write_text(json.dumps(result,indent=2)+'\n')
fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
for name,val in curves.items():axs[0,0].semilogy(range(31),val/val[0],label=name)
axs[0,0].set(title='Paired perturbation amplification',xlabel='Nondimensional time',ylabel='Velocity difference / initial difference');axs[0,0].legend()
axs[0,1].semilogy(range(31),[max(s['doubling_relative_error'],1e-16) for s in paired['coarse']['samples']]);axs[0,1].set(title='Does doubling the probe double the response?',xlabel='Nondimensional time',ylabel='Relative departure from linearity')
e=sorted([a for a in energy if a['case']=='original_3d'],key=lambda a:a['time']);x=np.arange(len(e));w=.34
axs[1,0].bar(x-w/2,[a['measured_energy_derivative'] for a in e],w,label='Recorded dE/dt');axs[1,0].bar(x+w/2,[a['predicted_energy_derivative'] for a in e],w,label='Native-operator audit');axs[1,0].set_xticks(x,[str(int(a['time'])) for a in e]);axs[1,0].set(title='Original 3D pilot: energy loss accounted for',xlabel='Nondimensional time',ylabel='Total-energy rate');axs[1,0].legend()
ax=axs[1,1];ax.plot([r['time'] for r in orig],[r['min_temperature'] for r in orig],color='tab:blue',label='Minimum T');ax.set(xlabel='Nondimensional time',ylabel='Minimum T',title='Rejected original pilot');ax2=ax.twinx();ax2.plot([r['time'] for r in orig],[r['max_mach'] for r in orig],color='tab:red',label='Maximum Mach');ax2.axhline(.5,color='tab:red',ls=':');ax2.set_ylabel('Maximum Mach');ax.axvline(51.5,color='gray',ls='--');ax.text(55,.8,'Acceptance limit crossed',fontsize=9)
fig.suptitle('Shear-pilot diagnostic audit — not a validated LSV experiment',fontsize=16)
fig.savefig(out/'audit-summary.png',dpi=150);plt.close(fig)
print(json.dumps({'energy':[(a['case'],a['time'],a['predicted_energy_derivative'],a['measured_energy_derivative']) for a in energy],'paired':{k:{'amplification':v['amplification_t30'],'fits':v['fits']} for k,v in paired.items()},'doubling_final':paired['coarse']['samples'][-1]['doubling_relative_error'],'forensic':result.get('forensic'), 'derivative_check':checks,'original_end':result['original_end']},indent=2))
