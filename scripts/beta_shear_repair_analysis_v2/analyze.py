"""Analyze finished, bounded repair diagnostics on Unity; publish only aggregates."""
import csv,json,sys,tomllib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True)
roots=[Path(p) for p in sys.argv[2:]]
def rows(p):return [{k:float(v) if k!='event' else v for k,v in r.items()} for r in csv.DictReader(p.open())]
def toml(p):return tomllib.loads(p.read_text())
result={'scope':'Bounded thermal/energy/positivity/jet-screen diagnostics; no production release.','roots':[str(r) for r in roots],'benchmarks':{},'replays':{},'screens':{},'profiles':{}}
for root in roots:
    if (root/'DRIVER_PASS.toml').exists():result.setdefault('validation',{})[str(root)]=toml(root/'DRIVER_PASS.toml')
    for file,key in [('NATIVE_ENERGY.toml','native_energy'),('JET_SCREENS.toml','screens'),('PROFILE_OUTCOMES.toml','profiles')]:
        if (root/file).exists():result.setdefault(key,{}).update(toml(root/file))
    if (root/'HEATED_AUDIT.toml').exists():result['heated_audit']=toml(root/'HEATED_AUDIT.toml')
    for case in ['cap_anticyclonic','cap_cyclonic','traditional_anticyclonic','traditional_cyclonic']:
        p=root/case
        if (p/'diagnostics.csv').exists():
            result['profiles'].setdefault(case,{}).update(samples=rows(p/'diagnostics.csv'))
            if (p/'GUARD_STOP.toml').exists():result['profiles'][case]['guard']=toml(p/'GUARD_STOP.toml')
    if (root/'spectrum.json').exists():result['spectrum']=json.loads((root/'spectrum.json').read_text())
    for case in ['weno','centered_momentum','centered_both']:
        p=root/case
        if (p/'diagnostics.csv').exists():
            rr=rows(p/'diagnostics.csv');result['benchmarks'][case]={'metadata':toml(p/'metadata.toml'),'samples':rr,'max_abs_energy_defect_over_initial_ke':max(abs(x['energy_defect_over_initial_ke']) for x in rr)}
    for case in ['original','restored_cap']:
        p=root/case
        if not (p/'stage_terms.csv').exists():continue
        rr=rows(p/'stage_terms.csv');item={'metadata':toml(p/'metadata.toml')}
        bad=next((i for i,r in enumerate(rr) if r['event']=='state' and (r['q']<=0 or r['rho']<=0 or not np.isfinite(r['T']))),None)
        if bad is not None:
            state=rr[bad];coord=tuple(state[k] for k in ['i','j','k'])
            prev=next((r for r in reversed(rr[:bad]) if r['event']=='tendency' and tuple(r[k] for k in ['i','j','k'])==coord),None)
            stage_dt=state['last_stage_dt']*{2:1.,3:3.,1:2.}[int(state['next_stage'])]
            detail={'state':state,'preceding_tendency':prev,'effective_stage_interval':stage_dt,
                    'slow_increment':stage_dt*state['G'],'acoustic_remainder':state['q']-state['U0']-stage_dt*state['G']}
            if prev:
                detail['slow_components_sum']=sum(prev[k] for k in ['advection','viscous_source','entropy_source','external_source'])
                detail['component_closure_error']=prev['G']-detail['slow_components_sum']
            item['first_bad_in_recorded_patch']=detail
        else:item['first_bad_in_recorded_patch']=None
        if (p/'GUARD_STOP.toml').exists():item['guard']=toml(p/'GUARD_STOP.toml')
        if (p/'diagnostics.csv').exists():item['samples']=rows(p/'diagnostics.csv')
        result['replays'][case]=item
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    return None if isinstance(x,float) and not np.isfinite(x) else x
result=clean(result);(out/'repair-summary.json').write_text(json.dumps(result,indent=2)+'\n')
plt.rcParams['svg.fonttype']='none';fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for name,v in result['benchmarks'].items():
    rr=v['samples'];axes[0,0].plot([r['time'] for r in rr],[100*r['energy_defect_over_initial_ke'] for r in rr],label=name.replace('_',' '))
axes[0,0].set(title='Unheated energy benchmark',xlabel='Nondimensional time',ylabel='Unaccounted energy / initial KE (%)');axes[0,0].legend(fontsize=8)
for sign,label in [(1,'Anticyclonic central shear'),(-1,'Cyclonic central shear')]:
    xx=[];yy=[]
    for speed in [.015625,.03125,.0625,.125,.25]:
        v=result['screens'].get(f'cap_{sign*speed}')
        if v:xx.append(speed);yy.append(v['cap_PV_min'])
    axes[0,1].plot(xx,yy,'o-',label=label)
axes[0,1].axhline(0,color='k',lw=.7);axes[0,1].set(title='Balanced cap jets: PV is not guaranteed',xlabel='Jet amplitude',ylabel='Minimum cap (curl u + 2Ω) · ∇ ln θ');axes[0,1].legend(fontsize=8)
for name,v in result['replays'].items():
    rr=v.get('samples',[])
    if rr:axes[1,0].plot([r['time'] for r in rr],[r['Tmin'] for r in rr],label=name.replace('_',' '))
axes[1,0].set(title='Late-state forensic replay (already rejected flow)',xlabel='Nondimensional time',ylabel='Minimum temperature');axes[1,0].legend(fontsize=8);axes[1,0].ticklabel_format(useOffset=False,axis='x')
sp=result.get('spectrum',{}).get('cases',{})
for name in ['reference','deep_anti','cap_anti','cap_cyclone','traditional_anti','traditional_cyclone']:
    pts=[]
    for ny,nz in [(128,80),(256,80)]:
        v=sp.get(f'{name}_{ny}_{nz}')
        if v and v['eigenvalues']:pts.append((ny,max(x['real'] for x in v['eigenvalues'])))
    if pts:axes[1,1].plot(*zip(*pts),'o-',label=name)
axes[1,1].set(title='Separate inviscid cross-section eigenvalue screen',xlabel='Horizontal cells (80 vertical cells)',ylabel='Largest real eigenvalue found');axes[1,1].legend(fontsize=7)
for ax in axes.flat:ax.grid(alpha=.2)
fig.suptitle('Repair validation — diagnostics, not an LSV/GRS production result')
fig.savefig(out/'repair-checks.svg');fig.savefig(out/'repair-checks.png',dpi=120)
print('ANALYSIS_COMPLETE',out,flush=True)
