"""Read existing fields on Unity; no trajectory is advanced or rewritten."""
from pathlib import Path
import sys, csv, json, hashlib, tomllib, time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'feedback_maintenance_v1'))
import analyze_existing as old
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TERMS=('advection','linear','quadratic','cubic','viscous')
integrate=np.trapezoid if hasattr(np,'trapezoid') else np.trapz

def main(root, out, gate_path):
    started=time.monotonic();gate=tomllib.loads(gate_path.read_text())
    assert gate['status']=='PASS' and gate['native_reconstruction_max_error']<1e-9
    out.mkdir(parents=True,exist_ok=False)
    old.self_checks();summary={};histories={};score=np.zeros((old.N,old.N));modal=np.zeros_like(score)
    idx=(-np.arange(old.N))%old.N;low_coords=np.argwhere(old.LOW)
    for case in old.CASES:
        frames,scalars,config,provenance=old.inventory(root,case)
        rows=[];cumulative={k:0. for k in TERMS};previous=None;previous_t=None
        threshold=2*np.sqrt(np.mean(old.load_field(frames[0.])**2));max_scalar_error=0.
        for t,path in sorted(frames.items()):
            q=old.load_field(path);qh=np.fft.fft2(q);u,v=old.velocity(q)
            au,av=old.centered_advection(u,v)
            sources={'advection':old.curl(au,av),'quadratic':config['b']*config['f']*q*q,'cubic':-config['c']*q*q*q}
            sh={k:np.fft.fft2(v-v.mean()) for k,v in sources.items()}
            sh['viscous']=-old.NU*old.LAM*qh
            sh['linear']=config['a']*old.K2*old.SHELL*qh
            work={k:(qh.conj()*v).real*old.INVLAM/old.NORM for k,v in sh.items()}
            ek=.5*np.abs(qh)**2*old.INVLAM/old.NORM
            row={'t':t,'energy':float(ek.sum()),'energy_k_le_2':float(ek[old.LOW].sum()),
                 'skewness':float(np.mean(q**3)/np.mean(q*q)**1.5),'mean_zeta':float(q.mean()),
                 'energy_weighted_wavelength':float(2*np.pi*np.sum(np.divide(ek,np.sqrt(old.K2),out=np.zeros_like(ek),where=old.K2>0))/ek.sum())}
            for k in TERMS:
                row[k+'_lowk_power']=float(work[k][old.LOW].sum())
                row[k+'_total_power']=float(work[k].sum())
            for cutoff in (2,4,8,12,16,32):
                inside=(old.K2>0)&(old.K2<=cutoff**2)
                row[f'advective_power_k_le_{cutoff}']=float(work['advection'][inside].sum())
            if previous is not None:
                for k in TERMS:cumulative[k]+=.5*(t-previous_t)*(row[k+'_lowk_power']+previous[k+'_lowk_power'])
            row.update({'integral_'+k:float(v) for k,v in cumulative.items()})
            row['lowk_change']=row['energy_k_le_2']-scalars[0.]['energy_k_le_2']
            row['sampling_residual']=row['lowk_change']-sum(cumulative.values())
            for key,reference in [('energy','energy'),('energy_k_le_2','energy_k_le_2'),('quadratic_lowk_power','quadratic_lowk_energy_rate'),('cubic_lowk_power','cubic_lowk_energy_rate')]:
                value=scalars[t][reference]
                assert np.isclose(row[key],value,rtol=2e-5,atol=2e-9),(case,t,key)
                max_scalar_error=max(max_scalar_error,abs(row[key]-value))
            if t%5==0:
                morphology=old.components(q,np.sign(config['f']),threshold)
                row['major_core_count']=morphology['count_each_carrying_at_least_five_percent_positive_circulation']
                row['largest_core_area_fraction']=morphology['largest_area_fraction']
            else:
                row['major_core_count']=float('nan');row['largest_core_area_fraction']=float('nan')
            if case=='biased_plus' and t<=200 and t%5==0:
                z=qh/old.N**2;modal+=ek
                reverse=z[np.ix_(idx,idx)]
                for i,j in low_coords:
                    pairs=config['b']*config['f']*np.real(np.conj(z[i,j])*z*np.roll(reverse,(int(i),int(j)),axis=(0,1)))*old.INVLAM[i,j]
                    score+=np.abs(pairs)
            rows.append(row);previous=row;previous_t=t
        with (out/(case+'_formation.csv')).open('w') as fp:
            w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        windows={}
        for start,end in ((0,80),(80,160),(160,240),(240,320),(320,400),(0,200),(200,400),(0,400)):
            sample=[r for r in rows if start<=r['t']<=end];times=[r['t'] for r in sample]
            values={k:float(integrate([r[k+'_lowk_power'] for r in sample],times)) for k in TERMS}
            change=sample[-1]['energy_k_le_2']-sample[0]['energy_k_le_2']
            # A second quadrature at interval 1 estimates temporal-sampling sensitivity.
            coarser=sample[::2];coarse={k:float(integrate([r[k+'_lowk_power'] for r in coarser],[r['t'] for r in coarser])) for k in TERMS}
            residual=change-sum(values.values())
            windows[f'{start}_{end}']={'energy_change':change,'integrated_work':values,'sampling_residual':residual,
                'residual_over_sum_abs_work':abs(residual)/max(sum(abs(v) for v in values.values()),1e-30),
                'quadrature_change_1_vs_half':{k:coarse[k]-values[k] for k in TERMS}}
        summary[case]={'windows':windows,'max_absolute_scalar_reconstruction_error':max_scalar_error,'provenance':provenance,
                       'final':{k:v for k,v in rows[-1].items() if np.isfinite(v)}}
        histories[case]=rows
        print('FORMATION_ANALYZED',case,json.dumps(windows['0_400']),flush=True)
    # Conjugate-pair ranking uses only t<=200. Later snapshots are held out.
    ranked=[]
    for i,x in enumerate(old.K.astype(int)):
        for j,y in enumerate(old.K.astype(int)):
            if (x>0 or (x==0 and y>0)) and 0<x*x+y*y<=64**2:
                ii,jj=(-i)%old.N,(-j)%old.N
                ranked.append({'kx':int(x),'ky':int(y),'absolute_quadratic_lowmode_work_score':float(score[i,j]+score[ii,jj]),
                               'energy_score':float(modal[i,j]+modal[ii,jj]),'forcing_band':bool(old.SHELL[i,j])})
    ranked.sort(key=lambda a:a['absolute_quadratic_lowmode_work_score'],reverse=True)
    (out/'mode_ranking.json').write_text(json.dumps({'training_times':[0,200],'held_out_times':[210,400],
       'ranking':'sum of absolute signed pair contributions to quadratic low-mode energy work; not a fitted closure coefficient',
       'max_candidate_wavenumber':64,'pairs':ranked},indent=2))
    fig,axs=plt.subplots(3,2,figsize=(13,11),layout='constrained')
    for case,rows in histories.items():
        tt=[r['t'] for r in rows]
        for ax,key in zip(axs[0],('energy_k_le_2','energy_weighted_wavelength')):ax.plot(tt,[r[key] for r in rows],label=case)
    rows=histories['biased_plus'];tt=[r['t'] for r in rows]
    for k in TERMS:
        axs[1,0].plot(tt,[r[k+'_lowk_power'] for r in rows],label=k)
        axs[1,1].plot(tt,[r['integral_'+k] for r in rows],label=k)
    axs[1,1].plot(tt,[r['lowk_change'] for r in rows],color='black',label='measured change')
    axs[2,0].plot(tt,[r['sampling_residual'] for r in rows],label='unassigned sampling residual')
    for cutoff in (2,4,8,12,16,32):axs[2,1].plot(tt,[r[f'advective_power_k_le_{cutoff}'] for r in rows],label=f'k≤{cutoff}')
    for ax,title in zip(axs.flat,('Absolute low-mode energy','Energy-weighted wavelength','Instantaneous low-mode work: +f','Cumulative work: +f','Budget residual: +f','Advective transfer into modes below cutoff')):
        ax.set(title=title,xlabel='Dimensionless time');ax.legend(fontsize=7);ax.grid(alpha=.2);ax.axhline(0,color='.5',lw=.5)
    fig.savefig(out/'formation_budget.png',dpi=150);plt.close(fig)
    bridge=json.loads(Path('rendered/fplane-sign-v1/analysis.json').read_text())
    audit={'status':'NOT_IDENTIFIABLE_FROM_EXISTING_ATTRIBUTION','native_budget_closed':all(c['validation']['max_step_relative_residual_full']<1e-8 for c in bridge['cases'].values()),
           'source_report_sha256':hashlib.sha256(Path('rendered/fplane-sign-v1/analysis.json').read_bytes()).hexdigest(),
           'cases':{k:{p:v[p] for p in ('f','unassigned_over_convective_central_moment','convective_sign_persists_across_blocks','physical_attribution_screen_pass')} for k,v in bridge['cases'].items()},
           'conclusion':'Numerical closure alone does not identify convective signed feedback. No polynomial coefficient fitted to these contaminated source partitions.',
           'next_physical_test_requirement':'Resolve the source-partition error relative to the proposed feedback, then predict signed scale-dependent work in a held-out interval or rotation case. The inherited-state f flip is not an independent formation test.',
           'new_3d_runs_launched':False}
    (out/'convection_bridge_audit.json').write_text(json.dumps(audit,indent=2))
    report={'status':'COMPLETE','compute_location':'Unity','native_operator_validation':gate,'cases':summary,
            'convection_bridge':audit,'elapsed_seconds':time.monotonic()-started,
            'limitations':['Float32 saved fields, native-validated staggered reconstruction.','Half-unit sampled work, not every-stage accounting; temporal quadrature sensitivity and residual explicitly reported.',
              'Late-time mode holdout tests a selected reduction; it is not an independent seed or a convection validation.','No inverse-cascade conclusion from morphology alone.']}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print('FORMATION_AUDIT_COMPLETE',flush=True)

if __name__=='__main__':main(*map(Path,sys.argv[1:]))
