"""Consolidate existing published convection evidence on Unity; no refitting."""
from pathlib import Path
import json,hashlib,sys

sign_path=Path('rendered/fplane-sign-v1/analysis.json')
fit_path=Path('rendered/fplane-closure-v1-1939734/report.json')
sign=json.loads(sign_path.read_text());fit=json.loads(fit_path.read_text())
assert sign['status']=='COMPLETE' and fit['status']=='COMPLETE'
models={}
for name,filtered in fit['filters'].items():
    models[name]={}
    for model in ('linear','polynomial_free','polynomial_sign_saturating'):
        d=filtered['models'][model]
        models[name][model]={'coefficients':d['coefficients'],
            'heldout_interior':d['heldout_interior'],
            'heldout_central_third_moment_source':d['heldout_source_contributions']['central_m3_source_dimless'],
            'heldout_low_mode_energy_work':d['heldout_source_contributions']['low_mode_energy_work_m2_s3']}
report={'status':'PHYSICAL_CLOSURE_NOT_VALIDATED','compute_location':'Unity','new_simulations':False,'new_fit':False,
        'source_reports':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (sign_path,fit_path)},
        'training_days':fit['training_days'],'heldout_days':fit['test_days'],'fit_region':fit['fit_region'],
        'full_channel_attribution':{k:{p:v[p] for p in ('unassigned_over_convective_central_moment','physical_attribution_screen_pass')} for k,v in sign['cases'].items()},
        'interior_native_audit':{k:{p:v[p] for p in ('native_advection_reconstruction_rms_over_C','interior_unassigned_over_convective_moment')} for k,v in fit['native_six_hour_audit'].items()},
        'heldout_models':models,
        'conclusions':['The existing held-out test is reused rather than repeating a fit.',
            'A quadratic term captures part of the central-third-moment production, but its spatial source prediction barely improves on the linear baseline.',
            'The unconstrained cubic coefficient has the opposite sign to damping; imposing cubic damping drives its fitted coefficient to approximately zero.',
            'The fitted models do not predict the measured low-mode source work accurately, including after the declared filter.',
            'Native budget closure and physical stress attribution are different tests; the latter still fails.',
            'No fitted coefficients are transplanted into the 2D control and no new 3D job is launched.'],
        'physical_next_gate':'An independently trustworthy signed and scale-dependent convective stress-work prediction, including its numerical-reconstruction uncertainty; reproduction of skewness alone is insufficient.'}
out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2,allow_nan=False))
print(json.dumps({'status':report['status'],'output':str(out),'new_fit':False}))
