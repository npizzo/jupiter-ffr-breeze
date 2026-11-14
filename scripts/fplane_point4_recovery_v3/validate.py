"""Finish the previously unreachable analysis gate, reusing verified GPU pilots."""
from pathlib import Path
import json,tomllib,hashlib,subprocess,sys,datetime
root=Path('data/fplane-point4-v2/preflight');source=Path('scripts/fplane_point4_v2')
assert not (root/'PASS.toml').exists(), 'Already released; do not repeat'
gate=tomllib.loads((root/'GPU_PASS.toml').read_text())
assert gate['status']=='PASS' and gate['compute_location']=='Unity'
assert gate['fine_ten_day_estimated_gpu_hours']<gate['maximum_estimated_hours']==12
assert gate['allocation_hours']==16
for name,h in gate['source_hashes'].items():assert hashlib.sha256((source/name).read_bytes()).hexdigest()==h,name
for name,h in gate['config_hashes'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h,name
assert hashlib.sha256((root/'day90_heated_bundle.jld2').read_bytes()).hexdigest()==gate['bundle_sha256']
subprocess.run([sys.executable,str(source/'test_diagnostics.py'),str(root)],check=True)
report=json.loads((root/'DIAGNOSTICS_PASS.json').read_text());assert report['status']=='PASS'
freeze=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True)
(root/'analysis-environment.txt').write_text(freeze)
receipt={'status':'PASS','compute_location':'Unity','GPU_preflight_job':1725453,'physics_rerun':False,'physical_source_unchanged':True,'analysis_environment':str(Path(sys.executable).parent.parent),'environment_sha256':hashlib.sha256(freeze.encode()).hexdigest(),'checks':report['checks'],'validated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(root/'RECOVERY_PASS.json').write_text(json.dumps(receipt,indent=2)+'\n')
# Do not bypass any gate: only release after both independent checks have passed.
(root/'PASS.toml').write_bytes((root/'GPU_PASS.toml').read_bytes())
print('ALL_POINT4_GATES_PASSED',json.dumps(receipt),flush=True)
