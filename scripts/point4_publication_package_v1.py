from pathlib import Path
import json,zipfile,shutil,hashlib
root=Path.cwd();out=root/'rendered/fplane-analysis-initial-v1';receipt=json.loads((out/'PASS.json').read_text());r=json.loads((out/'analysis.json').read_text())
with zipfile.ZipFile(out/'full_diagnostic_reports.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for n in r['completed_cases']:
  src=root/'rendered/fplane-point4-v2'/n
  z.writestr(n+'/diagnostics.json',(src/'diagnostics.json').read_text().replace(str(root)+'/',''))
  shutil.copy2(src/'vorticity_PDF.png',out/(n+'_PDF.png'))
for p in out.iterdir():
 if p.is_file() and p.name!='PASS.json':
  receipt['assets'][p.name]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
assert all(v['bytes']<24000000 for v in receipt['assets'].values())
(out/'PASS.json').write_text(json.dumps(receipt,indent=2));print('Ready',len(receipt['assets'])+1,'files; diagnostic archive bytes:',(out/'full_diagnostic_reports.zip').stat().st_size)
