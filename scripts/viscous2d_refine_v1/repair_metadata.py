"""Correct only printed RMS metadata; preserve original and all model fields.
Julia parses 2E0 and 2E1 as exponent literals, not 2*E0 and 2*E1.
The initializer's velocities and energies are unaffected by those two metadata expressions.
"""
from pathlib import Path
import sys,tomllib,re,math,json,hashlib
p=Path(sys.argv[1]);original=p.read_text();d=tomllib.loads(original)
backup=p.with_name('initialization.original.toml');assert not backup.exists();backup.write_text(original)
corrected=dict(coarse_urms=math.sqrt(2*d['coarse_initial_energy']),fine_urms=math.sqrt(2*d['fine_initial_energy']))
text=original
for key,value in corrected.items():text=re.sub(r'^'+key+r' = .*$',key+' = '+repr(value),text,flags=re.M)
p.write_text(text)
record={'status':'CORRECTED_METADATA_ONLY','reason':'Julia exponent-literal parsing in two displayed RMS expressions','original_sha256':hashlib.sha256(original.encode()).hexdigest(),'corrected_values':corrected,'model_fields_modified':False}
p.with_name('metadata_correction.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
