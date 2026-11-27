"""Run only on Unity. Inspect published vector artwork; no guessed profile fits."""
import sys,json,hashlib
from pathlib import Path
import fitz
root=Path.cwd();pdf=root/'references/zhang-marcus-2024/paper.pdf';out=root/'runs/marcus-cv1-v1/reference';out.mkdir(parents=True,exist_ok=True)
doc=fitz.open(pdf);result={'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'pages':{}}
for n in [9,11,14,15,16,17,18,38]:
 p=doc[n];paths=[]
 for d in p.get_drawings():
  col=d['color']
  if col and max(col)-min(col)>.15:
   paths.append(dict(color=col,rect=list(d['rect']),n=len(d['items']),items=[(x[0],*[list(v) if hasattr(v,'x') else v for v in x[1:]]) for x in d['items']]))
 result['pages'][str(n+1)]={'words':p.get_text('words'),'colored_paths':paths}
(out/'vector_artwork.json').write_text(json.dumps(result,default=str))
for n in [9,11]:
 p=doc[n];p.get_pixmap(matrix=fitz.Matrix(2,2)).save(out/f'page{n+1}.png')
print(json.dumps({'pdf_sha256':result['pdf_sha256'],'pages':{k:[(p['color'],p['rect'],p['n']) for p in v['colored_paths']] for k,v in result['pages'].items()}}))
