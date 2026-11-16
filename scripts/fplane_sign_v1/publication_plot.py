"""Bounded f-sign source audit. All analysis and rendering on Unity."""
from pathlib import Path
import sys,json,tomllib,subprocess,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio_ffmpeg
from gate import check_case,read
root,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
names=['mean transport','mean flux − transport','convective stress','native − centered advection','Coriolis','diffusion','horizontal forcing']
report={'scope':'Existing capped heated anelastic f-plane; same day-90 state, opposite f; six-hour response, not a fresh sign-selection ensemble.','compute_location':'Unity','cases':{},'criteria':{'native':'Every-step full-channel closure <=1e-8 of increment scale plus explicit roundoff allowance','attribution':'Unassigned numerical-advection and mean-flux corrections smaller than 10% of claimed signed convective moment contribution; sign stable across time blocks. Otherwise mechanism remains unresolved.'}}
fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')
for case,color in [('plus','tab:blue'),('minus','tab:orange')]:
 d=root/case;gate=check_case(d);meta=tomllib.loads((d/'metadata.toml').read_text());assert meta['status']=='COMPLETE'
 records=[read(p) for p in sorted((d/'native_vorticity').glob('budget_*.jld2'))];last=records[-1];t=np.array([float(h['t']) for h in records]);hours=(t-t[0])/3600
 raw=np.array([h['third_moment_storage'][0] for h in records]);central=np.array([h['central_third_storage'][0] for h in records]);parts=np.array([h['central_third_parts'][0] for h in records]);z=last['Z'];f=float(last['f']);duration=t[-1]-t[0]
 axs[0,0].plot(hours,raw,label=case,color=color);axs[0,1].plot(hours,central,label=case,color=color)
 x=np.arange(7);axs[0,2].bar(x+(-.18 if case=='plus' else .18),parts[-1]/duration,width=.36,label=case,color=color)
 q=last['conditional_sum'][:,:,0];count=last['conditional_count'][:,0];means=np.divide(q,count[:,None],out=np.full_like(q,np.nan),where=count[:,None]>0)
 edges=last['conditional_edges'];centers=(edges[1:]+edges[:-1])/2;valid=np.isfinite(centers)&(count>0)
 axs[1,0].plot(centers[valid],means[valid,2],'.-',label=case,color=color)
 even=.5*(means[:,2]+means[::-1,2]);valid_even=valid&np.isfinite(even)&(centers>0)
 axs[1,1].plot(centers[valid_even],even[valid_even],'.-',label=case,color=color)
 # Six consecutive one-hour blocks measure sign persistence, not independent-sample confidence.
 block=[]
 for i in range(6):
  ia=np.argmin(abs(hours-i));ib=np.argmin(abs(hours-(i+1)))
  contribution=parts[ib]-parts[ia];block.append(contribution.tolist())
 axs[1,2].plot(np.arange(6)+.5,np.array(block)[:,2]/3600,'.-',label=case,color=color)
 total=parts[-1];unassigned=abs(total[1])+abs(total[3]);den=abs(total[2]);ratio=float(unassigned/max(den,1e-40))
 sign_persistent=bool(np.all(np.sign(np.array(block)[:,2])==np.sign(total[2])) and den>0)
 report['cases'][case]={'validation':gate,'f':f,'days':[t[0]/86400,t[-1]/86400],'wall_seconds':meta['wall_seconds'],'central_third_term_changes_full':dict(zip(names,total.tolist())),'raw_third_term_changes_full':dict(zip(names,last['third_moment_parts'][0].tolist())),'regional_raw_moment_changes':last['third_moment_parts'].tolist(),'regional_central_moment_changes':last['central_third_parts'].tolist(),'one_hour_central_moment_term_changes':block,'unassigned_over_convective_central_moment':ratio,'convective_sign_persists_across_blocks':sign_persistent,'physical_attribution_screen_pass':ratio<.1 and sign_persistent,'conditional_edges':[float(v) if np.isfinite(v) else ('-inf' if v<0 else 'inf') for v in edges],'conditional_count':count.tolist(),'conditional_drift_terms':[[float(v) if np.isfinite(v) else None for v in row] for row in means],'interpretation':'Differences include transient adjustment of an inherited +f state. Symmetry check uses a separate reflected state; this physical pair is not expected to be pointwise mirrored.'}
axs[0,0].set(title='Raw third-moment change, full channel',xlabel='Hours after restart',ylabel='Change in mean(Z^3) [s^-3]')
axs[0,1].set(title='Central third-moment change, full channel',xlabel='Hours after restart',ylabel='Change in mean((Z - mean(Z))^3) [s^-3]')
axs[0,2].set(title='Central-moment contributions',ylabel='Mean tendency [s⁻⁴]',xticks=np.arange(7),xticklabels=names);axs[0,2].tick_params(axis='x',rotation=65)
axs[1,0].set(title='Conditional total convective-stress tendency',xlabel='Z / |f|',ylabel='E[C | Z] [s⁻²]')
axs[1,1].set(title='Even part: [C(+|Z|)+C(−|Z|)] / 2',xlabel='|Z| / |f|',ylabel='Even drift [s⁻²]')
axs[1,2].set(title='Convective contribution by time block',xlabel='Hours after restart',ylabel='Central third-moment tendency [s⁻⁴]')
for ax in axs.flat:ax.legend();ax.axhline(0,color='.5',lw=.5)
fig.savefig(out/'sign_budget.png',dpi=160);plt.close(fig)

import hashlib
p=out/"PASS.json"
m=json.loads(p.read_text())
f=out/"sign_budget.png"
m["assets"][f.name]={"bytes":f.stat().st_size,"sha256":hashlib.sha256(f.read_bytes()).hexdigest()}
m["publication_note"]="Budget figure axis labels redrawn with supported characters; all numerical values and movie unchanged."
p.write_text(json.dumps(m,indent=2)+"\n")
