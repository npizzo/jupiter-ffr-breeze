"""Render completed audit aggregates on Unity."""
import json,pathlib,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype']='none'
p=pathlib.Path('data/beta-shear-deep-audit');r=json.loads((p/'audit-summary.json').read_text())
fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
labels={'coarse':'Original y/z grid','fine':'2× finer y/z grid','traditional':'Horizontal rotation removed'}
for name,v in r['paired_tests'].items():
    s=v['samples'];a=[x['velocity_difference_rms']/s[0]['velocity_difference_rms'] for x in s]
    axs[0,0].plot([x['time'] for x in s],a,label=labels[name])
axs[0,0].set(title='Paired perturbations: growth survives refinement',xlabel='Nondimensional time',ylabel='Velocity difference / initial difference');axs[0,0].legend(fontsize=8)
s=r['paired_tests']['coarse']['samples'];axs[0,1].plot([x['time'] for x in s],[100*x['doubling_relative_error'] for x in s]);axs[0,1].set(title='Amplitude check: not perfectly linear',xlabel='Nondimensional time',ylabel='Departure from doubled response (%)');axs[0,1].set_ylim(bottom=0)
e=sorted([a for a in r['energy_audits'] if a['case']=='original_3d'],key=lambda a:a['time']);x=np.arange(len(e));w=.34
axs[1,0].bar(x-w/2,[a['measured_energy_derivative'] for a in e],w,label='Recorded dE/dt');axs[1,0].bar(x+w/2,[a['predicted_energy_derivative'] for a in e],w,label='Native-operator audit');axs[1,0].set_xticks(x,[str(int(a['time'])) for a in e]);axs[1,0].set(title='Original 3D pilot: energy loss accounted for',xlabel='Nondimensional time',ylabel='Total-energy rate');axs[1,0].legend(fontsize=8)
f=r['forensic'];h=f['health_last_rows'];axs[1,1].plot([x['time'] for x in h],[x['min_rhotheta'] for x in h],'o-');axs[1,1].axhline(0,color='black',ls=':');axs[1,1].ticklabel_format(useOffset=False,style='plain',axis='x');axs[1,1].set(title='Replay: thermodynamic variable loses positivity',xlabel='Nondimensional time',ylabel='Minimum density × potential temperature');axs[1,1].text(.04,.96,'First invalid cell: top layer\nDensity and velocities remain finite',transform=axs[1,1].transAxes,va='top',fontsize=9)
fig.suptitle('Shear-pilot diagnostic audit — rejected run, not an LSV result',fontsize=15)
for ax in axs.flat:ax.grid(alpha=.15)
fig.savefig(p/'audit-checks.svg');fig.savefig(p/'audit-checks.png',dpi=130);plt.close(fig)
print('Rendered audit-checks.svg and audit-checks.png')
