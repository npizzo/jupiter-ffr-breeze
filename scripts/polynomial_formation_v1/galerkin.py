"""A selected-mode Galerkin projection of the SAME 512-grid vorticity dynamics.

All coefficients are computed from the native centered-flux Fourier symbols;
none are fitted. Selection uses t<=200; t>=210 is held out. Unity only.
"""
from pathlib import Path
import sys,json,time,signal
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'feedback_maintenance_v1'))
import analyze_existing as ref

def bins(indices, values, size):
    return np.bincount(indices,weights=values.real,minlength=size)+1j*np.bincount(indices,weights=values.imag,minlength=size)

class Projection:
    def __init__(self,modes,B=.002,C=.0002):
        self.modes=np.array(sorted(set(map(tuple,modes))),int);self.M=len(self.modes);self.B=B;self.C=C
        assert all(tuple(-v) in set(map(tuple,self.modes)) for v in self.modes)
        assert not np.any(np.all(self.modes==0,axis=1))
        x,y=self.modes.T;d=2*np.pi/512
        self.lam=(4*np.sin(x*d/2)**2+4*np.sin(y*d/2)**2)/d**2
        self.k2=x*x+y*y;self.linear=.0002*self.k2*((self.k2>=144)&(self.k2<=256)&(x!=0))-.0001*self.lam
        self.low=self.k2<=4
        self.p=np.repeat(np.arange(self.M),self.M);self.q=np.tile(np.arange(self.M),self.M)
        sums=self.modes[self.p]+self.modes[self.q]
        self.sums,self.pair_output=np.unique(sums,axis=0,return_inverse=True);self.J=len(self.sums)
        index={tuple(v):i for i,v in enumerate(self.modes)};sindex={tuple(v):i for i,v in enumerate(self.sums)}
        self.retained=np.array([index.get(tuple(v),-1) for v in sums]);self.valid=self.retained>=0
        self.q2_on_modes=np.array([sindex.get(tuple(v),-1) for v in self.modes])
        self.zero=sindex.get((0,0),None)
        out=[];mid=[];last=[]
        for i,k in enumerate(self.modes):
            for j,c in enumerate(self.modes):
                h=sindex.get(tuple(k-c))
                if h is not None:out.append(i);mid.append(h);last.append(j)
        self.cout=np.array(out);self.cmid=np.array(mid);self.clast=np.array(last)
        U=(np.exp(1j*y*d)-1)/d/self.lam;V=-(np.exp(1j*x*d)-1)/d/self.lam
        Ux=.5*(1+np.exp(1j*x*d))*U;Vy=.5*(1+np.exp(1j*y*d))*V
        Cy=.5*(1+np.exp(-1j*y*d))*U;Cx=.5*(1+np.exp(-1j*x*d))*V
        sx,sy=sums.T;dxm=(1-np.exp(-1j*sx*d))/d;dym=(1-np.exp(-1j*sy*d))/d
        dxp=(np.exp(1j*sx*d)-1)/d;dyp=(np.exp(1j*sy*d)-1)/d
        cross=Cy[self.p]*Cx[self.q]
        au=-dxm*Ux[self.p]*Ux[self.q]-dyp*cross
        av=-dxp*cross-dym*Vy[self.p]*Vy[self.q]
        self.A=dxm*av-dym*au

    def terms(self,z):
        product=z[self.p]*z[self.q]
        q2=bins(self.pair_output,product,self.J)
        quad=bins(self.retained[self.valid],product[self.valid],self.M)*self.B
        adv=bins(self.retained[self.valid],(self.A*product)[self.valid],self.M)
        cubic=-self.C*bins(self.cout,q2[self.cmid]*z[self.clast],self.M)
        return {'advection':adv,'quadratic':quad,'cubic':cubic,'linear_and_viscous':self.linear*z}

    def rhs(self,t,z):return sum(self.terms(z).values())

    def extract(self,z):return z[self.modes[:,0]%512,self.modes[:,1]%512]

def reference_terms(q):
    z=np.fft.fft2(q)/512**2;u,v=ref.velocity(q);au,av=ref.centered_advection(u,v)
    return z,{'advection':np.fft.fft2(ref.curl(au,av))/512**2,
              'quadratic':.002*np.fft.fft2(q*q)/512**2,'cubic':-.0002*np.fft.fft2(q*q*q)/512**2,
              'linear_and_viscous':(.0002*ref.K2*ref.SHELL-.0001*ref.LAM)*z}

def validation():
    modes=[(1,0),(-1,0),(0,1),(0,-1),(12,1),(-12,-1),(13,1),(-13,-1),(2,0),(-2,0)]
    p=Projection(modes);rng=np.random.default_rng(42);grid=np.zeros((512,512),complex)
    for x,y in modes:
        if x>0 or (x==0 and y>0):
            v=rng.normal()+1j*rng.normal();grid[x%512,y%512]=v;grid[-x%512,-y%512]=v.conjugate()
    q=np.fft.ifft2(grid*512**2).real;z=p.extract(grid);_,full=reference_terms(q);reduced=p.terms(z)
    errors={k:float(np.linalg.norm(v-p.extract(full[k]))/max(np.linalg.norm(p.extract(full[k])),1e-12)) for k,v in reduced.items()}
    assert max(errors.values())<1e-9,errors
    work=float(np.sum((z.conj()*reduced['advection']).real/p.lam));assert abs(work)<1e-10
    return {'source_projection_relative_errors':errors,'advective_energy_work':work,'status':'PASS'}

def main(root,analysis,out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    checks=validation();ranking=json.loads((analysis/'mode_ranking.json').read_text())['pairs']
    frames,_,_,_=ref.inventory(root,'biased_plus');mandatory=[(x,y) for x in range(-2,3) for y in range(-2,3) if 0<x*x+y*y<=4]
    candidates=[];selected=None
    # Three bounded sizes; no automatic expansion if they fail the held-out test.
    for pairs in (16,32,64):
        selected_pairs=[]
        for band,number in ((True,pairs//2),(False,pairs//2)):
            eligible=[v for v in ranking if v['forcing_band']==band and v['kx']**2+v['ky']**2>4]
            selected_pairs.extend(eligible[:number])
        modes=mandatory.copy()
        for v in selected_pairs:modes.extend([(v['kx'],v['ky']),(-v['kx'],-v['ky'])])
        p=Projection(modes);errors={k:0. for k in ('total','quadratic','cubic','advection')};norms=errors.copy()
        energy_numerator=energy_denominator=0.
        for t in np.arange(210,401,10):
            q=ref.load_field(frames[t]);z,full=reference_terms(q);zz=p.extract(z);reduced=p.terms(zz)
            total_full=sum(full.values());total_reduced=sum(reduced.values())
            for k in errors:
                a=total_reduced if k=='total' else reduced[k];b=p.extract(total_full if k=='total' else full[k])
                errors[k]+=float(np.sum(np.abs((a-b)[p.low])**2/p.lam[p.low]))
                norms[k]+=float(np.sum(np.abs(b[p.low])**2/p.lam[p.low]))
            energy_numerator+=float(np.sum(np.abs(zz)**2/p.lam))
            energy_denominator+=float(np.sum(np.abs(z)**2*ref.INVLAM))
        errors={k:float(np.sqrt(errors[k]/max(norms[k],1e-30))) for k in errors}
        passed=all(errors[k]<=.2 for k in ('total','quadratic','cubic'))
        entry={'complex_modes_including_conjugates':p.M,'independent_real_dimensions':p.M,'selected_pairs_target':pairs,
               'modes':p.modes.tolist(),'held_out_relative_lowmode_tendency_errors':errors,
               'held_out_energy_fraction_retained':energy_numerator/energy_denominator,
               'prediction_gate_pass':passed,'gate':'relative low-mode total, quadratic and cubic RHS errors <=0.2 over withheld t=210..400'}
        candidates.append(entry);print('ROM_SCREEN',json.dumps(entry),flush=True)
        if passed and selected is None:selected=p
    forecast=None
    if selected is not None:
        # Only the smallest passing model is advanced, using the original projected IC.
        q0=ref.load_field(frames[0.]);z0=selected.extract(np.fft.fft2(q0)/512**2)
        solution=solve_ivp(selected.rhs,(0,400),z0,method='DOP853',rtol=1e-7,atol=1e-10,max_step=.5,t_eval=np.arange(0,401,5))
        assert solution.success and np.isfinite(solution.y).all()
        energies=.5*np.sum(np.abs(solution.y)**2/selected.lam[:,None],axis=0)
        low=.5*np.sum(np.abs(solution.y[selected.low])**2/selected.lam[selected.low,None],axis=0)
        np.savez(out/'reduced_trajectory.npz',times=solution.t,modes=selected.modes,zeta=solution.y,energy=energies,low_energy=low)
        # Dynamics prediction is assessed separately from instantaneous source projection.
        truth=np.array([.5*np.sum(np.abs(np.fft.fft2(ref.load_field(frames[t]))[ref.LOW])**2*ref.INVLAM[ref.LOW])/ref.NORM for t in solution.t])
        error=float(np.linalg.norm(low-truth)/max(np.linalg.norm(truth),1e-30))
        forecast={'time_range':[0,400],'rhs_evaluations':solution.nfev,'relative_low_energy_history_error':error,
                  'prediction_pass':error<=.2,'final_low_energy':float(low[-1]),'reference_final_low_energy':float(truth[-1])}
    report={'status':'COMPLETE','compute_location':'Unity','validation':checks,'selection_times':[0,200],'held_out_times':[210,400],
            'candidates':candidates,'forecast':forecast,'elapsed_seconds':time.monotonic()-started,
            'conclusion':'A reduction is not accepted merely because it produces a vortex. No coefficients were tuned to held-out data.',
            'limitations':['This is a Galerkin projection of the native 512-grid discretization, not a new convection reduction.',
                'Candidate maximum |k|=64; interactions do not alias at the original 512 grid. Omitted modes remain a closure error.',
                'Time holdout on one trajectory is weaker than independent-initial-condition validation.',
                'If no candidate passes, no reduced trajectory is integrated and no larger campaign is automatically launched.']}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False));print('ROM_AUDIT_COMPLETE',flush=True)

if __name__=='__main__':
    signal.signal(signal.SIGALRM,lambda *args:(_ for _ in ()).throw(TimeoutError('Bounded ROM audit CPU/wall budget exceeded')))
    signal.alarm(780)
    main(*map(Path,sys.argv[1:]))
