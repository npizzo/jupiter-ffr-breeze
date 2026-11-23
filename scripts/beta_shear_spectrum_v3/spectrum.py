"""Unity-only inviscid, adiabatic x-uniform compressible Euler spectrum.

This is a separate centered finite-volume linearization, not the WENO/SGS
time-stepping Jacobian. Source terms, diffusion, and mean-state evolution are
omitted explicitly. Report eigenvector residuals and grid sensitivity. Positive
PV alone and a sparse eigenvalue search do not certify full stability.
"""
from pathlib import Path
import sys,json,time
import numpy as np
import h5py
import scipy.sparse as sp
from scipy.sparse.linalg import eigs,ArpackNoConvergence

def operators(ny,nz,dy,dz):
    def one(n,d):
        rows=np.r_[np.arange(n-1),np.arange(1,n)]
        cols=np.r_[np.arange(n-1),np.arange(n-1)]
        D=sp.coo_matrix((np.r_[np.ones(n-1),-np.ones(n-1)]/d,(rows,cols)),shape=(n,n-1)).tocsr()
        A=sp.coo_matrix((np.ones(2*(n-1))*.5,(rows,cols)),shape=(n,n-1)).tocsr()
        return D,A
    d,a=one(ny,dy);Dy=sp.kron(sp.eye(nz),d,format='csr');Ay=sp.kron(sp.eye(nz),a,format='csr')
    d,a=one(nz,dz);Dz=sp.kron(d,sp.eye(ny),format='csr');Az=sp.kron(a,sp.eye(ny),format='csr')
    return Dy,Dz,Ay,Az

def matrix(r,T,U,dy,dz,f,h,beta,g=12.):
    nz,ny=r.shape;n=ny*nz;Dy,Dz,Ay,Az=operators(ny,nz,dy,dz);By=Ay.T;Bz=Az.T
    theta=((r*T)**.6/r).ravel();c=(5/3*(r*T)**.4).ravel();u=U.ravel()
    ff=np.tile(f+beta*(np.arange(1,ny)*dy-ny*dy/2),nz)
    I=sp.eye(n,format='csr');Z=sp.csr_matrix((n,n));Zy=sp.csr_matrix((Dy.shape[1],n));Zz=sp.csr_matrix((Dz.shape[1],n))
    # (delta rho, delta rho*theta, delta rho*u, delta rho*v, delta rho*w).
    return sp.bmat([
        [Z,Z,Z,-Dy,-Dz],
        [Z,Z,Z,-Dy@sp.diags(By@theta),-Dz@sp.diags(Bz@theta)],
        [Z,Z,Z,-Dy@sp.diags(By@u)+Ay@sp.diags(ff),-Dz@sp.diags(Bz@u)-h*Az],
        [Zy,Dy.T@sp.diags(c),-sp.diags(ff)@By,None,None],
        [-g*Bz,Dz.T@sp.diags(c),h*Bz,None,None]],format='csc')

def unit_checks():
    # Homogeneous compressible rotating fluid, zero gravity, closed normal walls:
    # the complete nonzero spectrum must be imaginary (acoustic/inertial waves).
    r=np.ones((6,8));J=matrix(r,r,np.zeros_like(r),.5,.3,5.54,2.30,.0032,g=0.)
    w=np.linalg.eigvals(J.toarray());err=float(max(abs(w.real)))
    assert err<1e-6,err
    # Constant mass is an exact left nullvector: wall-normal mass flux vanishes.
    mass=np.r_[np.ones(r.size),np.zeros(J.shape[0]-r.size)]
    conservation=float(np.max(abs(mass@J)));assert conservation<1e-12
    return {'uniform_wave_max_abs_real_eigenvalue':err,'mass_left_null_error':conservation}

out=Path(sys.argv[1]);out.mkdir(exist_ok=True)
source=Path(sys.argv[2]);result={"scope":"Inviscid adiabatic x-uniform compressible spectrum. Equivalent second-order cross-plane operator removes the stationary balanced nullspace. Not the forced WENO/SGS Jacobian.","unit_checks":unit_checks(),"cases":{}}
for p in sorted(source.glob("*.jld2")):
    t0=time.monotonic()
    with h5py.File(p) as q:a={k:q[k][()] for k in ["rho","T","u","dy","dz","f","h","beta"]}
    n=a["rho"].size;J=matrix(a.pop("rho"),a.pop("T"),a.pop("u"),**a)
    A=J[:3*n,3*n:];B=J[3*n:,:3*n];K=(B@A).tocsc();found=[];errors=[]
    for shift in [.25,1.]:
        try:vals,vecs=eigs(K,k=6,sigma=shift,which="LM",tol=1e-10,maxiter=500)
        except ArpackNoConvergence as e:vals,vecs=e.eigenvalues,e.eigenvectors;errors.append(str(e))
        for k,mu in enumerate(vals):
            lam=np.sqrt(complex(mu));y=vecs[:,k];z=np.r_[A@y/lam,y]
            res=np.linalg.norm(J@z-lam*z)/max(np.linalg.norm(J@z),abs(lam)*np.linalg.norm(z),1e-15)
            found.append({"real":float(lam.real),"imag":float(lam.imag),"relative_residual":float(res),"squared_shift":shift})
    result["cases"][p.stem]={"eigenvalues":found,"errors":errors,"seconds":time.monotonic()-t0,"size":J.shape[0]}
    (out/"spectrum.json").write_text(json.dumps(result,indent=2)+"\n")
    print(p.stem,"max_found_real",max(x["real"] for x in found),"maxres",max(x["relative_residual"] for x in found),"seconds",time.monotonic()-t0,flush=True)
