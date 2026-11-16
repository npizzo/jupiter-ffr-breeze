"""Centered diagnostic operators and explicitly conditional axisymmetric screening."""
import numpy as np
from scipy.ndimage import minimum_filter, map_coordinates

def dx(q,h): return (np.roll(q,-1,0)-np.roll(q,1,0))/(2*h)
def dy(q,h,parity=1):
    p=np.concatenate((parity*q[:,0:1,...],q,parity*q[:,-1:,...]),axis=1)
    return (p[:,2:,...]-p[:,:-2,...])/(2*h)
def centered(uf,vf,wf):
    return (uf+np.roll(uf,-1,0))/2,(vf[:,:-1,:]+vf[:,1:,:])/2,(wf[:,:,:-1]+wf[:,:,1:])/2
def curl(u,v,hx,hy): return dx(v,hx)-dy(u,hy)
def corner_curl(uf,vf,hx,hy):
    q=np.zeros_like(vf)
    q[:,1:-1,:]=((vf-np.roll(vf,1,0))/hx)[:,1:-1,:]-np.diff(uf,axis=1)/hy
    return q
def candidates(u,v,q,f,hx,hy,ymask,sign=-1,limit=3):
    # Complex eigenvalues of the horizontal velocity gradient; not an instability criterion.
    ux,uy,vx,vy=dx(u,hx),dy(u,hy),dx(v,hx),dy(v,hy,-1)
    rotation_dominated=q*q-(ux-vy)**2-(vx+uy)**2>0
    a=sign*q/f
    extrema=(a==-minimum_filter(-a,size=9,mode=('wrap','nearest')))&(a>=0.25)&rotation_dominated
    extrema&=ymask[None,:]
    ij=np.argwhere(extrema)
    if len(ij)==0:return []
    order=np.argsort(a[ij[:,0],ij[:,1]])[::-1]
    return [tuple(map(int,ij[k])) for k in order[:limit]]
def sample(q,xx,yy,hx,hy):
    # Cell-center coordinates, explicit periodic x padding; y coordinates stay inside walls.
    nx=q.shape[0];pad=np.concatenate((q[-1:,:],q,q[:1,:]),axis=0)
    ix=(xx/hx-0.5)%nx+1;iy=yy/hy-0.5
    return map_coordinates(pad,[ix.ravel(),iy.ravel()],order=1,mode='nearest').reshape(xx.shape)
def radial(u,v,q,ij,f,hx,hy):
    i,j=ij;h=max(hx,hy);r=np.arange(1.,13.)*h
    phi=np.linspace(0,2*np.pi,96,endpoint=False)
    x=(i+.5)*hx+r[:,None]*np.cos(phi);y=(j+.5)*hy+r[:,None]*np.sin(phi)
    assert y.min()>=hy/2 and y.max()<=(u.shape[1]-.5)*hy
    us=sample(u,x,y,hx,hy);vs=sample(v,x,y,hx,hy)
    # Remove uniform translation; complete-circle mean V is unchanged by this operation.
    us-=u[i,j];vs-=v[i,j]
    vt=-us*np.sin(phi)+vs*np.cos(phi);vr=us*np.cos(phi)+vs*np.sin(phi)
    V=vt.mean(1);zeta=sample(q,x,y,hx,hy).mean(1)
    zeta_V=np.gradient(r*V,r,edge_order=2)/r
    angular=f+2*V/r
    disc=angular*(f+zeta);disc_V=angular*(f+zeta_V)
    resolved=(r>=4*h)&(r<=10*h)
    coherence=float(np.sum(V[resolved]**2)/max(np.sum(np.mean(vt[resolved]**2+vr[resolved]**2,axis=1)),1e-30))
    negative=resolved&(disc/f**2 < -0.01)&(disc_V/f**2 < -0.01)
    contiguous=bool(np.any(negative[:-1]&negative[1:]))
    return {'r_m':r,'V_m_s':V,'azimuthal_zeta_over_f':zeta/f,'circulation_zeta_over_f':zeta_V/f,
            'Phi_over_f2':disc/f**2,'Phi_circulation_over_f2':disc_V/f**2,
            'azimuthal_velocity_std_m_s':vt.std(1),'coherence':coherence,
            'resolved_negative_both':contiguous,'axisymmetric_screen_pass':coherence>=.5 and contiguous,
            'min_Phi_over_f2_resolved':float(np.min(disc[resolved]/f**2)),
            'min_Phi_circulation_over_f2_resolved':float(np.min(disc_V[resolved]/f**2))}
