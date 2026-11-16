import numpy as np
from operators import dx,dy,curl,radial,sample,corner_curl

def checks():
    n=256;h=1.;x=(np.arange(n)+.5)*h;y=x.copy();X,Y=np.meshgrid(x,y,indexing='ij')
    # Smooth periodic derivative, wall parity and uniform flow.
    a=np.sin(2*np.pi*X/n);exact=2*np.pi/n*np.cos(2*np.pi*X/n)
    assert np.max(abs(dx(a,h)-exact))<3e-6
    assert np.max(abs(curl(np.ones((n,n)),np.ones((n,n)),h,h)))==0
    assert np.allclose(sample(a,X,Y,h,h),a)
    # Analytic solid-body discriminant (f+2Omega)^2, including a strong anticyclone.
    ij=(128,128);xx=X-X[ij];yy=Y-Y[ij];f=1.
    for om in [-1.2,-.2,.2]:
        u=-om*yy;v=om*xx;q=np.full_like(u,2*om)
        p=radial(u,v,q,ij,f,h,h)
        assert np.max(abs(p['Phi_over_f2']-(f+2*om)**2))<1e-12
        assert not p['axisymmetric_screen_pass']
    # Localized Gaussian anticyclone with a centrifugally unstable annulus.
    # V=Omega*r*exp(-r²/a²), zeta=2Omega*(1-r²/a²)*exp(-r²/a²).
    om=-1.2;radius=7.;r2=xx**2+yy**2;e=np.exp(-r2/radius**2)
    u=-om*yy*e;v=om*xx*e;q=2*om*(1-r2/radius**2)*e
    p=radial(u,v,q,ij,f,h,h);r=p['r_m'];e=np.exp(-r*r/radius**2)
    exact=(f+2*om*e)*(f+2*om*(1-r*r/radius**2)*e)
    resolved=(r>=4*h)&(r<=10*h)
    gaussian_error=float(np.max(abs(np.asarray(p['Phi_over_f2'])[resolved]-exact[resolved])))
    assert gaussian_error<.04, ('resolved Gaussian discriminant',gaussian_error)
    assert p['axisymmetric_screen_pass']
    # Mirroring the local velocity and f preserves the discriminant.
    p2=radial(-u,-v,-q,ij,-f,h,h)
    assert np.allclose(p['Phi_over_f2'],p2['Phi_over_f2'])
    # Native and centered curl commute with four-corner averaging in the interior.
    rng=np.random.default_rng(481);uf=rng.normal(size=(32,32,3));vf=rng.normal(size=(32,33,3));vf[:,[0,-1]]=0
    qc=corner_curl(uf,vf,1.,1.);qc=(qc+np.roll(qc,-1,0))/2;qc=(qc[:,:-1]+qc[:,1:])/2
    uc=(uf+np.roll(uf,-1,0))/2;vc=(vf[:,:-1]+vf[:,1:])/2
    assert np.allclose(qc,curl(uc,vc,1.,1.),rtol=0,atol=1e-14)
    return {'status':'PASS','tests':['periodic derivative','uniform-flow curl','interpolation identity',
            'three analytic solid-body vortices','localized Gaussian instability annulus','rotation/sign covariance',
            'independent face-curl and centered-curl equivalence'],
            'gaussian_max_discriminant_error_resolved_radii':gaussian_error,
            'gaussian_max_discriminant_error_all_radii':float(np.max(abs(p['Phi_over_f2']-exact))),
            'validation_radius_note':'The 0.04 discriminant-error tolerance applies at r=4–10 cells, matching the diagnostic acceptance region; smaller radii remain excluded.'}
if __name__=='__main__':print(checks())
