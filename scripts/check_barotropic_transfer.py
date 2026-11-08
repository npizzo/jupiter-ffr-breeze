import sys,numpy as np
sys.path.insert(0,str(__import__('pathlib').Path(__file__).parent))
from barotropic_transfer import transfer
n=64;x=np.arange(n)*2*np.pi/n;y=(np.arange(n)+.5)*2*np.pi/n;xx,yy=np.meshgrid(x,y)
u=np.sin(xx)*np.cos(yy);v=-np.cos(xx)*np.sin(yy);zero=np.zeros_like(u)
r=transfer(dict(ubar=u,vbar=v,Rxx=zero,Rxy=zero,Ryy=zero),2*np.pi,2*np.pi)
assert max(abs(r['self_flux']))<1e-14 and r['projection_energy_fraction']<1e-25
# Multiple interacting modes verify conservative self-transfer and reflection.
u=np.sin(xx)*np.cos(yy)+.3*np.cos(2*xx)*np.cos(3*yy)
v=-np.cos(xx)*np.sin(yy)+.2*np.sin(2*xx)*np.sin(3*yy)
r=transfer(dict(ubar=u,vbar=v,Rxx=zero,Rxy=zero,Ryy=zero),2*np.pi,2*np.pi)
assert abs(r['self_transfer_residual'])<1e-14
# Channel mode has v=0 at y=0,pi; centers avoid exact-wall samples.
y=(np.arange(n)+.5)*np.pi/n;xx,yy=np.meshgrid(x,y)
r=transfer(dict(ubar=np.sin(xx)*np.cos(yy),vbar=-np.cos(xx)*np.sin(yy),Rxx=zero,Rxy=zero,Ryy=zero),2*np.pi,np.pi,True)
assert max(abs(r['self_flux']))<1e-14 and r['projection_energy_fraction']<1e-25
print('Spectral diagnostic checks passed: analytic solenoidal modes, energy Parseval identity, conservative self-transfer and channel parity.')

# A closed non-collinear triad has nonzero shell transfer with zero total work.
y=(np.arange(n)+.5)*2*np.pi/n;xx,yy=np.meshgrid(x,y)
u=np.zeros_like(xx);v=np.zeros_like(xx)
for kx,ky,amplitude,phase in [(1,1,1,.2),(2,3,.3,.4),(3,4,.2,1.3)]:
    wave=np.sin(kx*xx+ky*yy+phase)
    u-=amplitude*ky*wave;v+=amplitude*kx*wave
r=transfer(dict(ubar=u,vbar=v,Rxx=zero,Rxy=zero,Ryy=zero),2*np.pi,2*np.pi)
assert max(abs(r['self_transfer']))>1e-6
assert abs(r['self_transfer_residual'])<1e-12
print('Nonzero triad transfer and zero total nonlinear work verified.')
