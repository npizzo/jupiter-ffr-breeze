# These are resolved, cell-centered diagnostics, NOT native WENO budget terms.
sx(a,dx)=(circshift(a,(-1,0,0)).-circshift(a,(1,0,0)))./(2dx)
function sy(a,dy)
 p=cat(a[:,1:1,:],a,a[:,end:end,:];dims=2)
 return (p[:,3:end,:].-p[:,1:end-2,:])./(2dy)
end
function sz(a,dz)
 n=size(a,3)
 return cat((a[:,:,2:2].-a[:,:,1:1])./dz,
    (a[:,:,3:n].-a[:,:,1:n-2])./(2dz),
    (a[:,:,n:n].-a[:,:,n-1:n-1])./dz;dims=3)
end
function layer_diagnostics(u,v,w,zeta,theta,rho,z,dx,dy,dz,f,zi)
 div=sx(u,dx).+sy(v,dy)
 stretch=.-(f.+zeta).*div
 tilt=sz(u,dz).*sy(w,dy).-sz(v,dz).*sx(w,dx)
 zp=zeta.-mean(zeta,dims=1)
 q=reshape(rho,1,1,:)
 n=size(u,2); core=ceil(Int,.1n):floor(Int,.9n)
 profile(a)=vec(mean(a;dims=(1,2)))
 sigma=sqrt.(profile(zp.^2))
 safe_skew(a)=mean(a.^3)/max(mean(a.^2)^1.5,eps(Float64)^3)
 layers=Dict{String,Any}()
 for (label,mask) in ["lower"=>(10000 .<= z .<= 70000),
                       "interface"=>(75000 .<= z .<= 90000),
                       "upper"=>(90000 .<= z .<= 105000),"full"=>trues(length(z))]
  any(mask)||error("Empty diagnostic layer")
  weights=reshape(rho[mask]./sum(rho[mask]),1,1,:)
  avg(a)=dropdims(sum(a[:,:,mask].*weights;dims=3);dims=3)
  zz=avg(zeta); pp=zz.-mean(zz;dims=1)
  uu=avg(u); vv=avg(v)
  layers[label]=Dict("z_min_m"=>minimum(z[mask]),"z_max_m"=>maximum(z[mask]),
     "u"=>Float32.(uu),"v"=>Float32.(vv),"zeta"=>Float32.(zz),"zeta_prime"=>Float32.(pp),
     "divergence"=>Float32.(avg(div)),"vertical_velocity"=>Float32.(avg(w)),
     "stretching_centered"=>Float32.(avg(stretch)),"tilting_centered"=>Float32.(avg(tilt)),
     "signed_skew_prime"=>sign(f)*safe_skew(pp),"signed_skew_prime_interior"=>sign(f)*safe_skew(pp[:,core]),
     "mean_zeta_over_f"=>mean(zz)/f,"rms_prime_over_abs_f"=>sqrt(mean(pp.^2))/abs(f),
     "zonal_u"=>vec(mean(uu;dims=1)),"zonal_vorticity_over_f"=>vec(mean(zz;dims=1))./f)
 end
 return Dict("layers"=>layers,"z_m"=>z,
    "signed_skew_prime_profile"=>sign(f).*profile(zp.^3)./max.(sigma.^3,eps(Float64)^3),
    "rms_prime_over_abs_f_profile"=>sigma./abs(f),
    "upflow_divergence_covariance_profile"=>profile(w.*div),
    "stretching_prime_M3_tendency_profile"=>3sign(f).*profile(zp.^2 .*(stretch.-mean(stretch;dims=1))),
    "tilting_prime_M3_tendency_profile"=>3sign(f).*profile(zp.^2 .*(tilt.-mean(tilt;dims=1))),
    "diagnostic_caveat"=>"Centered stretching/tilting are correlations, not a closed native tendency budget; zonal circulation and walls retained separately")
end
