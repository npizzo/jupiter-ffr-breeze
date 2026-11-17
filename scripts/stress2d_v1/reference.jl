# Independent allocating reference, used only in checks and sparse diagnostics on Unity.
function reference(u,v,c;part=:all)
    n=size(u,1);h=2pi/n;o=BandOperator(n,n;kmin=1.,kmax=KCUT,allow_zonal=true)
    P(a)=real.(FFTW.ifft(FFTW.fft(a).*(o.k2.<=KCUT^2)))
    ub,vb=P(u),P(v)
    Dxp(a)=(circshift(a,(-1,0)).-a)./h;Dyp(a)=(circshift(a,(0,-1)).-a)./h
    Dxm(a)=(a.-circshift(a,(1,0)))./h;Dym(a)=(a.-circshift(a,(0,1)))./h
    z=Dxm(vb).-Dym(ub);zc=(z.+circshift(z,(-1,0)).+circshift(z,(0,-1)).+circshift(z,(-1,-1)))./4
    sx=(Dxp(ub).-Dyp(vb))./2;sy=(Dym(ub).+Dxm(vb))./2
    nu(q)=(part in (:all,:linear) ? c.nu0 : 0.).+(part in (:all,:quadratic) ? c.nu1 : 0.).*q./c.f.+(part in (:all,:cubic) ? c.nu2 : 0.).*(q./c.f).^2
    nx,ny=nu(zc),nu(z)
    xx,xy=P(-2 .*nx.*sx),P(-2 .*ny.*sy)
    fu=-Dxm(xx).-Dyp(xy);fv=-Dxp(xy).+Dym(xx)
    curl=Dxm(fv).-Dym(fu)
    work=mean(u.*fu.+v.*fv);identity=-4mean(nx.*sx.^2 .+ ny.*sy.^2)
    (;fu,fv,curl,work,identity,z,sx,sy,ub,vb)
end
function diagnostics(u,v,c)
    n=size(u,1);o=BandOperator(n,n;kmin=1.,kmax=KCUT,allow_zonal=true)
    z=vorticity(o,u,v);zh=FFTW.fft(z);invlam=map(x->x==0 ? 0. : 1/x,o.lambda)
    low=(o.k2.>0).&(o.k2.<=4);norm=n^4;m2=mean(z.^2)
    d=Dict{String,Float64}("energy"=>.5mean(u.^2 .+ v.^2),"enstrophy"=>.5m2,"mean_zeta"=>mean(z),
      "skewness"=>mean(z.^3)/m2^1.5,"third_moment"=>mean(z.^3),
      "energy_k_le_2"=>.5sum((abs2.(zh).*invlam)[low])/norm,
      "highk_enstrophy_fraction"=>sum(abs2.(zh)[o.k2.>=(n/4)^2])/sum(abs2,zh),
      "viscous_energy_rate"=>-VISCOSITY*sum(abs2,zh)/norm)
    for part in (:linear,:quadratic,:cubic)
        r=reference(u,v,c;part)
        @assert isapprox(r.work,r.identity;rtol=1e-9,atol=1e-12)
        part==:cubic && @assert r.work<=1e-12
        ch=FFTW.fft(r.curl)
        d[string(part)*"_energy_rate"]=r.work
        d[string(part)*"_energy_identity"]=r.identity
        d[string(part)*"_lowk_energy_rate"]=sum((real.(conj.(zh).*ch).*invlam)[low])/norm
        d[string(part)*"_m3_rate"]=3mean(z.^2 .*r.curl)
        if part==:linear
            d["filtered_rossby_rms"]=sqrt(mean(abs2,r.z))/abs(c.f)
            d["filtered_rossby_max"]=maximum(abs,r.z)/abs(c.f)
            d["filtered_fraction_abs_ro_gt_02"]=mean(abs.(r.z./c.f).>.2)
            d["filtered_strain_max_over_f"]=maximum(sqrt.(2r.sx.^2 .+ 2r.sy.^2))/abs(c.f)
        end
    end
    d
end
