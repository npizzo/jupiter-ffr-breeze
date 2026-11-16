# Compatible native-grid modal kinetic energy; q is corner curl.
# FFTs are unnormalised: Parseval factor is (nx*ny)^(-2).
function spectral_budget(o,q,parts)
    qh=FFTW.fft(q);norm=(o.nx*o.ny)^2
    invlam=map(x->x==0 ? 0.0 : 1/x,o.lambda)
    emode=abs2.(qh).*invlam./(2norm)
    transfers=[real.(conj.(qh).*FFTW.fft(p)).*invlam./norm for p in parts]
    # Integer upper shell edges: shell K contains K-1 < |k| <= K.
    shell=ceil.(Int,sqrt.(o.k2).-1e-12);nk=maximum(shell)
    E=zeros(nk+1);T=zeros(nk+1,length(parts))
    for I in eachindex(q)
        k=shell[I]+1;E[k]+=emode[I]
        for j in eachindex(parts);T[k,j]+=transfers[j][I];end
    end
    flux=-cumsum(T[:,1])
    (;E,T,flux,k=collect(0:nk))
end

function write_spectra(path,clock,s)
    open(path,"a") do io
        for j in eachindex(s.k)
            println(io,join((clock.iteration,clock.time,0.6*clock.time,s.k[j],s.E[j],s.T[j,:]...,s.flux[j]),','))
        end
    end
end

function spectral_checks(o,q,parts,u,v)
    s=spectral_budget(o,q,parts);en=mean(u.^2 .+v.^2)/2
    eerr=abs(sum(s.E)-en)/max(en,1e-30)
    psi=spectral_curl(o,q,map(x->x==0 ? 0.0 : -1/x,o.lambda))
    work=[-mean(psi.*p) for p in parts]
    terr=maximum(abs.(vec(sum(s.T;dims=1)).-work))/max(maximum(abs,work),1e-30)
    # Native incompressible centered advection should conserve total energy.
    net=abs(sum(s.T[:,1]))/max(sum(abs,s.T[:,1]),1e-30)
    eerr<1e-9||error("Spectral/face energy mismatch: $eerr")
    terr<1e-9||error("Parseval work mismatch: $terr")
    # Near-zero transfers need an absolute tolerance as well.
    (net<1e-8 || abs(sum(s.T[:,1]))<1e-14)||error("Native advection has nonconservative net transfer")
    (;s,eerr,terr,net)
end
