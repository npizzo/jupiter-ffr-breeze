# CPU-only state preparation on Unity. Scalar cell integrals and periodic
# horizontal-momentum dual-volume integrals are preserved. Vertical face
# momentum is linearly interpolated, with both impermeable faces exactly zero.
@inline minmod(a,b)=a*b<=0 ? zero(a) : copysign(min(abs(a),abs(b)),a)
function remap_axis(A,dim,nnew;periodic=true,offset=.5,vertical_face=false)
    perm=Tuple(vcat(dim,[j for j=1:3 if j!=dim]));ip=invperm(collect(perm))
    B=permutedims(A,perm);ns=size(B);B=reshape(B,ns[1],:)
    nold=vertical_face ? ns[1]-1 : ns[1]
    C=zeros(Float64,vertical_face ? nnew+1 : nnew,size(B,2))
    if vertical_face
        for j=0:nnew
            x=j*nold/nnew;k=min(floor(Int,x),nold-1);q=x-k
            @views C[j+1,:].=(1-q).*B[k+1,:].+q.*B[k+2,:]
        end
        C[1,:].=0;C[end,:].=0
    else
        h=nold/nnew
        for col=axes(B,2),j=0:nnew-1
            a=(j+offset-.5)*h;b=a+h;total=0.
            # Old cell center k+offset; its dual cell starts at k+offset-.5.
            k0=floor(Int,a-offset+.5+1e-11)
            k1=ceil(Int,b-offset+.5-1e-11)-1
            for k=k0:k1
                lo=max(a,k+offset-.5);hi=min(b,k+offset+.5)
                ii=periodic ? mod(k,nold)+1 : clamp(k+1,1,nold)
                il=periodic ? mod(k-1,nold)+1 : max(ii-1,1)
                ir=periodic ? mod(k+1,nold)+1 : min(ii+1,nold)
                v=B[ii,col];s=minmod(v-B[il,col],B[ir,col]-v)
                total+=v*(hi-lo)+.5s*((hi-k-offset)^2-(lo-k-offset)^2)
            end
            C[j+1,col]=total/h
        end
    end
    permutedims(reshape(C,size(C,1),ns[2],ns[3]),ip)
end
function transfer(A,target;face=0)
    B=A
    for d=1:3
        B=remap_axis(B,d,target[d];periodic=d<3,offset=d==face ? 0. : .5,vertical_face=(d==3 && face==3))
    end
    B
end
function extract_state(m)
    (;rho=Array(interior(density_field(m))),rhotheta=Array(interior(m.formulation.potential_temperature_density)),
      mx=Array(interior(m.momentum.ρu)),my=Array(interior(m.momentum.ρv)),mz=Array(interior(m.momentum.ρw)),
      dye=Array(interior(m.tracers.dye)),t=Float64(m.clock.time))
end
function install_state!(m,s)
    set!(m;ρ=s.rho,ρθ=s.rhotheta,ρu=s.mx,ρv=s.my,ρw=s.mz,qᵗ=0.,dye=s.dye)
    m.clock.time=s.t
    Oceananigans.TimeSteppers.update_state!(m)
end
function mapped_state(s,target)
    (;rho=transfer(s.rho,target),rhotheta=transfer(s.rhotheta,target),
      mx=transfer(s.mx,target;face=1),my=transfer(s.my,target;face=2),mz=transfer(s.mz,target;face=3),
      dye=transfer(s.dye,target),t=s.t)
end
function state_error(a,b)
    maximum(maximum(abs,getproperty(a,k).-getproperty(b,k)) for k in (:rho,:rhotheta,:mx,:my,:mz,:dye))
end
