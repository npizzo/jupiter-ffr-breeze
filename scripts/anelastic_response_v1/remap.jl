# Conservative horizontal refinement of staggered mass fluxes on Unity.
# D_f P_F = P_C D_c, so the same linear center prolongation in transverse
# directions preserves discrete anelastic continuity without a new projection.
function prolong_center(a,dim;periodic=false,limited=false)
    n=size(a,dim);sz=collect(size(a));sz[dim]=2n;b=Array{eltype(a)}(undef,sz...)
    slab(j)=selectdim(a,dim,periodic ? mod1(j,n) : clamp(j,1,n))
    for j=1:n
        q=slab(j);qm=slab(j-1);qp=slab(j+1)
        delta = if limited
            left=q.-qm;right=qp.-q
            ifelse.(left.*right.>0,sign.(left).*min.(abs.(left),abs.(right))./4,zero(eltype(a)))
        else
            (qp.-qm)./8
        end
        selectdim(b,dim,2j-1).=q.-delta
        selectdim(b,dim,2j).=q.+delta
    end
    b
end
function prolong_faces(a,dim;periodic=false)
    n=size(a,dim)-(periodic ? 0 : 1);sz=collect(size(a));sz[dim]=2n+(periodic ? 0 : 1)
    b=Array{eltype(a)}(undef,sz...)
    function slab(j)
        periodic && return selectdim(a,dim,mod1(j,n))
        j==0 && return -selectdim(a,dim,2)
        j==n+2 && return -selectdim(a,dim,n)
        selectdim(a,dim,j)
    end
    for j=1:n
        selectdim(b,dim,2j-1).=slab(j)
        selectdim(b,dim,2j).=(slab(j).+slab(j+1))./2 .- (slab(j+2).-slab(j+1).-slab(j).+slab(j-1))./16
    end
    periodic || (selectdim(b,dim,2n+1).=slab(n+1))
    b
end
function average_pairs(a,dim)
    n=size(a,dim);iseven(n) || error("Pair average requires even size")
    (selectdim(a,dim,1:2:n).+selectdim(a,dim,2:2:n))./2
end
