
export InfDeme, randgen, evolve

function evolve(fun, x, n)
    for i=1:n
        x = fun(x)
    end
    return x
end

"""
    InfDeme

A finite random mating population consisting of a single cytotype with
individuals characterized by a one-dimensional trait obeying infinitesimal
quantitative genetics. Defined for haploid, diploid and tetraploid cytotypes.
"""
struct InfDeme{K,T}
    z::Vector{T}
    F::Matrix{T}
    V::T
end

InfDeme(z::Vector{T}, V::T; k=1, F=zeros(length(z),length(z))) where T = 
    InfDeme{Val{k},T}(z, F, V)

(d::InfDeme{K,T})(z, F; V=d.V) where {K,T} = InfDeme{K,T}(z, F, V)

Base.length(d::InfDeme) = length(d.z)
ploidy(d::InfDeme{K}) where K = K.parameters[1] 

function randgen_noinbreeding(d)
    z = similar(d.z)
    N = length(d)
    for k=1:N
        @inbounds z[k] = randoff(d, rand(1:N), rand(1:N))
    end
    return d(z, d.F)
end

function randgen(d::InfDeme)
    z = similar(d.z)
    P = zeros(size(d.F)) 
    N = length(d)
    Fd = zeros(length(z))
    for k=1:N
        i = rand(1:N)
        j = rand(1:N)
        P[k,i] += 0.5
        P[k,j] += 0.5
        z[k] = randoff(d, i, j)
        Fd[k] = Fdiag(d, i, j)
    end
    F = update_F(d, P)
    F[diagind(F)] .= Fd
    return d(z, F)
end

function randoff(d::InfDeme{Val{1}}, i, j)
    Cij = d.V * (1. - d.F[i,j]) 
    zmp = (d.z[i] + d.z[j]) / 2
    return rand(Normal(zmp, √Cij))
end

function randoff(d, i, j)
    Cij = d.V * (1. - (d.F[i,i] + d.F[j,j])/2)
    zmp = (d.z[i] + d.z[j]) / 2
    return rand(Normal(zmp, √Cij))
end

Fdiag(d::InfDeme{Val{1}}, i, j) = 1.
Fdiag(d::InfDeme{Val{2}}, i, j) = i == j ? 
    0.50 * (1 + d.F[i,i]) : d.F[i,j]
Fdiag(d::InfDeme{Val{4}}, i, j) = i == j ? 
    0.25 * (1 + 3d.F[i,i]) : (d.F[i,i] + d.F[j,j] + 4d.F[i,j])/6

update_F(d::InfDeme{Val{1}}, P) = P * d.F * P'
update_F(d::InfDeme{Val{K}}, P) where K = P * (d.F + (I - Diagonal(d.F))/K) * P' 


# Mixed-ploidy deme
# =================
struct InfDemeMixed{T}
    z::Vector{T}
    F::Matrix{T}
    U::Matrix{T}   # unreduced gamete rate matrix
    w::Vector{T}   # ploidy specific viability (Pr germination)
    k::Vector{Int} # ploidy levels 
    V::Vector{T}   # base 1/2 segregation variance for different ploidy levels
    β::Vector{T}   # allelic effect scalers
    ξ::T           # 1 - c - f + (3/2)cf
end

(d::InfDemeMixed)(z, F, k) = InfDemeMixed(z, F, d.U, d.w, k, d.V, d.β, d.ξ)

Base.length(d::InfDemeMixed) = length(d.z)

function InfDemeMixed(z, U, V0, β, ξ; w=[0., 1., 0., 1.], k=fill(2, length(z)))
    N = length(z)
    F = zeros(N, N)
    V = [(m/2)*β[m]*V0 for m=1:length(w)]
    InfDemeMixed(z, F, U, w, k, V, β, ξ) 
end

# ploidy level of a random gamete from i
randk(d::InfDemeMixed, i) = sample(1:4, Weights(d.U[d.k[i],:]))

# contribution to the segvar from parent i for offspring of ploidy ko
function segvar(d::InfDemeMixed, i, ko) 
    if 2d.k[i] == ko 
        d.V[ko] * d.ξ * (1 - d.F[i,i])
    else
        d.V[ko] * (1 - d.F[i,i])
    end
end

function offmean(d::InfDemeMixed, i, ko)
    if d.k[i] == ko
        d.z[i] / 2
    elseif 2d.k[i] == ko 
        d.β[ko] * d.z[i]
    end
end

function getFi(d::InfDemeMixed, i, ko)
    if 2d.k[i] == ko 
        d.F[i,i] * (1 - d.ξ) + d.ξ
    else
        d.F[i,i]
    end
end

function Fdiag(d::InfDemeMixed, i, j, ko)
    Fi = getFi(d::InfDemeMixed, i, ko)
    Fj = getFi(d::InfDemeMixed, j, ko)
    return if i == j  # selfing, but correct if selfing of unreduced gametes?
        (1/ko)*(1 + (ko - 1)*Fi)
    elseif ko == 2
        d.F[i,j]
    elseif ko == 4
        (Fi + Fj + d.F[i,j])/6
    end
end

function randoff(d::InfDemeMixed, i, j)
    ki = randk(d, i)
    kj = randk(d, j)
    ko = ki + kj
    rand() > d.w[ko] && return -1, -1, -1, -1.0
    β  = d.β[ko]
    Cij = segvar(d, i, ko) + segvar(d, j, ko)
    Zij = offmean(d, i, ko) + offmean(d, j, ko)
    ko, ki, kj, rand(Normal(Zij, √Cij))
end

function update_F(d::InfDemeMixed, P)
    P * (d.F + (I - Diagonal(d.F)) ./ d.k) * P'
end

function randgen(d)
    z = similar(d.z)
    k = similar(d.k)
    P = zeros(size(d.F)) 
    N = length(d)
    Fd = zeros(length(z))
    o = 1
    while o <= N
        i = rand(1:N)
        j = rand(1:N)
        ko, kki, kkj, zo = randoff(d, i, j)
        ko == -1 && continue
        P[o,i] += kki/ko
        P[o,j] += kkj/ko
        k[o] = ko
        z[o] = zo
        Fd[o] = Fdiag(d, i, j, ko)
        o += 1
    end
    F = update_F(d, P)
    F[diagind(F)] .= Fd
    return d(z, F, k)
end


# Tetraploids (obsolete and wrong)
# ================================
struct InfDemeTetra{T}
    z::Vector{T}
    G::Vector{T}
    F::Matrix{T}  # symmetric, use special type?
    V::T
end

function InfDemeTetra(z, V) 
    N = length(z) 
    InfDemeTetra(z, ones(N), zeros(N,N), V)
end

Base.length(d::InfDemeTetra) = length(d.z)

segvar(d::InfDemeTetra) = mean(d.G) * d.V

function randgen(d::InfDemeTetra) 
    N = length(d)
    z = similar(d.z)
    G = similar(d.G)
    P = zeros(N, N)
    F = similar(d.F)
    D = similar(F)
    Fd = Vector{eltype(F)}(undef, N) 
    for k=1:N, l=k:N
        @inbounds D[k,l] = D[l,k] = (d.F[k,k] + d.F[l,l] + 4d.F[k,l])/6
    end
    for k=1:N
        i = rand(1:N)
        j = rand(1:N)
        Cij = d.V * (d.G[i] + d.G[j]) / 2 
        zmp = (d.z[i] + d.z[j]) / 2
        P[k,i] += 0.5
        P[k,j] += 0.5
        z[k] = rand(Normal(zmp, √Cij))
        G[k] = 1 - D[i,j]
    end
    F = P * (d.F + (0.25I - Diagonal(d.F)/4)) * P'
    #F[diagind(F)] .= diag(P * D * P')
    #return InfDemeTetra(z, G, F, d.V), (z=z, F=F, D=D, P=P, G=G)
    return InfDemeTetra(z, G, F, d.V)
end
# I have the feeling I am doing some redundant stuff in here, but not sure
# exactly where


evolve_and_track_segvar(d::InfDemeTetra, t) = map(1:t) do gen
    d = generate_offspring(d)
    d.V * mean(d.G)
end



