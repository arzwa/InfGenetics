module InfGenetics

using LinearAlgebra, Distributions

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

function randgen(d)
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



end # module
