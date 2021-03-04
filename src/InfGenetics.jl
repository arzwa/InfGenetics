module InfGenetics

using LinearAlgebra, Distributions

# Haploids and diploids
# =====================
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

# haploids, neutral generation, constant popsize
function generate_offspring(d)
    z = similar(d.z)
    P = zeros(size(d.F)) 
    N = length(d)
    for k=1:N
        i = rand(1:N)
        j = rand(1:N)
        P[k,i] += 0.5
        P[k,j] += 0.5
        z[k] = randoff(d, i, j)
    end
    F = update_F(d, P)
    return d(z, F) 
end

function randoff(d::InfDeme{Val{1}}, i, j)
    Cij = d.V * (1. - d.F[i,j]) 
    zmp = (d.z[i] + d.z[j]) / 2
    return rand(Normal(zmp, Cij))
end

function update_F(d::InfDeme{Val{1}}, P)
    F = P * d.F * P'
    F[diagind(F)] .= 1.
    return F
end

function randoff(d::InfDeme{Val{2}}, i, j)
    Cij = d.V * (1. - (d.F[i,i] + d.F[j,j])/2)
    zmp = (d.z[i] + d.z[j]) / 2
    return rand(Normal(zmp, Cij))
end

function update_F(d::InfDeme{Val{2}}, P)
    P * (d.F + (0.5I - Diagonal(d.F)/2)) * P' 
end

# Tetraploids
# ===========

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

function generate_offspring(d::InfDemeTetra) 
    z = similar(d.z)
    G = similar(d.G)
    F = similar(d.F)
    P = zeros(size(d.F))
    N = length(d)
    D = similar(F)
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
        z[k] = rand(Normal(zmp, Cij))
        G[k] = 1 - D[i,j]
    end
    F = P * (d.F + (0.25I - Diagonal(d.F)/4)) * P'
    F[diagind(F)] .= diag(P * D * P')
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
