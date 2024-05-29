
"""
    InfDeme

A finite random mating population consisting of a single cytotype with
individuals characterized by a one-dimensional trait obeying infinitesimal
quantitative genetics. Defined for haploid, diploid and tetraploid cytotypes.
"""
struct InfDeme{K,T}  # K is ploidy level
    z::Vector{T}  # trait value (scalar, at least for now)
    F::Vector{T}  # inbreeding coefficients
    Φ::Matrix{T}  # coancestry coefficients
    V::T          # segregation variance
    α::T          # probability of double reduction
end

# constructor
InfDeme(z::Vector{T};
        k = 1,
        V = 1.,
        F = zeros(length(z)),
        Φ = zeros(length(z),length(z)),
        α = 0. 
       ) where T = 
    InfDeme{Val{k},T}(z, F, Φ, V, α)

(d::InfDeme{K,T})(z, F, Φ; V=d.V, α=d.α) where {K,T} = 
    InfDeme{K,T}(z, F, Φ, V, α)

popsize(d::InfDeme) = length(d.z)
ploidy(d::InfDeme{K}) where K = K.parameters[1] 

function wfgen(d::InfDeme)
    z = similar(d.z)
    N = popsize(d)
    P = zeros(N,N) 
    F = zeros(length(z))
    for k=1:N
        i = rand(1:N)
        j = rand(1:N)
        P[k,i] += 0.5
        P[k,j] += 0.5
        z[k] = randoff(d, i, j)
        F[k] = Fcoeff(d, i, j) 
    end
    Φ = P * d.Φ * P'
    setdiagΦ!(d, Φ, F)
    return d(z, F, Φ)
end

function randoff(d, i, j)
    Cij = segvar(d, i, j)
    zmp = (d.z[i] + d.z[j]) / 2
    return rand(Normal(zmp, √Cij))
end

Fcoeff(d::InfDeme{Val{1}}, i, j) = 1.
Fcoeff(d::InfDeme{Val{2}}, i, j) = d.Φ[i,j]
Fcoeff(d::InfDeme{Val{4}}, i, j) = (2d.α + (1-d.α) * (d.F[i] + d.F[j]) + 4d.Φ[i,j])/6

segvar(d::InfDeme{Val{1}}, i, j) = d.V * (1. - d.Φ[i,j])
segvar(d::InfDeme{Val{2}}, i, j) = d.V * (1. - (d.F[i] + d.F[j])/2)
segvar(d::InfDeme{Val{4}}, i, j) = d.V * (1. - (d.F[i] + d.F[j])/2) * (1 + 2d.α)

setdiagΦ!(d::InfDeme{Val{1}}, Φ, F) = Φ[diagind(Φ)] .= 1.
setdiagΦ!(d::InfDeme{Val{2}}, Φ, F) = Φ[diagind(Φ)] .= 0.5  .+ 0.5  .* F
setdiagΦ!(d::InfDeme{Val{4}}, Φ, F) = Φ[diagind(Φ)] .= 0.25 .+ 0.75 .* F

popstats(d::InfDeme) = (mz=mean(d.z), vz=var(d.z), mF=mean(d.F))

