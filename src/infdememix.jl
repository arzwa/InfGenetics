# Think of an environmental context (natural selection, density regulation, ...)
abstract type PopulationContext end

# ~ Barton & Etheridge
struct DirectionalSelection{T} <: PopulationContext
    b::T  # selection gradient
    # add density dependence?
end

# ~ Polechova & Barton
struct StabilizingSelection{T} <: PopulationContext
    θ::T    # phenotypic optimum
    ω::T    # variance of stabilizing selection
    r::T    # growth rate
    K::T    # carrying capacity (can be Inf when no density dependence)
end

(c::StabilizingSelection)(;θ=c.θ, ω=c.ω, r=c.r, K=c.K) = 
    StabilizingSelection(θ, ω, r, K)

@with_kw struct InfDemeMix{T}
    z::Vector{T}   # trait values (scalars, for now)
    F::Vector{T}   = zeros(length(z)) # inbreeding
    Φ::Matrix{T}   = zeros(length(z), length(z)) # coancestry
    m::Vector{Int} = fill(2, length(z)) # ploidy levels
    U::Matrix{T}   = [0.95 0.05 ; 0.0 0.0 ; 0.0 0.95]  # cytotype × gamete 
    w::Vector{T}   = [1., 1., 1.]  # cytotype specific viability 
    β::Vector{T}   = [1., 0.75, 0.5]     # variance scaler, will be in [0,1]
    γ::Vector{T}   = [1., √(2*β[2]/3), √(β[3]/2)]  # trait scalers
    V::T = 0.5     # half the diploid segvar
    α::T = 0.0     # double reduction rate
    ξ::T = 0.5     # probability of an IBD unreduced gamete
end
# U is a matrix with for each cytotype (rows) the probability of generating a
# particular (euploid) gamete.
# Note that `U` can be seen as capturing a fertility component of fitness,
# whereas `w` captures a viability component of fitness
# Cytotype-specific assortativity is one additional thing we may be interested
# in.

(d::InfDemeMix)(;z=d.z, F=d.F, Φ=d.Φ, m=d.m, U=d.U) = 
    InfDemeMix(z=z, F=F, Φ=Φ, m=m, U=U, w=d.w, V=d.V, β=d.β, α=d.α, ξ=d.ξ)

Base.getindex(d::InfDemeMix, i) = (d.z[i], d.m[i], d.F[i])
Base.length(d::InfDemeMix) = length(d.z)

function getβ(d, m)
    m == 2 && return 1. 
    m == 3 ? 4d.β/3 : d.β  # β4 = 0.5 -> β3 = 2/3
end

# compute the offspring distribution for the family of <i,j>
# this is, for each cytotype, compute the frequency (before selection on
# trait values), offspring mean (conditional on cytotype) and offspring
# variance (conditional on cytotype). 
struct Family{T}
    i::Int
    j::Int
    X::Matrix{T}
end

# first column of X are the frequencies for each gametic ploidy combination
# second column collects mean trait values (if freq != 0)
# third column collects segregation variances
freq(x::Family, i, j) = x.X[_toindex(i,j),1]
avg( x::Family, i, j) = x.X[_toindex(i,j),2]
var( x::Family, i, j) = x.X[_toindex(i,j),3]

# these two are isomorphisms
function _toindex(i, j)
    (i == 1 && j == 1) && return 1
    (i == 2 && j == 1) && return 2
    (i == 1 && j == 2) && return 3
    (i == 2 && j == 2) && return 4
    return 0  # will cause error further
end

function _toploidy(k)
    k == 1 && return 1, 1
    k == 2 && return 2, 1
    k == 3 && return 1, 2
    k == 4 && return 2, 2
    return 0, 0
end

# constructs a `Family`
function offspring_distribution(d::InfDemeMix{T}, i, j) where T
    @unpack β, γ = d
    zi, mi, Fi = d[i]
    zj, mj, Fj = d[j]
    X = fill(NaN, 4, 3)
    for ii=1:2, jj=1:2
        k = _toindex(ii, jj)
        m = ii + jj  # zygote ploidy
        X[k,1] = d.U[mi-1,ii] * d.U[mj-1,jj] * d.w[m-1]
        X[k,1] == 0. && continue
        #β = getβ(d, m)
        #X[k,2] = √β * ((ii/mi) * zi + (jj/mj) * zj)
        X[k,2] =  gameticvalue(γ, m, ii, mi, zi) + gameticvalue(γ, m, jj, mj, zj)
        X[k,3] =  β[m-1] * (segvar(d, ii, mi, Fi) + segvar(d, jj, mj, Fj))
    end
    return Family(i, j, X)
end

# mo = ploidy offspring, mg = ploidy gamete, mp = ploidy parent
gameticvalue(γ, mo, mg, mp, zp) = (γ[mo-1]/γ[mp-1]) * (mg/mp) * zp

# constructs all families
function offspring_distribution(d)
    return [offspring_distribution(d, i, j) for i=1:length(d), j=1:length(d)]
end

# a random offspring from a family
function randoff(f::Family{T}, d::InfDemeMix) where T
    mi, mj = _toploidy(sample(1:4, Weights(vec(f.X[:,1]))))
    z = rand(Normal(avg(f, mi, mj), √var(f, mi, mj)))
    F = Fcoeff(d, mi, mj, f.i, f.j)
    z, (mi + mj)::Int, F::T
end

# apply directional selection: changes the frequencies of each component
# (cytotype) by multiplying with environmental fitness, and shifts the trait
# distribution according to Gaussian-exponential convolution
# note that the fitness function is w(z) = exp(b*z) (z > 0 => w > 1)
#function directional_selection!(f::Family, b)
#    for i=1:2, j=1:2
#        f.X[i,j,1] == 0. && continue  # prevent NaN issues
#        f.X[i,j,1] *= (exp(b * avg(f, i, j) + (b^2/2) * var(f, i, j)))
#        f.X[i,j,2] += b * var(f, i, j)
#    end
#    return f
#end

𝔼fitness(family) = sum(family.X[:,1])

# Selection happens on a collection of families. We return a set of (possible
# non-unique) families, and the offspring generation can be obtained by
# sampling one offpsring individual from each.
function selection(d::InfDemeMix, context::StabilizingSelection)
    @unpack r, K = context
    fs  = offspring_distribution(d)
    fs′ = map(f->selection(f, context), fs)
    w̄ij = 𝔼fitness.(fs′)
    N   = length(d)
    w̄g  = sum(w̄ij)/N  # genetic component
    w̄e  = exp(r*(1-(N/K)))
    N′  = rand(Poisson(w̄e * w̄g))
    S   = sample(vec(fs′), Weights(vec(w̄ij)), N′)
end

function selection(f::Family, context::StabilizingSelection)
    @unpack θ, ω = context
    X = fill(NaN, 4, 3)
    for i=1:2, j=1:2  # each gametic ploidy combination
        k = _toindex(i, j)
        if f.X[k,1] == 0. 
            X[k,1] = 0.
            continue
        end
        z̄ = avg(f, i, j)
        v = var(f, i, j)
        ϕ = 1/((1/ω) + (1/v))
        w̄g = √(ω/(ω+v)) * exp(-(θ - z̄)^2/(2*(ω + v)))
        X[k,1] = w̄g * f.X[k,1]
        X[k,2] = (z̄/v + θ/ω)*ϕ
        X[k,3] = ϕ
    end
    return Family(f.i, f.j, X)
end

function segvar(d::InfDemeMix, m, mi, Fi)
    # mi = ploidy of individual i
    # Fi = inbreeding coeff of individual i
    # m  = ploidy of gamete
    if 2m == mi == 2  # diploid reduced
        return (1-Fi) * d.V
    elseif m == mi == 2  # diploid unreduced
        return 4 * d.ξ * (1-Fi) * d.V
    elseif 2m == mi == 4  # tetraploid reduced
        return (1-Fi) * (1+2d.α) * d.V
    elseif mi == 3 && m == 1
        return (1-Fi) * (8d.V/9)
    elseif mi == 3 && m == 2
        return (1-Fi) * (1/3 + d.ξ) * 8d.V/3
    else
        @error "meiosis not implemented $mi → $m" 
    end
end

function Fcoeff(d::InfDemeMix, mi, mj, i, j)
    # m ploidy of zygote
    # mi ploidy of the i-contributed gamete
    # mj ploidy of the j-contributed gamete
    _, mmi, Fi = d[i]
    _, mmj, Fj = d[j]
    m = mi + mj
    if m == 2
        return d.Φ[i,j]
    elseif m == 3
        return ((mi-1) * Fstar(mmi, Fi, d.ξ, d.α) + 
                (mj-1) * Fstar(mmj, Fj, d.ξ, d.α) +
                2d.Φ[i,j])/3
    elseif m == 4
        return (Fstar(mmi, Fi, d.ξ, d.α) + 
                Fstar(mmj, Fj, d.ξ, d.α) + 
                4d.Φ[i,j])/6
    end
end

# probability of sampling a diploid gamete with two IBD alleles
function Fstar(m, F, ξ, α)
    m == 2 && return F*(1-ξ) + ξ
    m == 3 && return F*(1-ξ) + ξ
    m == 4 && return F*(1-α) + α
    return NaN
end

# a single generation of random mating with selection
function generation(d::InfDemeMix, context::PopulationContext)
    N = length(d)
    N == 0 && return deepcopy(d)
    S = selection(d, context)
    N′= length(S)
    O = randoff.(S, Ref(d))
    z = first.(O)
    m = map(x->x[2], O)
    F = last.(O)
    P = zeros(N′,N) 
    for k=1:N′
        P[k,S[k].i] += 0.5
        P[k,S[k].j] += 0.5
    end
    Φ = P * d.Φ * P'
    setdiagΦ!(Φ, m, F)
    return d(z=z, F=F, Φ=Φ, m=m)
end

setdiagΦ!(Φ, m, F) = Φ[diagind(Φ)] .= (1 .+ (m .- 1) .* F) ./ m

