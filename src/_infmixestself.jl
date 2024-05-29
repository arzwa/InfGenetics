# Establishment model (migration + directional selection)
@with_kw struct InfDemeMixEstSelf{T}
    U::Matrix{T} = [0.95 0.05 ; 0.05 0.05 ; 0.0 0.95]  # cytotype × gamete 
    w::Vector{T} = ones(3)              # cytotype specific viability 
    β::Vector{T} = [1, √(2/3), √(1/2)]  # variance scaler, will be in [0,1]
    α::Vector{T} = zeros(3)             # double reduction rate
    σ::Vector{T} = zeros(3)             # selfing rates
    γ::T = 0.25                         # selection gradient
    m::T = 1.0                          # 𝔼 number of migrants/gen                          
    V::T = 0.5                          # Vx (haplotype variance in diploids)
    θ::Vector{T} = fill(-3.0, 3)        # trait mean of migrants
    c::Vector{T} = cytotype_equilibrium(U)
end

function fitnesses(z, c, M::InfDemeMixEstSelf)
    @unpack γ, U, w = M
    return [(U[k-1,1] + U[k-1,2])*w[k-1]*exp(γ*zk) for (zk,k) in zip(z,c)]  
end

# As in Abu Awad et al. 2014
function selfingp(w, c, σ)
    N = length(w)
    x = [(1-σ[ci-1])*wi for (wi, ci) in zip(w, c)]
    Z = sum(x)
    [σ[ci-1]*wi/(σ[ci-1]*wi + (Z - xi)/(N-1)) for (wi, xi, ci) in zip(w, x, c)]
end

function samplediff(ws, i)
end

function generation(pop::InfPop{T}, M::InfDemeMixEstSelf) where T
    @unpack θ, α, β, m, V = M
    # Poisson number of migrants
    nr = length(pop.z)
    nm = rand(Poisson(m))
    cm = rand(Multinomial(nm, M.c))
    cm = vcat([fill(i+1, k) for (i, k) in enumerate(cm)]...)
    zm = [rand(Normal(θ[m-1], β[m-1]*√V)) for m in cm]
    z  = [pop.z ; zm]
    c  = [pop.c ; cm]
    F  = [pop.F ; zeros(nm)]
    Φ  = [pop.Φ zeros(nr, nm); zeros(nm, nr) zeros(nm, nm)]
    # fitness, this is the contribution to the gamete pool
    ws = fitnesses(z, c, M)
    sp = selfingp(ws, c, M.σ)
    # the number of female gametes produced is proportional to fitness
    # we assume an infinite pool of male gametes, to which each individual
    # contributes proportional to fitness
    w̄  = length(ws) == 0 ? 0.0 : mean(ws)
    # number of individuals in the next generation = number of fertilized ovules
    N  = rand(Poisson(w̄*length(z)))   
    # sample N mothers and N potential fathers proportional to fitness
    # XXX 'wasteful' random numbers when there's a lot of samples, see below
    ix = sample(1:length(ws), Weights(ws), N, replace=true) 
    P  = zeros(N,nr+nm)
    pop_ = InfPop(z=zeros(N), F=zeros(N), Φ=zeros(N,N), c=zeros(Int, N))
    for k=1:N
        i = ix[k]  # mother individual
        j = rand() < sp[i] ? i : sample(1:length(ws), Weights(ws))  # samplediff(ws, i)
        # XXX would be better to have strict outcrossing vs. selfing...
        # But than `N` may not be right --> assign fitnesses to pairs of
        # individuals after all...
        zi, Fi, ci = z[i], F[i], c[i]
        zj, Fj, cj = z[j], F[j], c[j]
        gi, yi, Vi = gamete(M, zi, Fi, ci)
        gj, yj, Vj = gamete(M, zj, Fj, cj)
        cij = gi + gj
        βij = β[cij-1]
        βi  = β[ci-1]
        βj  = β[cj-1]
        Yi  = rand(Normal(yi*βij/βi, βij*√Vi)) 
        Yj  = rand(Normal(yj*βij/βj, βij*√Vj)) 
        pop_.z[k] = Yi + Yj
        pop_.c[k] = cij
        P[k,i] += gi/cij
        P[k,j] += gj/cij
        if cij == 2
            pop_.F[k] = Φ[i,j]
        elseif cij == 3
            (h, ch) = gi == 2 ? (i, ci) : (j, cj)
            Fh = (1-α[ch-1])*F[h] + α[ch-1]
            pop_.F[k] = (Fh + 2Φ[i,j])/3
        elseif cij == 4
            Fi = (1-α[ci-1])*F[i] + α[ci-1]
            Fj = (1-α[cj-1])*F[j] + α[cj-1]
            pop_.F[k] = (Fi + Fj + 4Φ[i,j])/6
        end
    end
    Φ_ = P * Φ * P'
    for i=1:N
        Φ_[i,i] = (1/pop_.c[i]) * (1 + (pop_.c[i] - 1) * pop_.F[i]) 
    end
    pop_.Φ .= Φ_
    return pop_
end

