# Establishment model (migration + directional selection)
# With fitness determined for pairs of individuals, this allows to better
# incorporate selfing and assortative mating... (different sexes requires a bit
# more work...)
# XXX Should study more carefully how assortative mating/selfing are usually
# modelled.
@with_kw struct InfDemeMixEstPair{T}
    U::Matrix{T} = [0.95 0.05 ; 0.05 0.05 ; 0.0 0.95]  # cytotype × gamete 
    w::Vector{T} = ones(3)              # cytotype specific viability 
    β::Vector{T} = [1, √(2/3), √(1/2)]  # variance scaler, will be in [0,1]
    α::Vector{T} = zeros(3)             # double reduction rate
    σ::Vector{T} = zeros(3)             # selfing rate by ploidy
    ρ::Matrix{T} = ones(3,3)            # assortativity by ploidy
    γ::T = 0.25                         # selection gradient
    m::T = 1.0                          # 𝔼 number of migrants/gen                          
    V::T = 0.5                          # Vx (haplotype variance in diploids)
    θ::Vector{T} = fill(-3.0, 3)        # trait mean of migrants
    c::Vector{T} = cytotype_equilibrium(U)
end

# But how exactly do we do this?
# 1. We let the number of offspring in the next generation be Poisson(NW̄) where
#    W̄ is the viability component. But note this is not what we've done
#    earlier, we also included the fertility component in the offspring
#    number (i.e. `uk1 + uk2`)
# 2. We take W̄ = E[wᵢⱼ] = ∑∑wᵢⱼ/N², hence NW̄ = ∑∑wᵢⱼ/N, where wᵢⱼ includes both
#    fertility, viability of offspring/parents, ...
# Note, B&E use the expected viability fitness of the offspring, but one could
# also take the viability of the parents? That is somewhat easier, as we don't
# need to consider the ploidy levels of the offspring etc.

# XXX: ∑ⱼ Wᵢⱼ = expected number of offspring from individual i as mother
# so sum(W, dims=2) gives the expected number of fertilized ovules for each
# individual.
# Here we assume that the selfing rate parameter is the realized selfing rate.
# This may not make so much sense when we are thinking of a SI system.
function fitnesses(z, c, M)
    @unpack σ, ρ, γ, w, U = M
    N = length(z)
    W = zeros(N,N)
    # fitness values
    Y = [w[c[j]-1]*exp(γ*z[j])*(U[c[j]-1,1] + U[c[j]-1,2]) for j=1:N]
    Z = sum(Y)
    for i=1:N
        Zi = (Z - Y[i])  # mean fitness of others * (N-1)
        for j=1:N
            ci = c[i]-1; cj=c[j]-1
            W[i,j] = if i == j 
                σ[ci]*Y[i]
            else
                (1-σ[ci])*ρ[ci,cj]*Y[i]*Y[j]/Zi
                # P(father = j|outcrossing, mother i)P(outcrossing| mother i) P(mother i)
            end
        end
    end
    return W
end

# XXX Here we assume a SI system: if self-compatible, we have a proportion of
# self-fertilizing ovules, and a proportion of randomly fertilized ovules, but
# random fertilization may also result in selfing (i.e. we draw pollen from an
# infinite gamete pool).
function fitnessesSI(z, c, M)
    @unpack σ, ρ, γ, w, U = M
    N = length(z)
    W = zeros(N,N)
    # fitness values
    Y = [w[c[j]-1]*exp(γ*z[j])*(U[c[j]-1,1] + U[c[j]-1,2]) for j=1:N]
    Z = sum(Y)
    for i=1:N
        for j=1:N
            ci=c[i]-1; cj=c[j]-1
            W[i,j] = if i == j && σ[ci] == 0.0  # self-incompatible
                0.0
            elseif i == j && σ[ci] > 0.0  # self-compatible
                σ[ci]*Y[i] + (1-σ[ci])*ρ[ci,cj]*Y[i]^2/Z
            else
                (1-σ[ci])*ρ[ci,cj]*Y[i]*Y[j]/Z
            end
        end
    end
    return W
end

#function pairfitness(zi, zj, ci, cj, σ, γ, w, ρ, U)
#    wi  = w[ci-1]*exp(γ*zi)*(U[ci-1,1] + U[ci-1,2])
#    wj  = w[cj-1]*exp(γ*zj)*(U[cj-1,1] + U[cj-1,2])
#    wij = wi*wj*ρ[ci-1,cj-1]
#    return wij
#end

function generation(pop::InfPop{T}, M::InfDemeMixEstPair) where T
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
    N  = length(z)
    # fitness for parental *pairs*
    W  = fitnessesSI(z, c, M)
    NW̄  = N == 0 ? 0.0 : sum(W)
    # number of individuals in the next generation = sum(W)/N = mean(W)*N
    N_ = rand(Poisson(NW̄))   
    # N_ parental pairs
    ws = vec(W)
    ix = sample(1:length(ws), Weights(ws), N_, replace=true) 
    P  = zeros(N_, N)
    pop_ = InfPop(z=zeros(N_), F=zeros(N_), Φ=zeros(N_,N_), c=zeros(Int, N_))
    cart = CartesianIndices((1:N, 1:N))
    for k=1:N_
        # to cartesian index...
        i, j = cart[ix[k]].I
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
    for i=1:N_
        Φ_[i,i] = (1/pop_.c[i]) * (1 + (pop_.c[i] - 1) * pop_.F[i]) 
    end
    pop_.Φ .= Φ_
    return pop_
end

