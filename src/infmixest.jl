# Establishment model (migration + directional selection)
# Note that while it is straightforward to include density regulation, this
# does not make much sense with directional selection. So density regulation
# should go hand in hand with stabilizing selection. 
@with_kw struct InfDemeMixEst{T}
    U::Matrix{T} = [0.95 0.05 ; 0.05 0.05 ; 0.0 0.95]  # cytotype × gamete 
    w::Vector{T} = ones(3)              # cytotype specific viability 
    β::Vector{T} = [1, √(2/3), √(1/2)]  # variance scaler, will be in [0,1]
    α::Vector{T} = zeros(3)             # double reduction rate
    σ::Vector{T} = zeros(3)             # selfing rate -> zero means random selfing, NaN means SI
    ρ::Vector{T} = zeros(3)             # assortative mating by ploidy
    γ::T = 0.25                         # selection gradient
    m::T = 1.0                          # 𝔼 number of migrants/gen                          
    V::T = 0.5                          # Vx (haplotype variance in diploids)
    θ::Vector{T} = ones(3)              # trait value at which growth rate becomes positive
    c::Vector{T} = cytotype_equilibrium(U)
end

function parameters(M::InfDemeMixEst)
    @unpack U, β, α, σ, ρ, γ, m, V, θ, c, w = M 
    w2, w3, w4 = w
    b2, b3, b4 = β
    s2, s3, s4 = σ
    a2, a3, a4 = α
    r2, r3, r4 = ρ
    z2, z3, z4 = θ
    c2, c3, c4 = c
    (u21=U[1,1], u22=U[1,2], u31=U[2,1], u32=U[2,2], u41=U[3,1], u42=U[3,2],
     w2=w2, w3=w3, w4=w4, b2=b2, b3=b3, b4=b4, s2=s2, s3=s3, s4=s4, r2=r2,
     r3=r3, r4=r4, z2=z2, z3=z3, z4=z4, c2=c2, c3=c3, c4=c4, V=V, g=γ, 
     a2=a2, a3=a3, a4=a4, m=m)
end

function cytotype_equilibrium(U)
    u21, u22 = U[1,:]
    u31, u32 = U[2,:]
      _, u42 = U[3,:]
    g1 = (u21 - 4*u31 - 2*u32 + 2*u42)/(
        2*(u21 + u22 - 2*u31 - 2*u32 + u42)) - sqrt(
            u21^2 - 4*u21*u32 + 8*u22*u31 - 4*u22*u42 + 4*u32^2)/(
                2*(-u21 - u22 + 2*u31 + 2*u32 - u42))
    [g1^2, 2g1*(1-g1), (1-g1)^2] 
end

Ew(z̄, V, γ, θ) = exp(γ*(z̄ - θ) + (γ^2/2)*V)
hapsegvar(c, F, V)    = 1 < c < 4 ? (c-1)/c * (1-F) * V           : 0.0
dipsegvar(c, α, F, V) = 1 < c ≤ 4 ? (2/c)*(c*(1 + α) - 2)*(1-F)*V : 0.0
# NOTE: dipsegvar:
# tet (2/4)(2 + 4α)(1-F)V 
# tri (2/3)(1 + 3α)(1-F)V
# dip (2/2)(0 + 2α)(1-F)V
# hence
#    (2/c)(c-2 + cα)(1-F)V
#   =(2/c)(c(1+α)-2)(1-F)V

"""
    Zᵢⱼ ∼ N{βᵢⱼ[(gᵢ/cᵢ)(zᵢ/βᵢ) + (gⱼ/cⱼ)(zⱼ/βⱼ)], Vᵢ + Vⱼ}
"""
function family_distribution(M::InfDemeMixEst, xi, xj)
    @unpack U, β, γ, V, α, w, θ = M 
    zi, Fi, ci = xi
    zj, Fj, cj = xj
    Vi1 = hapsegvar(ci, Fi, V)
    Vi2 = dipsegvar(ci, α[ci-1], Fi, V)
    Vj1 = hapsegvar(cj, Fj, V)
    Vj2 = dipsegvar(cj, α[cj-1], Fj, V)
    # diploid offspring
    z2 = β[1] * ((1/ci)*(zi/β[ci-1]) + (1/cj)*(zj/β[cj-1]))
    V2 = β[1]^2 * (Vi1 + Vj1)
    w2  = Ew(z2, V2, γ, θ[1]) * U[ci-1,1] * U[cj-1,1] * w[1]
    # triploid offspring 1
    z31 = β[2] * ((2/ci)*(zi/β[ci-1]) + (1/cj)*(zj/β[cj-1]))
    V31 = β[2]^2 * (Vi2 + Vj1)
    w31 = Ew(z31, V31, γ, θ[2]) * U[ci-1,2] * U[cj-1,1] * w[2]
    # triploid offspring 1
    z32 = β[2] * ((1/ci)*(zi/β[ci-1]) + (2/cj)*(zj/β[cj-1]))
    V32 = β[2]^2 * (Vi1 + Vj2)
    w32 = Ew(z32, V32, γ, θ[2]) * U[ci-1,1] * U[cj-1,2] * w[2]
    # tetraploid offspring
    z4 = β[3] * ((2/ci)*(zi/β[ci-1]) + (2/cj)*(zj/β[cj-1]))
    V4 = β[3]^2 * (Vi2 + Vj2)
    w4  = Ew(z4, V4, γ, θ[3]) * U[ci-1,2] * U[cj-1,2] * w[3]
    [w2, w31, w32, w4], [z2, z31, z32, z4], [V2, V31, V32, V4]
end

function families(M::InfDemeMixEst, pop::InfPop)
    # Fitness/family distributions for all parental pairs
    N = popsize(pop) 
    W = Array{Float64}(undef, N, N, 4)
    Z = Array{Float64}(undef, N, N, 4)
    V = Array{Float64}(undef, N, N, 4)
    for i=1:N
        for j=i:N
            ws, zs, Vs = family_distribution(M, pop[i], pop[j])
            W[i,j,:] .= ws
            W[j,i,:] .= ws
            Z[i,j,:] .= zs
            Z[j,i,:] .= zs
            V[i,j,:] .= Vs
            V[j,i,:] .= Vs
        end
    end
    W ./ N, Z, V
end

"""
    selfing!(W, σ, c)

- σ[i] = NaN  => self incompatible without intrinsic disadvantage 
- σ[i] = -Inf => self-incompatible with intrinsic disadvantage
- σ[i] = 0    => random selfing (probability 1/N)
- σ[i] > 0    => partial selfing through female function, pollen is random
                 irrespective of selfing rates
"""
function selfing!(W, σ, c)
    # (i,j) entry has i as mother and j as father
    # selfing is modeled as proportion of self-fertilized ovules, we assume the
    # pollen pool is unaffected by selfing
    # so  (i,j) i!=j entry (1-σᵢ)*E[wij]/N
    # and (i,i) entry σᵢE[wij] + (1-σᵢ)E[wii]/N
    N = size(W,1)
    s = [σ[c[i]-1] for i=1:N]  # selfing rates
    for i=1:N
        if isnan(s[i])  # self-incompatible without intrinsic disadvantage
            # If we don't want an intrinsic disadvantage to SI (except failure
            # if alone), then we have to account for N=1 case...
            r = N == 1 ? 0.0 : N/(N-1)
            W[i,:,:] .*= r
            W[i,i,:] .= 0.0
        elseif !isfinite(s[i])  # self-incompatible with intrinsic disadvantage
            W[i,i,:] .= 0.0
        else  # partial selfing, outcrossing pollen is random
            Wiself = W[i,i,:] .* s[i] .* N
            W[i,:,:] .= W[i,:,:] .* (1 - s[i])
            W[i,i,:] .+= Wiself
        end
    end
    return W
end

function assortative_mating!(W, ρ, c)
    N = size(W,1)
    x = N ./ counts(c, 2:4)
    for i=1:N
        k = c[i]-1
        for j=1:N
            if c[i] == c[j]
                W[i,j,:] .= W[i,j,:] .* (ρ[k] * x[k]) .+ W[i,j,:] .* (1-ρ[k])
            else
                W[i,j,:] .*= (1 - ρ[k])
            end
        end
    end
    return W
end

function migration(M::InfDemeMixEst, pop::InfPop)
    @unpack m, θ, V, c, β, α = M
    @unpack z, c, F, Φ = pop
    nr = popsize(pop)
    nm = rand(Poisson(m))
    add_n_migrants(M, nm, pop)
end

function add_n_migrants(M::InfDemeMixEst, nm, pop::InfPop)
    @unpack m, θ, V, c, β, α = M
    @unpack z, c, F, Φ = pop
    cm = rand(Multinomial(nm, M.c))
    cm = vcat([fill(i+1, k) for (i, k) in enumerate(cm)]...)
    # Vzₖ/Vz₂ = (k/2)βₖ² => Vzₖ = kβₖ²V
    # But, this is not really the case: this is an equilibrium result for an
    # isolated monocytotypic population. It does hold for α = 0...
    # Here we use the small u approximation
    Vs = (β .^ 2) .* [2,3,4] .* V .* [1, 1+2α[1]/3, 1+α[1]]  
    zm = [rand(Normal(0.0, √Vs[m-1])) for m in cm]
    add_unrelated_individuals(pop, zm, cm, zeros(nm))
end

function add_unrelated_individuals(pop, zm, cm, Fm)
    nr = popsize(pop)
    nm = length(zm)
    z  = [pop.z ; zm]
    c  = [pop.c ; cm]
    F  = [pop.F ; Fm]
    Φm = diagm(fill(0.5, nm))  # coancestry with self is 0.5!
    Φ  = [pop.Φ zeros(nr, nm); zeros(nm, nr) Φm]
    InfPop(z=z, c=c, F=F, Φ=Φ)
end

#function _otherfitness(M, pop)
#    @unpack θ, α, β, γ, w, U, m, V = M
#    [(U[k-1,1] + U[k-1,2])*w[k-1]*exp(γ*zk) for (zk,k) in zip(pop.z,pop.c)]  
#end

function generation(M::InfDemeMixEst, pop::InfPop{T}; Nmax=100000) where T
    @unpack θ, α, β, γ, σ, ρ, U, m, V = M
    # pop_ is after migration
    # _pop is next generation
    pop_ = migration(M, pop)
    N_   = popsize(pop_)
    N_ == 0 &&  return InfPop()  # XXX popsize after migration = 0 ? Nothing can happen
    W, Z, V = fitnesses(M, pop_)
#   𝔼N = isempty(W) ? 0.0 : sum(W) 
    𝔼N = sum(W)  
    # number of individuals in the next generation
    N = min(rand(Poisson(𝔼N)), Nmax)
    N == 0 && return InfPop()    # XXX popsize next generation is zero
    # sample parental pairs x offspring ploidy combinations
    P = zeros(N, N_)
    C = CartesianIndices((1:N_, 1:N_, 1:4))
    idx = sample(C, Weights(vec(W)), N)
    _pop = InfPop(z=zeros(N), F=zeros(N), Φ=zeros(N,N), c=zeros(Int, N))
    @unpack F, Φ = pop_
    for (k,ij) in enumerate(idx)
        (i, j, x) = Tuple(ij)
        ci = pop_[i].c
        cj = pop_[j].c
        cij = x == 1 || x == 2 ? x + 1 : x   # offspring ploidy
        z̄ij = Z[ij]
        Vij = V[ij]
#        zij = rand(Normal(z̄ij, √Vij))
        zij = rand(Normal(z̄ij + γ*Vij, √Vij))  # XXX 2025-03-28, need trait value after selection...
        _pop.z[k] = zij
        _pop.c[k] = cij
        if x == 1
            P[k,i] += 1/2
            P[k,j] += 1/2
            _pop.F[k] = Φ[i,j]
        elseif x == 2
            P[k,i] += 2/3
            P[k,j] += 1/3
            Fi = (1 - α[ci-1])*F[i] + α[ci-1]
            _pop.F[k] = (Fi + 2Φ[i,j])/3
        elseif x == 3
            P[k,i] += 1/3
            P[k,j] += 2/3
            Fj = (1 - α[cj-1])*F[j] + α[cj-1]
            _pop.F[k] = (Fj + 2Φ[i,j])/3
        elseif x == 4
            P[k,i] += 1/2
            P[k,j] += 1/2
            Fi = (1 - α[ci-1])*F[i] + α[ci-1]
            Fj = (1 - α[cj-1])*F[j] + α[cj-1]
            _pop.F[k] = (Fi + Fj + 4Φ[i,j])/6
        end
    end
    # inplace? little benefit
    Φ_ = P * Φ * P'
    #mul!(_pop.Φ, P, Φ_)
    for i=1:N
        #_pop.Φ[i,i] = (1/_pop.c[i]) * (1 + (_pop.c[i] - 1) * _pop.F[i]) 
        Φ_[i,i] = (1/_pop.c[i]) * (1 + (_pop.c[i] - 1) * _pop.F[i]) 
    end
    _pop.Φ .= Φ_
    return _pop
end

function fitnesses(M, pop)
    @unpack σ, ρ = M 
    N = popsize(pop)
    W, Z, V = families(M, pop)
    W = selfing!(W, σ, pop.c)
    W = assortative_mating!(W, ρ, pop.c)
    return W, Z, V
end

