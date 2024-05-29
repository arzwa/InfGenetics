# Deme with hard selection, stabilizing selection, and density regulation
@with_kw struct InfDemeHS{T}
    U::Matrix{T} = [0.95 0.05 ; 0.05 0.05 ; 0.0 0.95]  # cytotype × gamete 
    w::Vector{T} = ones(3)              # cytotype specific Malthusian growth rate
    β::Vector{T} = [1, √(2/3), √(1/2)]  # variance scaler, will be in [0,1]
    α::Vector{T} = zeros(3)             # double reduction rate
    θ::Vector{T} = zeros(3)             # optimum
    γ::T = 0.25                         # selection gradient
    V::T = 0.5                          # Vx (haplotype variance in diploids)
    K::T = 200.                         # carrying capacity
    c::Vector{T} = cytotype_equilibrium(U)
end

function fitnesses(z, c, M::InfDemeHS)
    @unpack w, U, θ, γ, K = M
    N = length(z)
    return [(U[k-1,1] + U[k-1,2])*exp(w[k-1]*(1 - N/K) - γ*(zk - θ[k-1])^2) 
        for (zk,k) in zip(z,c)]  
end

function generation(pop::InfPop{T}, M::InfDemeHS) where T
    @unpack θ, α, β, γ, w, U, K, V = M
    @unpack z, F, Φ, c = pop
    N = length(z)
    # fitness, this is the contribution to the gamete pool
    ws = fitnesses(z, c, M)
    w̄  = length(ws) == 0 ? 0.0 : mean(ws)
    # number of individuals in the next generation
    N_ = rand(Poisson(w̄*length(z)))   
    # 2N gametes
    ix = sample(1:length(ws), Weights(ws), 2N_, replace=true) 
    P  = zeros(N_,N)
    pop_ = InfPop(z=zeros(N_), F=zeros(N_), Φ=zeros(N_,N_), c=zeros(Int, N_))
    for k=1:N_
        i = ix[k]; j = ix[N_+k]
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

