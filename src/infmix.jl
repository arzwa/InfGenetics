"""
    InfDemeMix

Parameterization for a single deme evolving according to the mixed-ploidy
infinitesimal model. This assumes a diploid-triploid-tetraploid complex, where
only haploid and diploid gametes are formed. 

- `U` is the cytotype × gametic type matrix
- `w` are the intrinsic fitnesses of each cytotype
- `β` are the allelic effect scalers across cytotypes
- `α` are the probabilities that a diploid gamete from each cytotype contains
   two copies of the same parental allele (e.g. double reduction in tetraploids)
- `V` is the genetic variance associated with a haplotype in the hypothetical
   diploid reference population

`U` can be seen as capturing the fertility component of intrinsic (cytotype
specific) fitness, whereas `w` captures the viability component. 
"""
@with_kw struct InfDemeMix{T}
    U::Matrix{T} = [0.95 0.05 ; 0.05 0.05 ; 0.0 0.95]  # cytotype × gamete 
    w::Vector{T} = ones(3)              # cytotype specific viability 
    β::Vector{T} = [1, √(2/3), √(1/2)]  # variance scaler, will be in [0,1]
    α::Vector{T} = zeros(3)             # double reduction rate
    V::T = 0.5                          # Vx (haplotype variance in diploids)
end

"""
    InfPop

A mixed-ploidy population with trait values `z`, inbreeding coefficients `F`,
coancestry coefficients `Φ` and ploidy levels `c`.
"""
@with_kw struct InfPop{T}
    z::Vector{T}   = Float64[] # trait values (scalars, for now)
    F::Vector{T}   = zeros(length(z)) # inbreeding
    Φ::Matrix{T}   = diagm(fill(0.5, length(z))) # coancestry
    c::Vector{Int} = fill(2, length(z)) # ploidy levels
end

popsize(M::InfPop) = length(M.z)
Base.similar(X::InfPop) = InfPop(
    similar(X.z), similar(X.F), similar(X.Φ), similar(X.c))
Base.getindex(X::InfPop, i) = (z=X.z[i], F=X.F[i], c=X.c[i])

# We can have selection on parents, or as in B&E selection on offspring when we
# assume directional selection. The latter is mostly relevant for studying
# establishment etc.

# Let's start by implementing WF-like dynamics for the mixed-ploidy model.
# But note that the variance will erode if there's no mutation...
function fitnesses(pop, U, w, fitfun)
    [(U[k-1,1] + U[k-1,2])*w[k-1]*fitfun(z) for (z,k) in zip(pop.z,pop.c)]  
end

function generation(pop::InfPop, M::InfDemeMix, fitness::Function=z->1.0)
    @unpack F, Φ = pop
    @unpack α, β, w, U = M
    N = popsize(pop)
    # fertilities, i.e. contribution to gamete pool
    ws = fitnesses(pop, U, w, fitness)
    idx = sample(1:N, Weights(ws), 2N, replace=true) 
    pop_ = similar(pop)
    P = zeros(N,N)
    for k=1:N
        i = idx[k]; j = idx[N+k]
        zi, Fi, ci = pop[i]
        zj, Fj, cj = pop[j]
        gi, yi, Vi = gamete(M, zi, Fi, ci)
        gj, yj, Vj = gamete(M, zj, Fj, cj)
        cij = gi + gj
        βij = β[cij-1]
        βi  = β[ci-1]
        βj  = β[cj-1]
        Yi  = rand(Normal(yi*βij/βi, √(βij^2*Vi))) 
        Yj  = rand(Normal(yj*βij/βj, √(βij^2*Vj))) 
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

# this is conditional on a succesful gamete
function gamete(M, z, F, c)
    u = M.U[c-1,1]/(M.U[c-1,1] + M.U[c-1,2])
    rand() < u ? 
        (1, haploidgamete(M, z, F, c)...) : 
        (2, diploidgamete(M, z, F, c)...)
end

# Note that we can rescale only once we know the offspring cytotype, so here
# everything is wrt to diploid reference pop.
function haploidgamete(M, z, F, c)
    if c == 2
        z/2, (1-F)*M.V/2
    elseif c == 3
        #z/3, (1-F)*M.V/2  # XXX had this before, but wrong?  
        z/3, 2*(1-F)*M.V/3
    else
        @error "Shouldn't happen"
    end
end

function diploidgamete(M, z, F, c)
    if c == 2
        z, 2M.α[1]*(1-F)*M.V
    elseif c == 3
        2z/3, (2/3)*(1+3M.α[2])*(1-F)*M.V
    elseif c == 4
        z/2, (1+2M.α[3])*(1-F)*M.V
    else
        @error "Shouldn't happen"
    end
end

function neutral_noinbreeding_generation(
        pop::InfPop, M::InfDemeMix)
    @unpack F, Φ = pop
    @unpack α, β, w, U = M
    N = popsize(pop)
    ws = fitnesses(pop, U, w, z->1.0)
    ws ./= mean(ws)
    idx = sample(1:N, Weights(ws), 2N, replace=true) 
    pop_ = similar(pop)
    for k=1:N
        i = idx[k]; j = idx[N+k]
        zi, Fi, ci = pop[i]
        zj, Fj, cj = pop[j]
        gi, yi, Vi = gamete(M, zi, Fi, ci)
        gj, yj, Vj = gamete(M, zj, Fj, cj)
        cij = gi + gj
        βij = β[cij-1]
        βi  = β[ci-1]
        βj  = β[cj-1]
        Yi  = rand(Normal(yi*βij/βi, √(βij^2*Vi))) 
        Yj  = rand(Normal(yj*βij/βj, √(βij^2*Vj))) 
        pop_.z[k] = Yi + Yj
        pop_.c[k] = cij
    end
    return pop_
end

