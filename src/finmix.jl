# Finite number of unlinked loci under Mendelian inheritance
# Use the same model object as the infinitesimal model (`InfDemeMix`)

"""
    offspring(p1, p2, g1, g2, model)

Generate an offspring genotype from parents `p1` and `p2`, contributing a
gamete of ploidy level `g1` and `g2` respectively.
"""
function offspring(p1, p2, g1, g2, model)
    @unpack α = model
    L, c1 = size(p1)
    L, c2 = size(p2)
    c = g1 + g2
    x = Matrix{Float64}(undef, L, c)
    @inbounds for i=1:L
        segregate!(x, p1, g1, i,    1, α[c1-1]) 
        segregate!(x, p2, g2, i, g1+1, α[c2-1]) 
    end
    return x
end

function segregate!(x, p, g, i, o, α)
    if rand() < α 
        x[i,o:o+g-1] .= rand(p[i,:])
    else
        x[i,o:o+g-1] .= sample(p[i,:], g, replace=false)
    end
end

struct FinPop{T}
    z::Vector{T}
    g::Vector{Matrix{T}}
    c::Vector{Int}
end

Base.similar(X::FinPop) = FinPop(similar(X.z), similar(X.g), similar(X.c))
popsize(m::FinPop) = length(m.z)

function gamete(M::InfDemeMix, c)
    u = M.U[c-1,1]/(M.U[c-1,1] + M.U[c-1,2])
    rand() < u ? 1 : 2
end

function generation(pop::FinPop, M::InfDemeMix, fitness::Function)
    @unpack α, β, w, U = M
    N = popsize(pop)
    # fertilities, i.e. contribution to gamete pool
    ws = InfGenetics.fitnesses(pop, U, w, fitness)
    idx = sample(1:N, Weights(ws), 2N, replace=true) 
    pop_ = similar(pop)
    for k=1:N
        i = idx[k]; j = idx[N+k]
        gi = gamete(M, pop.c[i])
        gj = gamete(M, pop.c[j])
        gij = offspring(pop.g[i], pop.g[j], gi, gj, M)
        cij = gi + gj
        pop_.g[k] = gij
        pop_.c[k] = cij
        pop_.z[k] = sum(gij) * β[cij-1]
    end
    return pop_
end

