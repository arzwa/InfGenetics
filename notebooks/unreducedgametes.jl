
using Pkg; Pkg.activate("/home/arzwa/dev/InfGenetics/")
using Distributed, StatsBase, Plots, Serialization, Distributions, Optim
addprocs(8)
@everywhere using InfGenetics

# estimate parameter of a censored geometric distribution by ML
function logpdf_censored_geometric(p, c, x)
    x < c ? log(p) + (x-1)*log(1-p) : (c-1)*log(1-p)
end

function mle_censored_geometric(c, xs)
    fun = p -> -sum(logpdf_censored_geometric.(Ref(p), Ref(c), xs))
    optimize(fun, 0.0, 1.0).minimizer
end

us = 0.02:0.01:0.1
m  = 0.75
n  = 20
γ  = 0.125
z  = 1.5
sims = pmap(us) do u
    U = Umat(u, u)
    M = InfDemeMixEst(m=m, θ=fill(z, 3), U=U, γ=γ)
    ks = Int64[]
    K  = 0
    while length(ks) < n 
        est, ngen, c = simest2(M)
        K += 1
        if c[end] > c[1]
            push!(ks, K)
            K = 0
        end
    end
    @info u, ks
    u, ks
end

map(sims) do (u, ks)
    p = mle_censored_geometric(Inf, ks)
    u, log(p/(1-p))
end |> z->plot(z, ms=2, marker=true)
    
