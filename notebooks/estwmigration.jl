using Pkg; Pkg.activate("/home/arzwa/dev/InfGenetics/")
using Distributed, StatsBase, Plots, Serialization, Distributions, Optim
using DataFrames, CSV
addprocs(6)
@everywhere using InfGenetics, Distributions

# estimate parameter of a censored geometric distribution by ML
function logpdf_censored_geometric(p, c, x)
    x < c ? log(p) + (x-1)*log(1-p) : (c-1)*log(1-p)
end

function mle_censored_geometric(c, xs)
    fun = p -> -sum(logpdf_censored_geometric.(Ref(p), Ref(c), xs))
    optimize(fun, 0.0, 1.0).minimizer
end

function make_a_df(sims, nmax)
    map(eachcol(sims)) do col
        map(col) do xs 
            ps, reps = xs
            ts = getindex.(reps, 2)
            tt = 1/mle_censored_geometric(nmax, ts)
            nsucc = sum(first.(reps))
            xx = sum(map(x->x[3][3] > x[3][1] && x[1], reps))
            p = xx/nsucc
            e1,e2 = jeffreys_interval(xx, nsucc, 0.05)
            merge(ps, (nsucc=nsucc, n=xx, t=tt, p=p, p1=p-e1, p2=e2-p))
        end 
    end |> x->vcat(x...) |> DataFrame
end

# Basic model (m -> 0 limit)
# ==========================
nrep = 500_000
Nest = 100
zs   = [2.5, 2.0, 1.5]
γ    = 0.25
U    = Umat(0.05, 0.05)

# baseline
base = pmap(zs) do z
    @info z
    M  = InfDemeMixEst(m=0.0, θ=fill(z, 3), α=maxα, U=U, γ=γ)
    zsource = Normal(0,1)
    x2 = map(_->simest2(M, x0=InfPop(z=rand(zsource, 1), c=[2]), Nest=Nest), 1:nrep)
    x3 = map(_->simest2(M, x0=InfPop(z=rand(zsource, 1), c=[3]), Nest=Nest), 1:nrep)
    x4 = map(_->simest2(M, x0=InfPop(z=rand(zsource, 1), c=[4]), Nest=Nest), 1:nrep)
    (z, x2, x3, x4)
end

M  = InfDemeMixEst(m=0.0, θ=fill(2.0, 3), U=U, γ=γ)
m0 = map(base) do y
    p4 = map(y[2:end]) do xs
        length(filter(x->x[1] && (argmax(x[3]) == 3), xs))/nrep
    end
    p2 = map(y[2:end]) do xs
        length(filter(x->x[1] && (argmax(x[3]) == 1), xs))/nrep
    end
    y[1], sum(p4 .* M.c)/(sum(p2 .* M.c) + sum(p4 .* M.c))
end

# very weak migration, for comparison
res_ = pmap(zs) do z
    M = InfDemeMixEst(m=1e-6, θ=fill(z, 3), U=U, γ=0.25)
    x = map(_->InfGenetics.simest3(M, Nest=100), 1:nrep)
    x4 = filter(y->argmax(y[3])==3, x) 
    z, length(x4)/nrep
end


# With recurrent migration
# ========================
nrep = 50_000
ms   = 10 .^ range(log10(0.01), stop=log10(3), length=15)
zs   = [2.5, 2.0, 1.5]
γ    = 0.25
U    = Umat(0.05, 0.05)
nmax = 10_000
Nest = 100

sims1 = pmap(Iterators.product(ms, zs)) do (m, z)
    @info (m, z)
    model = InfDemeMixEst(m=m, θ=fill(z, 3), U=U)
    reps = map(1:nrep) do _
        sim = simest2(model, nmax=nmax, Nest=Nest)
    end
    InfGenetics.parameters(model), reps
end

df1 = make_a_df(sims1, nmax)
CSV.write("data/sims1c_.csv", df1)


# With double reduction etc. (maxα)
# =================================
sims2 = pmap(Iterators.product(ms, zs)) do (m, z)
    @info (m, z)
    model = InfDemeMixEst(m=m, θ=fill(z, 3), U=U, α=maxα)
    reps = map(1:nrep) do _
        sim = simest2(model, nmax=nmax, Nest=Nest)
    end
    InfGenetics.parameters(model), reps
end

df2 = make_a_df(sims2, nmax)
CSV.write("data/sims2c_.csv", df2)


# With selfing
# =============
# This is not just with selfing, but with loss of SI...
nrep = 25_000
zs   = 2.0
ms   = 10 .^ range(log10(0.01), stop=log10(3), length=15)
σs   = [[NaN, s, s] for s=0:0.2:1.0]
nmax = 10000
Nest = 100

sims3 = pmap(Iterators.product(ms, σs)) do (m, σ)
    @info (m, σ)
    model = InfDemeMixEst(m=m, θ=fill(zs, 3), U=U, σ=σ)
    reps = map(1:nrep) do _
        sim = simest2(model, nmax=nmax, Nest=Nest)
    end
    InfGenetics.parameters(model), reps
end

#sims3 = deserialize("data/_sims3.jls")

df3 = make_a_df(sims3, nmax)
CSV.write("data/sims3c_.csv", df3)


# With increased selfing (but no SI in diploids)
# ================================================
# random selfing assumed in diploids
nrep = 25_000
zs   = 2.0
stet = 0:0.2:1.0
σs   = [[0.0, st, st] for st=stet]
nmax = 10000
Nest = 100

sims4 = pmap(Iterators.product(ms, σs)) do (m, σ)
    @info (m, σ)
    model = InfDemeMixEst(m=m, θ=fill(zs, 3), U=U, σ=σ)
    reps = map(1:nrep) do _
        sim = simest2(model, nmax=nmax, Nest=Nest)
    end
    InfGenetics.parameters(model), reps
end

df4 = make_a_df(sims4, nmax)
CSV.write("data/sims4c_.csv", df4)

# With assortative mating
# =======================
nrep = 25_000
zs   = 2.0
ρs   = [fill(r, 3) for r=0:0.2:1.0]
γ    = 0.25
U    = Umat(0.05, 0.05)
nmax = 10000
Nest = 100

sims5 = pmap(Iterators.product(ms, ρs)) do (m, ρ)
    @info (m, ρ)
    model = InfDemeMixEst(m=m, θ=fill(zs, 3), U=U, ρ=ρ)
    reps = map(1:nrep) do _
        sim = simest2(model, nmax=nmax, Nest=Nest)
    end
    InfGenetics.parameters(model), reps
end

df5 = make_a_df(sims5, nmax)
CSV.write("data/sims5c_.csv", df5)
