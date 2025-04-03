using Pkg; Pkg.activate("/home/arzwa/dev/InfGenetics/")
using Distributed, StatsBase, Plots, Serialization, Distributions, Optim
using DataFrames, CSV
addprocs(6)
@everywhere using InfGenetics

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

# Basic model
# ===========
nrep = 50_000
ms   = 10 .^ range(log10(0.01), stop=log10(3), length=15)
zs   = [3.0, 2.5, 2.0, 1.5]
γ    = 0.25
U    = Umat(0.05, 0.05)
nmax = 10000
Nest = 100
M  = InfDemeMixEst(m=0.0, θ=fill(2.5, 3), U=U, γ=γ)

# baseline
base = pmap(zs) do z
    @info z
    M  = InfDemeMixEst(m=0.0, θ=fill(z, 3), U=U, γ=γ)
    nr = 2_500_000
    x2 = map(_->simest2(M, x0=InfPop(z=randn(1), c=[2]), Nest=Nest), 1:nr)
    x3 = map(_->simest2(M, x0=InfPop(z=randn(1), c=[3]), Nest=Nest), 1:nr)
    x4 = map(_->simest2(M, x0=InfPop(z=randn(1), c=[4]), Nest=Nest), 1:nr)
    p22 = sum(first.(x2) .&& map(x->argmax(x[3]) == 1, x2)) / nr
    p24 = sum(first.(x2) .&& map(x->argmax(x[3]) == 3, x2)) / nr
    p32 = sum(first.(x2) .&& map(x->argmax(x[3]) == 1, x2)) / nr
    p34 = sum(first.(x2) .&& map(x->argmax(x[3]) == 3, x2)) / nr
    p4  = sum(first.(x4) .&& map(x->argmax(x[3]) == 3, x4)) / nr
    (z=z, p2=p2, p4=p4, p32=p32, p34=p34)
end

map(base) do (z, p2, p4, p32, p34)
    Z = sum([p2, p32 + p34, p4] .* M.c)
    (z, (p4 * M.c[3] + p34 * M.c[2])/Z)
end
    
M = InfDemeMixEst(m=0.001, θ=fill(2.5, 3), U=U, γ=0.25)
x = map(_->simest2(M, x0=InfPop(z=Float64[], c=Int64[]), Nest=100), 1:10000)
filter(y->argmax(y[3])==3, x) 


sims1 = pmap(Iterators.product(ms, zs)) do (m, z)
    @info (m, z)
    model = InfDemeMixEst(m=m, θ=fill(z, 3), U=U)
    reps = map(1:nrep) do _
        sim = simest2(model, nmax=nmax, Nest=Nest)
    end
    InfGenetics.parameters(model), reps
end


df1 = make_a_df(sims1, nmax)
CSV.write("data/sims1b-df-z3.csv", df1)


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
CSV.write("data/sims2b-df-z3.csv", df2)


# With selfing
# =============
# This is not just with selfing, but with loss of SI...
nrep = 25_000
zs   = 2.0
σs   = [[NaN, s, s] for s=0:0.2:1.0]

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
CSV.write("data/sims3b-df.csv", df3)


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
CSV.write("data/sims4b-df.csv", df4)

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
CSV.write("data/sims5b-df.csv", df5)
