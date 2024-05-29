using InfGenetics, Plots, PlotThemes, StatsBase, Parameters; theme(:hokusai)
using Distributions, ThreadTools

function simest(M, minn=100, maxgen=Inf)
    pops = [InfPop(z=zeros(0))]
    t = 0
    while true
        t += 1
        pop = generation(pops[end], M)
        push!(pops, pop)
        (length(pop.z) > minn || t > maxgen) && (return t, pops)
    end
end

function simuntil(simfun, condition)
    while true
        x = simfun()
        condition(x) && return x
    end
end

function ρmat(a)
    ρ = [1.0 a a; a 1.0 a; a a 1.0]
    ρ ./= mean(ρ)
end

# As in Abu Awad
# assume loss of self-incompatibility
M = InfDemeMixEstSelf(σ=[0.0, 0.5, 0.5])
t, pops = simest(M); @info t
plot(permutedims(hcat(map(x->counts(x.c, 2:4), pops[end-30:end])...)))

xs = map(1:1000) do _
    t, pops = simest(M)
    x = (t, counts(pops[end].c, 2:4))
    @info x; x
end

filter(x->x[3]>x[1], last.(xs))

nrep = 10000
θs = -2.5:0.5:-1.5
ss = 0:0.2:1
res = map(θs) do θ
    map(ss) do σ
        @info σ, θ
        M = InfDemeMixEstSelf(σ=[0.0, σ, σ], θ=fill(θ, 3))
        map(1:nrep) do _
            t, pops = simest(M)
            counts(pops[end].c, 2:4)
        end
    end
end

y = map(res) do x
    map(xs->length(filter(x->x[3]>x[1], xs))/nrep, x)
end
plot(ss, y, 
    label=reshape(map(x->"\$\\theta=$x\$", -3:0.5:-1.5), 1, 4), 
    marker=true, xlabel="\$\\sigma\$", ylabel="\$P\$", size=(320,260), 
    legend=:topleft)

function tetswin(pop)
    xs = counts(pop.c, 2:4)
    xs[3] > xs[1]
end

function extract(pops)
    (N = map(x->length(x.z), pops),
     zs= map(x->x.z, pops),
     z = map(x->mean(x.z), pops),
     vz= map(x->var(x.z), pops),
     F = map(x->mean(x.F), pops),
     c = permutedims(mapreduce(x->counts(x.c, 2:4), hcat, pops)))
end


ex1 = simuntil(()->simest(InfDemeMixEstSelf(σ=[0.0, 1/3, 1/3])), x->tetswin(x[2][end]))
ex2 = simuntil(()->simest(InfDemeMixEstSelf(σ=[0.0, 1/3, 1/3])), x->!tetswin(x[2][end]))

x1 = extract(ex1[2])
x2 = extract(ex2[2])

plot(x1.F[end-20:end,:])
plot!(x2.F[end-20:end,:])



# ----------
#
M = InfDemeMixEstPair(σ=[0.0, 0.5, 0.5], θ=fill(-3.0,3), m=1.0)
t, pops = simest(M); @info t

nrep = 1000
θs = -2.0:1.0:-1.0
ss = 0:0.20:1
ms = [1.0, 5.0, 10.0]
res = map(ms) do m
    map(θs) do θ
        tmap(ss) do σ
            @info σ, θ, m 
            M = InfDemeMixEstPair(σ=[0.0, σ, σ], θ=fill(θ, 3), m=m)
            map(1:nrep) do _
                t, pops = simest(M)
                counts(pops[end].c, 2:4)
            end
        end
    end
end
    

xs = map(1:1000) do _
    t, pops = simest(M)
    x = (t, counts(pops[end].c, 2:4))
    @info x; x
end

map(enumerate(res)) do (k,xx)
    y = map(enumerate(xx)) do (k,x)
        map(xs->length(filter(x->x[3]>x[1], xs))/nrep, x)
    end
    plot(ss, y, color=[1 2 3],
        label=reshape(map(x->"\$\\theta=$x\$", -3:0.5:-1.5), 1, 4), 
        marker=true, xlabel="\$\\sigma\$", ylabel="\$P\$", size=(320,260), 
        legend=:topleft, title="\$M=$(ms[k])\$")
end |> x->plot(x..., layout=(1,3), size=(700,200), ylim=(0,1), margin=5Plots.mm)

savefig(joinpath(imgpth(), "selftest2.svg")) |> mdfig


