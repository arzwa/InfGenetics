using Pkg; Pkg.activate(@__DIR__)

# some checks of infinitesimal predictions against an actual Mendelian system
using Random, StatsBase, Distributions
include("src/unlinkedloci.jl")

ϕ = [1., 0., 0., 0., 0.]
V = 1.2
v2 = map(1:1000) do _
    G = randgenome(100, 2, ϕ, V)
    var(map(_->sum(gamete(G, 0.)), 1:200))
end |> mean
v4 = map(1:1000) do _
    G = randgenome(100, 4, ϕ, V)
    var(map(_->sum(gamete(G, 0.)), 1:200))
end |> mean
v2, V/2/2, v4, V/2/2

# check whether the segregation variance matches predictions
ϕ = [0.1, 0.54, 0.07, 0.25, 0.04]
α = 0.3
V = 3  # phenotypic var, segvar=V/2, gameticvar = V/4
vs = map(1:1000) do _
    G = randgenome(100, 4, ϕ, V)
    var(map(_->sum(gamete(G, α)), 1:100))
end
mean(vs), Φ(ϕ) * (V/4) * (1+2α), (1 - F(ϕ)) * (V/4) * (1+2α)

# the segregation variance for unreduced gametes
ξ = 0.3
V = 1  # segvar V0 = V/2, gameticvar = V0/2 = V/4
ϕ = [0.8, 0.2]
vs = map(1:1000) do _
    G = randgenome(100, 2, ϕ, V)
    v1 = var(map(_->sum(ugamete(G, ξ)), 1:100))
    v2 = var(map(_->sum( gamete(G, 0)), 1:100))
    v1, v2
end
mean(first.(vs)), mean(last.(vs)), 4*(V/4)*ϕ[1]*ξ, ϕ[1]*V/4

# triploids
ξ = 0.3
V = 1.2  # segvar V0 = V/2, gameticvar = V/3 
ϕ = [.5, 0.3, 0.2]
f = ϕ[3] + ϕ[2]/3
vs = map(1:1000) do _
    G = randgenome(100, 3, ϕ, V)
    v1 = var(map(_->sum(tripdipgamete(G, ξ)), 1:100))
    v2 = var(map(_->sum(triphapgamete(G)), 1:100))
    v1, v2
end
vs2 = first.(vs)
vs1 = last.(vs)
mean(vs2), (4/3)*(V/2)*(1/3 + ξ) * (1-f)
# (4/3)(V0/3 + ξV0)(1-F)

mean(vs1), (4*(V/2)/9)* (1-f)
# 4V0/9 * (1-F)

var(map(_->sum(randgenome(100, 2, ϕ, V)), 1:10000))
var(map(_->sum(randgenome(100, 4, ϕ, V)), 1:10000))

# Individual based infinitesimal model simulations
N = 200
t = 2000
α = 0.3
sim1 = iterate(wfgen, InfDeme(randn(N), k=2, V=1/2), t, callback=x->popstats(x[2]))
sim2 = iterate(wfgen, InfDeme(randn(N), k=4, V=1/2), t, callback=x->popstats(x[2]))
sim3 = iterate(wfgen, InfDeme(randn(N), k=4, V=1/2, α=α), t, callback=x->popstats(x[2]))

default(fg_legend=:transparent, bg_legend=:transparent, framestyle=:box,
        gridstyle=:dot)
plot(size=(√2*250,250), ylabel="\$V_{z}\$", xlabel="\$t\$")
plot!(getindex.(sim1, 2), color=:gray, label="diploids", alpha=0.5)
plot!(x->exp(-x/(N*2)), color=:black, lw=2, label="")
plot!(getindex.(sim2, 2), color=:cornflowerblue, label="tetraploids \$\\alpha=0\$", alpha=0.5)
plot!(x->exp(-x/(N*4)), color=:black, ls=:dash, lw=2, label="")
plot!(getindex.(sim3, 2), color=:salmon, label="tetraploids \$\\alpha = 0.3\$", alpha=0.5)
#savefig("doc/img/124.pdf")

plot(log.(getindex.(sim1, 2)))
plot!(log.(getindex.(sim2, 2)))
plot!(log.(getindex.(sim3, 2)))
