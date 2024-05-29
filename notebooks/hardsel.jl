using InfGenetics, Plots, PlotThemes, StatsBase, Parameters; theme(:hokusai)
using Distributions

N = 500
M = InfDemeHS()

# equilibrate
xs = let pop = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
    map(1:2000) do i
        pop = generation(pop, M)
        pop.F .= 0.0
        pop.Φ .= 0.0
        pop
    end
end

# environmental change
ys = let M = InfDemeHS(β=ones(3))
    os = map(x->4*sin(0.05x), 0:251)
    map(1:100) do j
        pop = deepcopy(xs[end])
        res = map(enumerate(os)) do (i,o)
            M = reconstruct(M, θ=fill(o, 3))
            pop = generation(pop, M)
        end
        cs = counts(pop.c, 1:4)
        @info j, cs
        (sum(cs) > 0, cs[4] > cs[2], res)
    end
end;

y = ys[10][3]
P1 = plot(permutedims(hcat(map(x->counts(x.c, 2:4), y)...)), ylabel="\$N\$", label=["diploid" "triploid" "tetraploid"], legend=:topleft)
P2 = plot(map(x->5*sin(0.05x), 0:251), ylabel="\$z\$", label="\$\\theta\$", color=:gray, legend=:topright)
plot!(map(x->mean(x.z), y), ribbon=map(x->2std(x.z), y), color=2, label="\$\\bar{z}\$")
plot(P1, P2, size=(550,200), layout=(1,2), xlabel="generation", margin=3Plots.mm)


