
using InfGenetics, StatsBase, Distributions, Parameters, Plots

MM = InfDemeMix(U=Umat(0.0, 0.0), V=0.5)
N  = 1000
map([2,4]) do c
    pop = InfPop(z=zeros(N), c=fill(c,N))
    for i=1:20
        pop = generation(pop, MM)
        pop.F .= 0.0
    end
    mean(pop.z), var(pop.z)
end



# have a closer look at the trait distributions
MM = InfDemeMix(U=Umat(0.05, 0.05), V=0.5, β=[1., 1, 1], α=maxα)
N  = 20000
pop = InfPop(z=randn(N), c=fill(2,N))
Zs  = [Float64[] for i=1:3]
for i=1:200
    pop = InfGenetics.neutral_noinbreeding_generation(pop, MM)
    pop.F .= 0.0
    (i % 1000 == 0 && @info i)
    for k=2:4
        length(Zs[k-1]) > 100000 && continue
        zs = pop.z[pop.c .== k]
        Zs[k-1] = [Zs[k-1] ; zs]
    end
end
# seems that the higher ploidy levels converge very slowly to the equilibrium?
plot()
for Z in Zs
    stephist!(Z, bins=(-6:0.2:6),normalize=true)
end
plot!()

Vpred_ = MM.β .^ 2 .* (2:4) .* MM.V
Vpred  = MM.β .^ 2 .* (2:4) .* MM.V
Vpred[2] *= (1 + (2/3)*MM.α[1])
Vpred[3] *= (1 + MM.α[1])
map(1:3) do k
    stephist(Zs[k], normalize=true, color=:black, 
        legend=k==3 ? :outertopright : false,
        title="$(strs[k+1])", label="simulation")
    plot!(Normal(0,1), fill=true, color=:gray, alpha=0.2, label="Normal(0,2V)")
    plot!(Normal(0, √(Vpred[k])), label="approximate eq.")
    plot!(Normal(0, √(Vpred_[k])), label="monocytotypic eq.")
end |> x->plot(x..., size=(800,180), xlabel="Trait value \$z\$", 
    layout=grid(1,3,widths=[0.26,0.26,0.48]), margin=4Plots.mm)

# only newly formed
n = 1000000
zs = map(1:n) do _
    zi = randn()
    zj = randn()
    Vi = 2MM.α[1]*MM.V
    Vj = 2MM.α[1]*MM.V
    zij = MM.β[3]*(zi + zj)
    vij = MM.β[3]^2*(Vi + Vj)
    rand(Normal(zij, √vij))
end
stephist(zs, norm=true)
vv = 4MM.β[3]^2*MM.V*(1+MM.α[1])
plot!(Normal(0, √vv))


