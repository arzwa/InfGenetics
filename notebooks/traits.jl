using InfGenetics, Distributions, Plots, StatsPlots, Parameters
using InfGenetics: haploidgamete, diploidgamete
plotsdefault()

N = 1000

xs = map(-2:1.0:2) do z̄
    # set u=1, all offspring is tetraploid
    M = InfDemeMix(U=Umat(0.0, 0.0))
    pop = InfPop(z=rand(Normal(z̄, √(2M.V)), N), c=fill(2,N))
    pop2 = generation(pop, M)
    pop4 = generation(pop, reconstruct(M, U=Umat(1.0, 0.0)))
    z̄, pop2.z, pop4.z
end

map(xs) do (z̄, z2, z4)
    bins = -6:0.25:6
    stephist(z2, fill=true, fillalpha=0.2, label="diploids", normalize=true, bins=bins)
    stephist!(z4, fill=true, fillalpha=0.2, label="tetraploids", normalize=true, bins=bins)
    vline!([z̄], label="", legend=:outertopright, size=(400,200)) 
end |> x->plot(x..., size=(900,300))

function familydist(zi, zj, ki, kj, ci, cj, Fi, Fj, M)
    @unpack β = M
    k   = ki + kj
    yi  = (ki/ci)*(zi/β[ci-1])
    yj  = (kj/cj)*(zj/β[cj-1])
    yi, Vi = ki == 1 ? haploidgamete(M, zi, Fi, ci) : diploidgamete(M, zi, Fi, ci)
    yj, Vj = kj == 1 ? haploidgamete(M, zj, Fj, cj) : diploidgamete(M, zj, Fj, cj)
    z̄ij = (β[k-1]/β[ci-1])*yi + (β[k-1]/β[cj-1])*yj
    Vij = β[k-1]^2*(Vi + Vj)
    (z=z̄ij, V=Vij)
end

M = InfDemeMix(U=Umat(0.05, 0.05), α=[1/6, 1/6, 1/6])

familydist(0., 0., 2, 2, 2, 2, 0.0, 0.0, M)

# trait value of newly formed tetraploid, with diploid parents is Gaussian with
# mean `z̄ = (2zi + 2zj)/β₄` and variance `β₄²(Vik + Vjl)`.

# Source population (mean zero)
βs = [[1, 1, 1], .√([1,2/3,1/2]), [1, 2/3, 1/2]]
map(βs) do β
    M = InfDemeMix(U=Umat(0.05, 0.05), α=[1/6, 1/6, 1/6], β=β)
    plot()
    for k=2:4
        plot!(Normal(0, k*M.V*M.β[k-1]^2))
    end
    plot!()
end |> x-> plot(x..., layout=(3,1))



# Fitness on the island exp(γ*(z-θ))

γ = 0.125
θ = 2
P1 = plot(z->exp(γ*(z-θ)), xlim=(-2, 4),
    legend=false, ylabel="Fitness", 
    xlabel="Trait value (\$z\$)", color=:black)
vline!([θ], color=:black, ls=:dash)
hline!([1], color=:black, ls=:dash)
vline!([0], color=:lightgray, legend=false)
β = [1.0, 0.9, 0.85]
M = InfDemeMix(U=Umat(0.05, 0.05), α=[1/6, 1/6, 1/6], β=β)
P2 = plot(legend=:topright)
for k=2:4
    d = ["","di","tri","tetra"]
    plot!(Normal(0, k*M.V*M.β[k-1]^2), label="$(d[k])ploid")
end
plot!(xlabel="Trait value (\$z\$)", xlim=(-2,4), ylabel="pdf")
plot(P1,P2, layout=(2,1), size=(220,370))

pp = InfGenetics.cytotype_equilibrium(Umat(0.05,0.05))
P3 = bar(pp,
    yscale=:log10, color=:white, ylim=(0.001,1), 
    xticks=([1,2,3],["diploid","triploid","tetraploid"]), 
    ylabel="proportion (\$\\log\\pi\$)", legend=false)
map(1:3) do k
    annotate!(k, 0.025, text(@sprintf("\$%.1f\\%%\$", pp[k]*100), 9, :center))
end
plot!(size=(230,140))

plot(P3, P2, P1, layout=(1,3), size=(750,180), margin=4Plots.mm)

