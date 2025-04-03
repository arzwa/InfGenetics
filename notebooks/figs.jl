using CSV, DataFrames, Plots, Serialization, InfGenetics
plotsdefault()

ms = 10 .^ range(log10(0.01), stop=log10(3), length=15)

df1 = CSV.read("data/sims1b-df.csv", DataFrame)
df2 = CSV.read("data/sims2b-df.csv", DataFrame)
sort!(df1, [:z2, :m], rev=true)
sort!(df2, [:z2, :m], rev=true)

#base = @NamedTuple{z::Float64, p2::Float64, p4::Float64}[(z = 2.0, p2 = 6.46e-5, p4 = 0.0014164), (z = 1.5, p2 = 0.0004414, p4 = 0.0046594), (z = 1.0, p2 = 0.0023506, p4 = 0.0139396)]
base = deserialize("data/baseline3.jls")

P1 = plot(title="(A) Time to establishment")
for gdf in groupby(df1, :z2)
    plot!(gdf[:,:m], gdf[:,:t], marker=true, ms=2, label="\$\\theta=$(gdf[1,:z2])\$")
end
P1 = plot(P1, marker=true, ms=2, yscale=:log10, xlabel="\$m\$",
    ylabel="\$\\overline{T}\$", xscale=:log10, legend=:bottomleft, 
    ylim=(10,500_000), yticks=[10,100,1000,10000,100000])
P2 = plot(title="(B) Tetraploid establishment", xscale=:log10)
for (df, ls, label) in zip([df1, df2], [:solid, :dot], [true, false])
    for (i,gdf) in enumerate(groupby(df, [:z2]))
        gdf = sort(gdf, :m)
        plot!(gdf[:,:m], gdf[:,:p], ribbon=(gdf[:,:p1], gdf[:,:p2]),
            marker=true, ms=2, ls=ls, color=i,
            label=label ? "\$\\theta=$(gdf[1,:z2])\$" : "", fillalpha=0.2)
    end
end
P2 = plot(P2, marker=true, ms=2, xlabel="\$m\$",
    ylabel="\$P_4\$", xscale=:log10, legend=false)
ppi = InfGenetics.cytotype_equilibrium(Umat(0.05,0.05))
hline!(ppi[[3]], label="",
    color=:black, xlabel="\$m\$", ylabel="\$P_4\$", ls=:dash)
for i=1:4
    px = (ppi[3]*base[i].p4)/(ppi[1]*base[i].p2 + ppi[3]*base[i].p4)
    @info px
    hline!([px ], color=i, label="")
end
plot!(bg_legend=:transparent, ylim=(0,0.17))
plot(P1, P2, size=(275*2,230), margin=3Plots.mm, titlefont=9)
savefig("doc/img_/fig3.pdf")


# Selfing/assortative mating
df3 = CSV.read("data/sims3b-df.csv", DataFrame)
df4 = CSV.read("data/sims4b-df.csv", DataFrame)
df5 = CSV.read("data/sims5b-df.csv", DataFrame)
sort!(df3, [:s4, :m])
sort!(df4, [:s4, :m])
sort!(df5, [:r4, :m])

P3 = plot(title="(A)", xscale=:log10)
for (df, ls, alpha, label) in zip([df3, df4], [:solid, :dot], [1, 0.7], [true, false])
    for (i,gdf) in enumerate(groupby(df, [:s4]))
        plot!(gdf[:,:m], gdf[:,:p], alpha=alpha,
            marker=true, ms=2, ls=ls, color=i,
            label=label ? "\$\\sigma=$(gdf[1,:s4])\$" : "")
    end
end
P3 = plot(P3, marker=true, ms=2, xlabel="\$m\$", ylim=(0,1),
    ylabel="\$P_4\$", xscale=:log10, legend=:outertopright)
P4 = plot(title="(B)", xscale=:log10)
for (i,gdf) in enumerate(groupby(df5, :r4))
    plot!(gdf[:,:m], gdf[:,:p], ylim=(0,0.5),
        marker=true, ms=2, color=i,
        label="\$\\rho=$(gdf[1,:r4])\$")
end
P4 = plot(P4, marker=true, ms=2, xlabel="\$m\$",
    ylabel="\$P_4\$", xscale=:log10, legend=:outertopright)
plot(P3, P4, size=(300*2,210), margin=3Plots.mm, titlefont=8)

savefig("doc/img_/fig5.pdf")
