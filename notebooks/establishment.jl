# Consider first, as in Barton & Etheridge 2018, for a given source population
# the probability of a single diploid/tetraploid migrant establishing.
using Distributed, Plots
#addprocs(6)
@everywhere using Parameters, InfGenetics
@everywhere M = InfDemeMixEst(U=Umat(0.0,0.0), γ=0.25, m=0.0, θ=[2.0,2.0,2.0], V=0.5)
gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(2), length=10)]

res = pmap(gs) do g
    nr = 500_000
    x2 = map(_->simest2(reconstruct(M, γ=g), x0=InfPop(z=[0.0], c=[2])), 1:nr)
    x4 = map(_->simest2(reconstruct(M, γ=g), x0=InfPop(z=[0.0], c=[4])), 1:nr)
    (sum(first.(x2)) / nr, sum(first.(x4)) / nr)
end

using Serialization
#serialize("data/est1.jls", res)
res = deserialize("data/est1.jls")

p2 = first.(res)
p4 = last.(res)
x2 = floor.(Int, p2 .* 500000)
x4 = floor.(Int, p4 .* 500000)
e2 = jeffreys_interval.(x2, 500000, 0.01)
e4 = jeffreys_interval.(x4, 500000, 0.01)

P1 = plot(gs, p2, marker=true, ms=2, label="diploid",
    ribbon=(p2 .- first.(e2), last.(e2) .- p2), fillalpha=0.4)
plot!(gs, p4 , marker=true, ms=2, label="tetraploid",
    ribbon=(p4 .- first.(e4), last.(e4) .- p4), fillalpha=0.4)
plot!(title="(A) \$\\theta=2, \\beta_4 = \\sqrt{1/2}\$", 
    ylabel="\$P\$", xlabel="\$\\gamma\$", ylim=(0,0.01))

bs = 0.25:0.05:1.0
resb = pmap(bs) do b
    nr = 500_000
    x4 = map(_->simest2(
        reconstruct(M, β=[1.0, 1.0, √b]), 
        x0=InfPop(z=[0.0], c=[4])), 1:nr)
    sum(first.(x4)) / nr
end

#serialize("data/est1b.jls", resb)
resb = deserialize("data/est1b.jls")
    
#ps = map(_->simest2(M, x0=InfPop(z=[0.0], c=[2])), 1:500_000)
#p2 = sum(first.(ps)) / length(ps)
p2 = 0.001172

P2 = scatter(.√bs , first.(resb) ./ p2, legend=false, color=:black, ms=2.5,
		ylabel="\$P_4/P_2\$", xlabel="\$\\beta_4^2\$", 
        ylim=(0,12.5), size=(260,210),
        #grid=true, 
        yticks=0:1:12, xlim=(0.48,1.02), xticks=0.5:0.1:1,
        title="(B) \$\\theta=2, \\gamma=0.25\$")
hline!([1], color=:firebrick, alpha=0.4)
vline!([√(0.5)], color=:gray, ls=:dash)

plot(P1, P2, size=(260*2,210), margin=3Plots.mm)
#savefig("doc/img/est1b.pdf")

nr = 500_000
V2 = 2M.β[1]^2 * M.V 
V4 = 4M.β[3]^2 * M.V 
x2 = map(_->simest2(M, x0=InfPop(z=[rand(Normal(0,√V2))], c=[2])), 1:nr)
x4 = map(_->simest2(M, x0=InfPop(z=[rand(Normal(0,√V4))], c=[4])), 1:nr)
(sum(first.(x2)) / nr, sum(first.(x4)) / nr)

# Full space...
gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(2), length=10)]
bs = 0.25:0.05:1.0
gb = Iterators.product(gs, bs)

res = pmap(gb) do (g, b)
    @info (g, b)
    nr = 500_000
    x2 = map(_->simest2(reconstruct(M, γ=g), 
        x0=InfPop(z=[0.0], c=[2])), 1:nr)
    x4 = map(_->simest2(reconstruct(M, γ=g, β=[1.0, 1.0, √b]), 
        x0=InfPop(z=[0.0], c=[4])), 1:nr)
    (g, b, sum(first.(x2)) / nr, sum(first.(x4)) / nr)
end

#serialize("data/est1-2d.jls", res)

resc = deserialize("data/est1-2d.jls")

rp = map(resc[2:end,:]) do (g, b, p2, p4)
    log10(p4/p2)
end 

cc = maximum(abs, filter(isfinite, rp))
P3 = heatmap(bs, log10.(gs[2:end]), rp, color=:balance, clims=(-cc, cc), 
#    xticks=(1:2:length(bs), bs[1:2:end]),
#    yticks=(1:length(gs)-1, map(x->@sprintf("%.2f", x), gs[2:end]))
    ylabel="\$\\log_{10}\\gamma\$",
    xlabel="\$\\beta_4\$",
    margin=4Plots.mm,
    title="(C) \$\\log_{10}\\left(P_4/P_2\\right)\$",
    size=(320,280)
)
vline!([√1/2], color=:black, lw=2, legend=false)

plot(P1, P2, P3, size=(800,200), margin=5Plots.mm,
    layout=grid(1,3,widths=[0.3,0.3,0.4]))
savefig("doc/img/est1b.pdf")
