# Consider first, as in Barton & Etheridge 2018, for a given source population
# the probability of a single diploid/tetraploid migrant establishing.
using Distributed, Plots, Serialization, Parameters
addprocs(8)
@everywhere using InfGenetics, Parameters

# Define the model
M = InfDemeMixEst(
    U=Umat(0.0,0.0),  # no unreduced gametes 
    γ=0.25,           # selection strength
    m=0.0,            # no migration
    θ=fill(2.5, 3),   # maladaptation
    V=0.5)            # segregation variance

ps = pmap(_->simest2(M, x0=InfPop(z=[0.0], c=[2])), 1:1_000_000)
pp = sum(first.(ps)) / length(ps)
#pp = 0.002928


# 1. Examine the effect of selection strength (γ)
# -----------------------------------------------
# Define a range of selection strengths to examine
gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(2), length=15)]

res = pmap(gs) do g
    @info g
    nr = 1_000_000  # number of replicate simulations
    M = InfDemeMixEst(
        U=Umat(0.0,0.0),  # no unreduced gametes 
        γ=g,              # selection strength
        m=0.0,            # no migration
        θ=fill(2.5, 3),   # maladaptation
        V=0.5)            # segregation variance
    x2 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[2])), 1:nr)
    x4 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[4])), 1:nr)
    (g, sum(first.(x2)) / nr, sum(first.(x4)) / nr)
end

# Store the results
serialize("data/est1-t2.5.jls", res)

# Make a plot
res0 = deserialize("data/est1-t2.jls")
res1 = deserialize("data/est1-t2.5.jls")
res2 = deserialize("data/est1-t3.jls")

P1 = plot()
map([res1]) do res
    nr = 1_000_000
    p2 = getindex.(res, 2)
    p4 = getindex.(res, 3)
    x2 = floor.(Int, p2 .* nr)
    x4 = floor.(Int, p4 .* nr)
    e2 = jeffreys_interval.(x2, nr, 0.01)
    e4 = jeffreys_interval.(x4, nr, 0.01)
    plot!(gs, p2, marker=true, ms=2, label="diploid",
        ribbon=(p2 .- first.(e2), last.(e2) .- p2), fillalpha=0.4, color=1)
    plot!(gs, p4 , marker=true, ms=2, label="tetraploid", color=2,
        ribbon=(p4 .- first.(e4), last.(e4) .- p4), fillalpha=0.4)
    plot!(title="(A) \$\\theta=2.5, \\beta_4 = \\sqrt{1/2}\$", 
        ylabel="\$P\$", xlabel="\$\\gamma\$", ylim=(0,0.0101))
end
plot!()

P1 = plot()
map(zip([res0, res1], [2,2.5])) do (res, θ)
    nr = 1_000_000
    i  = 1
    p2 = getindex.(res, 2)[i:end]
    p4 = getindex.(res, 3)[i:end]
    x2 = floor.(Int, p2 .* nr)
    x4 = floor.(Int, p4 .* nr)
    e2 = jeffreys_interval.(x2, nr, 0.01)
    e4 = jeffreys_interval.(x4, nr, 0.01)
    plot!(gs[i:end], p2, marker=true, ms=2, label=res==res0 ? "diploid" : "", 
        ribbon=(p2 .- first.(e2), last.(e2) .- p2), fillalpha=0.4, color=1)
    plot!(gs[i:end], p4, marker=true, ms=2, label=res==res0 ? "tetraploid" : "", 
        ribbon=(p4 .- first.(e4), last.(e4) .- p4), fillalpha=0.4, color=2)
    plot!(title="(A) \$2V=1, \\beta_4 = \\sqrt{1/2}\$", yscale=:log10, 
        legend=:bottomright)
    x = 2.2
    y = (last(p4) + last(p2))/2
    annotate!(x, y, text("\$\\theta=$θ\$", :left, 9))
end
plot!(xlabel="\$\\gamma\$", ylabel="\$P\$")

getindex.(res1,3) ./ getindex.(res1,2)

# 2. Examine the effect of the variance scaling (β)
# -------------------------------------------------
bs = 0.25:0.05:1.0
resb = pmap(bs) do b
    nr = 1_000_000
    M = InfDemeMixEst(
        U=Umat(0.0,0.0),  # no unreduced gametes 
        γ=0.25,           # selection strength
        m=0.0,            # no migration
        θ=fill(2.5, 3),   # maladaptation
        β=[1.0, 1.0, √b], 
        V=0.5)            # segregation variance
    x4 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[4])), 1:nr)
    b, sum(first.(x4)) / nr
end

serialize("data/est1b-t2.5.jls", resb)

resb = deserialize("data/est1b.jls")

P2 = scatter(bs , last.(resb) ./ pp, 
    legend=false, color=:black, ms=2.5,
		ylabel="\$P_4/P_2\$", xlabel="\$\\beta_4^2\$", 
        #yscale=:log10, left_margin=5Plots.mm,
        xticks=0.25:0.25:1,
        title="(B) \$\\theta=2.5, \\gamma=0.25\$")
hline!([1], color=:firebrick, alpha=0.4)
vline!([0.5], color=:gray, ls=:dash)

plot(P1, P2, size=(260*2,210), margin=3Plots.mm)
#savefig("doc/img/est1b.pdf")


# 3. Examine the γ by β space
# ---------------------------
gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(2), length=15)]
bs = 0.25:0.05:1.0
gb = Iterators.product(gs, bs)

res = pmap(gb) do (g, b)
    @info (g, b)
    nr = 500_000
    M = InfDemeMixEst(
        U=Umat(0.0,0.0),  # no unreduced gametes 
        γ=g,           # selection strength
        m=0.0,            # no migration
        θ=[2.5,2.5,2.5],  # maladaptation
        β=[1.0, 1.0, √b],
        V=0.5)            # segregation variance
    x4 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[4])), 1:nr)
    @info (g,b,"tetraploid done", sum(first.(x4))/nr)
    (g, b, sum(first.(x4)) / nr)
end

res2 = pmap(gs) do g
    nr = 2_500_000
    M = InfDemeMixEst(
        U=Umat(0.0,0.0),  # no unreduced gametes 
        γ=g,           # selection strength
        m=0.0,            # no migration
        θ=[2.5,2.5,2.5],  # maladaptation
        V=0.5)            # segregation variance
    x2 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[2])), 1:nr)
    @info (g,"diploid done", sum(first.(x2))/nr)
    sum(first.(x2)) / nr
end

serialize("data/est1-2d-30032025.jls", (res, res2))

res, res2 = deserialize("data/est1-2d-30032025.jls")


rp = log10.(last.(res) ./ res2)

cc = maximum(abs, filter(isfinite, rp))
P3 = heatmap(bs, log10.(gs[2:end]), rp[2:end,:], color=:balance, clims=(-cc, cc), 
#    xticks=(1:2:length(bs), bs[1:2:end]),
#    yticks=(1:length(gs)-1, map(x->@sprintf("%.2f", x), gs[2:end]))
    ylabel="\$\\log_{10}\\gamma\$",
    xlabel="\$\\beta_4\$",
    margin=4Plots.mm,
    title="(C) \$\\log_{10}\\left(P_4/P_2\\right)\$",
    size=(320,280)
)
vline!([√(1/2)], color=:black, lw=2, legend=false)



# combined plot
P2 = plot!(P2, left_margin=14Plots.mm)
plot(P1, P2, P3, size=(800,200), margin=5Plots.mm,
    layout=grid(1,3,widths=[0.3,0.30,0.40]))

savefig("doc/img/fig2.pdf")


# 4.
# --

gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(2), length=15)]
θs = range(1.0, stop=3.0, length=11)
gt = Iterators.product(gs, θs)

res = pmap(gt) do (g, t)
    @info (g, t)
    nr = 500_000
    M = InfDemeMixEst(
        U=Umat(0.0,0.0),  # no unreduced gametes 
        γ=g,           # selection strength
        m=0.0,         # no migration
        θ=[t,t,t],     # maladaptation
        V=0.5)         # segregation variance
    x4 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[4])), 1:nr)
    @info (g,t,"tetraploid done", sum(first.(x4))/nr)
    x2 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[2])), 1:nr)
    @info (g,t,"diploid done", sum(first.(x2))/nr)
    (g, t, sum(first.(x2)) / nr, sum(first.(x4)) / nr)
end

rp = log10.(last.(res) ./ getindex.(res,3))

cc = maximum(abs, filter(isfinite, rp))
P4 = heatmap(θs, log10.(gs[2:end]), rp[2:end,:], color=:balance, clims=(-cc, cc), 
#    xticks=(1:2:length(bs), bs[1:2:end]),
#    yticks=(1:length(gs)-1, map(x->@sprintf("%.2f", x), gs[2:end]))
    ylabel="\$\\log_{10}\\gamma\$",
    xlabel="\$\\beta_4\$",
    margin=4Plots.mm,
    title="(C) \$\\log_{10}\\left(P_4/P_2\\right)\$",
    size=(320,280)
)

#serialize("data/est1-4-02042025.jls", res)
res = deserialize("data/est1-4-02042025.jls")

Ew = map(gt) do (g, t)
    exp(-g*t + g^2/2)
end 
plot(gs, Ew, legend=:outertopright, label=reshape(θs, 1,11), ylabel="\$E[w]\$", xlabel="\$\\gamma\$")

P2 = getindex.(res,3)
P4 = getindex.(res,4)
X = (P2 ./ Ew)
heatmap(θs, gs, Ew, )

plot()
for i=2:length(gs), j=1:length(θs)
    scatter!([gs[i]^2/θs[j]], [P2[i,j]], color=1)
    scatter!([gs[i]^2/θs[j]], [P4[i,j]], color=2)
end
plot!(legend=false)

# ---
gs = 0:0.05:1.5  # g/θ
θs = [2.0, 2.5]
gt = Iterators.product(gs, θs)

res = pmap(gt) do (gx, t)
    g = gx*t
    @info (g, t)
    nr = 500_000
    M = InfDemeMixEst(
        U=Umat(0.0,0.0),  # no unreduced gametes 
        γ=g,           # selection strength
        m=0.0,         # no migration
        θ=[t,t,t],     # maladaptation
        V=0.5)         # segregation variance
    x4 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[4])), 1:nr)
    @info (g,t,"tetraploid done", sum(first.(x4))/nr)
    x2 = map(_->simest2(M, x0=InfPop(z=[0.0], c=[2])), 1:nr)
    @info (g,t,"diploid done", sum(first.(x2))/nr)
    (t, gx, sum(first.(x2)) / nr, sum(first.(x4)) / nr)
end

plot(gs, last.(res))
plot!(gs, getindex.(res, 3))

