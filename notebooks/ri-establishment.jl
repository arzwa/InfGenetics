# Consider for the case of succesful establishment (without/with migration),
# how relative fitness of migrants evolves
using Plots
using Parameters
using Distributed
using Serialization
addprocs(5)
@everywhere using InfGenetics
@everywhere using StatsBase

@everywhere const Sims = Vector{InfPop{Float64}}

# collect n2 and n4 establishment replicates for both diploid and tetraploids
@everywhere function collect_est(f, n2, n4)
    x2 = Vector{Sims}(undef, n2)
    x4 = Vector{Sims}(undef, n4)
    c2 = c4 = 1
    while c2 <= n2 || c4 <= n4
        res = simuntil(f)
        cs = counts(res[end][2].c, 2:4)
        if cs[3] > cs[1] && c4 <= n4
            x4[c4] = last.(res)
            c4 += 1
            @info c4, c2
        elseif cs[1] > cs[3] && c2 <= n2
            x2[c2] = last.(res)
            c2 += 1
            @info c4, c2
        end
    end
    return x2, x4
end

function avg(xs)
    map(i->mean(filter(
        y->y>0, [x[i] for x in xs if length(x) >= i])), 
        1:maximum(length, xs))
end

# initial migrant z=0
M = InfDemeMixEst(U=Umat(0.00,0.00), γ=0.25, m=0.0, θ=[2.0,2.0,2.0], V=0.5)
f2=()->simest(reconstruct(M, γ=0.25), x0=InfPop(z=zeros(1), c=Int64[2]))
f4=()->simest(reconstruct(M, γ=0.25), x0=InfPop(z=zeros(1), c=Int64[4]))
x2, _ = collect_est(f2, 1000, 0) 
_, x4 = collect_est(f4, 0, 1000) 

# initial migrant random
M = InfDemeMixEst(U=Umat(0.00,0.00), γ=0.25, m=0.0, θ=[2.0,2.0,2.0], V=0.5)
f2=()->simest(reconstruct(M, γ=0.25), x0=InfPop(z=randn(1), c=Int64[2]))
f4=()->simest(reconstruct(M, γ=0.25), x0=InfPop(z=randn(1), c=Int64[4]))
x2b, _ = collect_est(f2, 1000, 0) 
_, x4b = collect_est(f4, 0, 1000) 

# nonzero migration rates
ms = 10 .^ collect(range(-2, 0, 5))
Xs = pmap(ms) do m
    M = InfDemeMixEst(U=Umat(0.05,0.05), γ=0.25, m=m, θ=[2.0,2.0,2.0], V=0.5)
    f = ()->simest(M, x0=InfPop(z=zeros(0), c=Int64[]))
    y2, y4 = collect_est(f, 1000, 1000) 
    @info m
    m, y2, y4
end
    
M = InfDemeMixEst(U=Umat(0.05,0.05), γ=0.25, m=2.0, θ=[2.0,2.0,2.0], V=0.5)
f = ()->simest(M, x0=InfPop(z=zeros(0), c=Int64[]))
y2, y4 = collect_est(f, 100, 100) 

Ys = [(0.0, x2b, x4b) ; Xs]
## serialize("data/1000reps.jls", Ys)
#Ys = deserialize("data/1000reps.jls")

function bypopsize(x::Vector{T}, Nmax=100) where T
    d = Dict{Int64,T}()
    for rep in x
        for gen in rep
            N = length(gen.z)
            (N >= Nmax || N <= 0) && continue
            if haskey(d, N) 
                push!(d[N], gen)
            else
                d[N] = [gen]
            end
        end
    end
    return d
end

function bysubpopsize(x::Vector{T}; Nmax=100, k=2) where T
    d = Dict{Int64,T}()
    for rep in x
        for gen in rep
            N = counts(gen.c, 1:4)[k]
            (N >= Nmax || N <= 0) && continue
            if haskey(d, N) 
                push!(d[N], gen)
            else
                d[N] = [gen]
            end
        end
    end
    return d
end

Zs = map(Ys) do (m, y2, y4)
    m, bypopsize(y2), bypopsize(y4)
end
Zs = Zs[[1,2,4,6]]

Ps = map(Zs) do (m, X2, X4)
    P1 = scatter([(k,mean(mapreduce(x->x.F, vcat, v))) for (k,v) in X2], 
        xscale=:log10, ms=1.5, label="diploid", title="\$m=$(@sprintf "%.2f" m)\$")
    scatter!([(k,mean(mapreduce(x->x.F, vcat, v))) for (k,v) in X4], 
        xscale=:log10, ms=1.5, label="tetraploid")
    plot!(legend=m == 0.0 ? :bottomright : false, 
        xlabel="", ylabel="\$\\bar{F}\$", ylim=(0,1))
    P1b = scatter([(k,mean(mapreduce(x->x.z, vcat, v))) for (k,v) in X2], 
        xscale=:log10, ms=1.5, label="diploid", legend=false)
    scatter!([(k,mean(mapreduce(x->x.z, vcat, v))) for (k,v) in X4], 
            xscale=:log10, ms=1.5, label="tetraploid", ylim=(1,6))
    plot!(xlabel="\$N\$", ylabel="\$\\bar{z}\$")
    hline!([2], color=:black, alpha=0.2)
    plot(P1, P1b, size=(600,240),layout=(2,1))
end
plot(Ps..., layout=(1,length(Zs)), size=(700,280), 
    bottom_margin=5Plots.mm, left_margin=4Plots.mm)

savefig("doc/img/fig4.pdf")

z2b = bypopsize(x2b)
z4b = bypopsize(x4b)
scatter!(Ps[1][1], [(k,mean(mapreduce(x->x.F, vcat, v))) for (k,v) in z2b], ms=1.5, label="", 
    markerstrokewidth=1, markerstrokecolor=1, color=:white,  alpha=0.6)
    #color=1, alpha=0.3)
scatter!(Ps[1][1], [(k,mean(mapreduce(x->x.F, vcat, v))) for (k,v) in z4b], ms=1.5, label="", 
    markerstrokewidth=1, markerstrokecolor=2, color=:white,  alpha=0.6)
    #color=2, alpha=0.3)
scatter!(Ps[1][2], [(k,mean(mapreduce(x->x.z, vcat, v))) for (k,v) in z2b], ms=1.5, label="", 
    markerstrokewidth=1, markerstrokecolor=1, color=:white,  alpha=0.6)
    #color=1, alpha=0.3)
scatter!(Ps[1][2], [(k,mean(mapreduce(x->x.z, vcat, v))) for (k,v) in z4b], ms=1.5, label="",
    markerstrokewidth=1, markerstrokecolor=2, color=:white,  alpha=0.6)
    #color=2, alpha=0.3)
plot(Ps..., layout=(1,length(Zs)), size=(700,280), 
    bottom_margin=5Plots.mm, left_margin=4Plots.mm)

savefig("doc/img/fig4.pdf")

# This is hard to interpret? Why do tetraploids start of slower and lower?
# Well, it may be that, conditional on establishment, tetraploids need not
# change the mean so much in the first gens to eventually establish, this would
# lead to a 'bias' like the one observed... Tetraploids need not adapt as fast
# in order to establish eventually, hence conditional on establishment one
# observes slower change in trait mean.
#
C = map(Zs) do (m, z2, z4)
    map(y->mean(mapreduce(x->proportions(x.c, 2:4)[end], vcat, y)), 
        collect(values(sort(z4))))
end

plot(C[2:end], xscale=:log10)

xs = map(C[2:end]) do c
    findfirst(x->x>0.75, c)
end
vline!(xs)



# m = 0 case ---------------------------------------------------------

Zs = map(Ys) do (m, y2, y4)
    m, bypopsize(y2), bypopsize(y4)
end

map([1,5,10,20,40,80]) do k
    z2est = mapreduce(x->x.z, vcat, Zs[1][2][k])
    z4est = mapreduce(x->x.z, vcat, Zs[1][3][k])
    bins = -1.5:0.2:5
    stephist(vcat(z2est...), bins=bins, label="diploid", title="\$N=$k\$", normalize=true)
    stephist!(vcat(z4est...), bins=bins, label="tetraploid", normalize=true)
    vline!([0,2], label="", legend=k==20 ? :topleft : false)
end |> x->plot(x..., size=(700,300))


map(Zs) do Z
    v2 = map(sort(collect(keys(Z[2])))) do k
        z2est = map(x->var(x.z), Z[2][k])
        k, mean(z2est)
    end
    v4 = map(sort(collect(keys(Z[3])))) do k
        z4est = map(x->var(x.z), Z[3][k])
    k, mean(z4est)
    end
    plot(v2, xscale=:log10)
    plot!(v4)
end |> x->plot(x..., size=(700,300), ylim=(0,1.3))

