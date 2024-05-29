using InfGenetics, Distributions, StatsBase

N = 500
M = InfDemeMix(α=zeros(3), β=ones(3))
pop = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))

@btime generation(pop, M, z->1.0)

# In an infinite population (i.e. no inbreeding) we should get to m*M.V
# quickly, where m is the ploidy level, at least in the absence of double
# reduction.
for i=1:100
    pop = generation(M, z->1, pop)
    pop.F .= 0.0
    pop.Φ .= 0.0
    @info var(pop.z), 2M.V
end

function gen(pop, M)
    pop = generation(pop, M)
    pop.z .-= mean(pop.z)
    pop.F .= 0.0
    pop.Φ .= 0.0
    return pop
end
    
function cb(x)
    _, pop = x
    z2 = pop.z[pop.c .== 2]
    z3 = pop.z[pop.c .== 3]
    z4 = pop.z[pop.c .== 4]
    z2, z3, z4
end


bs = [ones(3), [1.0, 0.75, 0.5], [1.0, √(2/3), √(1/2)]]
ps = map(bs) do b
    M = InfDemeMix(α=zeros(3), β=b)
    N = 500
    zinit = rand(Normal(0, √(2M.V)), N)
    pop = InfPop(z=zinit, c=fill(2,N))
    res = fiterate((pop,i)->gen(pop, M), pop, 1000, callback=cb)
    zs = map(i->vcat(getindex.(res, i)...), 1:3) 
    plot(title="\$\\beta = $(round.(b, digits=2))\$")
    map(i->stephist!(zs[i], norm=true, bins=-4:0.3:4,
        label=["diploid", "triploid", "tetraploid"][i]), 1:3)
    map(zip(b,2:4)) do (bb, k)
        plot!(Normal(0, bb*√(k*M.V)), color=k-1, 
            label="",
            fill=true, alpha=0.2, linealpha=0.2)
    end
    plot!()
end

plot(ps..., ylabel="density", legend=:topright, layout=(1,3), size=(700,200), xlabel="\$z\$", margin=3Plots.mm)

M = InfDemeMix()
N = 500
zinit = rand(Normal(0, √(2M.V)), N)
pop = InfPop(z=zinit, c=fill(2,N))
res = fiterate((pop,i)->gen(pop, M), pop, 1000, callback=x->proportions(x[2].c, 2:4))

plot(permutedims(hcat(res...)), label="", alpha=0.3, legend=:right,
    size=(300,240), xlabel="generation", ylabel="proportion", title="\$u=v=0.05\$")
g1 = InfGenetics.mpeq(0.05, 0.05)
hline!([g1^2], color=1, label="diploids")
hline!([2g1*(1-g1)], color=2, label="triploids")
hline!([(1-g1)^2], color=3, label="tetraploids")

u = M.U[1,2]; v = M.U[2,1]
hline!([1-2u-u^2-4u*v], color=1, label="diploids")
hline!([4u*v + 2u], color=2, label="triploids")
hline!([u^2], color=3, label="tetraploids")


# ------------
for i=1:100
    pop = generation(pop, M, z->exp(-2z^2))
    pop.F .= 0.0
    pop.Φ .= 0.0
    @info var(pop.z), 2M.V
end

function simshift(M, N, γ, θ, n1=100, n2=100)
    pop = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
    # equilibrate, assume the population is very large
    ys = [(z=pop.z, c=pop.c)]
    for _=1:n1
        pop = generation(pop, M, z->exp(-γ*z^2))
        pop.F .= 0.0
        pop.Φ .= 0.0
        push!(ys, (z=pop.z, c=pop.c))
    end
    for _=1:n2
        pop = generation(pop, M, z->exp(-γ*(z - θ)^2))
        push!(ys, (z=pop.z, c=pop.c))
    end
    return ys
end

function tetraploid_est(M, N, γ, θ, n1=100, n2=100)
    k = 0
    while true
        ys = simshift(M, N, γ, θ, n1, n2)  
        d = counts(ys[end].c, 1:4)
        d[4] > d[2] && return (k+1, ys)
        k += 1
    end
end

N = 200
γ = 0.5
θ = 5
M = InfDemeMix(β=ones(3))
xs1 = simshift(M, N, γ, θ, 100, 100)
tw, xs2 = tetraploid_est(M, N, γ, θ, 100, 100)

[exp(-γ*θ + γ^2/2*m*M.β[m-1]^2*M.V) for m=2:4]


xs0 = simshift(M, N, γ, θ, 100, 0)

mean(map(z->exp.(γ*(z-θ)), xs0[end].z))


plot( map(x->mean(x.z), xs1))
plot!(map(x->mean(x.z), xs2))

plot( map(x->var( x.z), xs1))
plot!(map(x->var( x.z), xs2))

map(1:20) do _
    xs = tetraploid_est(M, N, γ, θ)
    @info xs[1]; xs
end

bs = 0.60:0.05:1.3
res = map(bs) do β4
    βs = [1.0, (1+β4)/2, β4]
    @info βs
    M = InfDemeMix(α=fill(0.2, 3), β=βs)
    xs = [tetraploid_est(M, N, γ, θ) for i=1:40]
    ys = first.(xs)
    @info 1/mean(ys); ys
end

dd = map(x->Beta(length(x), sum(x) - length(x)), res)
q1 = quantile.(dd, 0.025)
q2 = quantile.(dd, 0.975)
ms = mean.(dd)
plot(bs, ms, ribbon=(ms .- q1, q2 .- ms), size=(300,250), legend=false, 
    title="\$\\gamma=$γ, \\theta=$θ, N=$N\$\n\$\\alpha=\\xi=0.2, u=$(M.U[1,2]), v=$(M.U[2,1])\$",
    xlabel="\$\\beta_4\$", ylabel="\$P\$", ylim=(0,0.9))


# establishment
M = InfGenetics.InfDemeMixEst(α=zeros(3), β=ones(3))
N = 0
t = 0
pop = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
while true
    t += 1
    pop = generation(pop, M)
    length(pop.z) > 50 && break
end
(t, pop.z, pop.c)

function simest(M, minn=50)
    pop = InfPop(z=zeros(0))
    t = 0
    while true
        t += 1
        pop = generation(pop, M)
        length(pop.z) > minn && (return t, pop)
    end
end

function simtetest(M, minn=50)
    T = 0
    while true
        T += 1
        t, pop = simest(M, minn)
        d = counts(pop.c, 1:4)
        d[4] > d[2] && (return T, pop)
    end
end

function simdipest(M, minn=50)
    T = 0
    while true
        T += 1
        t, pop = simest(M, minn)
        d = counts(pop.c, 1:4)
        d[2] > d[4] && (return T, pop)
    end
end


# migration focus
ms = 10 .^ collect(-1:0.1:0.4)
res = map(ms) do m
    M = InfGenetics.InfDemeMixEst(α=fill(0.2, 3), β=ones(3), m=m)
    y = [simest(M) for i=1:500]
    ts = mean(first.(y))
    cs = map(x->counts(x.c, 1:4), last.(y))
    @info m, mean(ts)
    ts, cs
end

ptet = map(x->mean(map(y->y[4] > y[2], x)), last.(res))
tws  = first.(res)


bs = [ones(3), [1.0, 0.75, 0.5], [1.0, √(2/3), √(1/2)]]

bs = [ones(3)]
res = map(bs) do b
    map(10 .^ collect(-2:0.5:0.5)) do m
        @info b, m
        M = InfGenetics.InfDemeMixEst(α=fill(0.2, 3), β=b, m=m)
        t = [simest(M) for i=1:100]
        mean(first.(t))
    end
end

xs = res[1]
plot(-2:0.5:0.5, xs, yscale=:log10, marker=true, ms=4)

plot(1 ./ first.(xs), yscale=:log10)
plot!(1 ./ last.(xs))
