using Distributed, StatsBase, Plots, Serialization, Distributions, Optim
addprocs(6)
@everywhere using InfGenetics

# estimate parameter of a censored geometric distribution by ML
function logpdf_censored_geometric(p, c, x)
    x < c ? log(p) + (x-1)*log(1-p) : (c-1)*log(1-p)
end

function mle_censored_geometric(c, xs)
    fun = p -> -sum(logpdf_censored_geometric.(Ref(p), Ref(c), xs))
    optimize(fun, 0.0, 1.0).minimizer
end


# Basic model
# ===========
nrep = 1000
ms   = 10 .^ range(log10(0.1), stop=log10(3), length=10)
zs   = [2.0, 1.5, 1.0]
γ    = 0.25
U    = Umat(0.05, 0.05)
nmax = 10000
Nest = 100

sims = pmap(Iterators.product(ms, zs)) do (m, z)
    @info (m, z)
    reps = map(1:nrep) do _
        sim = simest2(InfDemeMixEst(m=m, θ=fill(z, 3), U=U), 
            nmax=nmax, Nest=Nest)
    end
    (m, z, reps)
end

sims = deserialize("data/sims1.jls")

P1 = plot(title="(A) Time to establishment")
for col in eachcol(sims)
    t = map(col) do reps
        ts = getindex.(reps[3], 2)
        1/mle_censored_geometric(nmax, ts)
    end
    z = col[1][2]
    plot!(ms, t, marker=true, ms=2, label="\$\\bar{z}_s = $z\$")
end
P1 = plot(P1, marker=true, ms=2, yscale=:log10, xlabel="\$m\$",
    ylabel="\$\\overline{T}\$", xscale=:log10, legend=:topleft)

P2 = plot(title="(B) Tetraploid establishment", xscale=:log10)
map(eachcol(sims)) do col
    xx = map(col) do reps
        nsucc = sum(first.(reps[3]))
        sum(map(x->x[3][3] > x[3][1] && x[1], reps[3])), nsucc
    end
    x = first.(xx)
    n = last.(xx)
    @info n
    e = jeffreys_interval.(x, n, 0.05)
    p = x ./ n
    z = col[1][2]
    plot!(ms, p, 
        ribbon=(p .- first.(e), last.(e) .- p), fillalpha=0.2,
        marker=true, ms=2, label="\$\\bar{z}_s = $z\$")
end
hline!([InfGenetics.cytotype_equilibrium(U)[3]], label="\$\\pi_4\$",
    color=:black, xlabel="\$m\$", ylabel="\$P_4\$", ls=:dash)

plot(P1, P2, size=(500,200), margin=3Plots.mm, titlefont=8)
#savefig("doc/img/estwmigration1.pdf") |> texfig


# With selfing
# =============
# This is not just with selfing, but with loss of SI...
nrep = 5000
ms   = 10 .^ range(log10(0.1), stop=log10(3), length=10)
zs   = -1.5
σs   = [[NaN, s, s] for s=0:0.2:1.0]
γ    = 0.25
u    = 0.05
U    = Umat(u, u)
nmax = 20000
Nest = 100

sims = pmap(Iterators.product(ms, σs)) do (m, σ)
    @info (m, σ)
    reps = map(1:nrep) do _
        sim = simest2(InfDemeMixEst(m=m, θ=fill(zs, 3), U=U, σ=σ), 
            nmax=nmax, Nest=Nest)
    end
    (m, σ[3], reps)
end

sims = deserialize("data/sims-si1.jls")

P1 = plot(title="(A) Time to establishment")
for col in eachcol(sims)
    t = map(col) do reps
        1/mle_censored_geometric(nmax, getindex.(reps[3], 2))
    end
    s = col[1][2]
    plot!(ms, t, marker=true, ms=2, label="\$\\sigma = $s\$")
end
P1 = plot(P1, marker=true, ms=2, xlabel="\$m\$",
    ylabel="\$\\overline{T}\$", xscale=:log10, legend=false, yscale=:log10)

P2 = plot(title="(B) Tetraploid establishment", xscale=:log10)
map(eachcol(sims)) do col
    xx = map(col) do reps
        nsucc = sum(first.(reps[3]))
        sum(map(x->x[3][3] > x[3][1] && x[1], reps[3])), nsucc
    end
    x = first.(xx)
    n = last.(xx)
    @info n
    e = jeffreys_interval.(x, n, 0.05)
    p = x ./ n
    s = col[1][2]
    plot!(ms, p, 
        ribbon=(p .- first.(e), last.(e) .- p), fillalpha=0.2,
        marker=true, ms=2, label="\$\\sigma = $s\$")
end
hline!([InfGenetics.cytotype_equilibrium(U)[3]], label="\$\\pi_4\$",
    color=:black, xlabel="\$m\$", ylabel="\$P_4\$", ls=:dash, 
    legend=:outertopright, size=(350,220))

plot(P1, P2, size=(540,200), layout=grid(1,2,widths=[0.415,0.585]),
    margin=3Plots.mm, titlefont=8)
savefig("doc/img/estwmigration-si.pdf") |> texfig

P1 = plot()
for (j,col) in enumerate(eachcol(sims))
    t = map(col) do reps
        1/mle_censored_geometric(nmax, getindex.(reps[3], 2))
    end
    x = map(col) do reps
        nsucc = sum(first.(reps[3]))
        sum(map(x->x[3][3] > x[3][1] && x[1], reps[3])), nsucc
    end
    p = first.(x) ./ last.(x)
    s = col[1][2]
    scatter!(t[1:1], p[1:1], color=j, label="", ms=2)
    for i=1:length(col)-1
        plot!(t[i:i+1], p[i:i+1], color=j, marker=false, ms=2, 
            label=i==1 ? "\$\\sigma = $s\$" : "", arrow=(:closed, 2.0))
    end
end
P1 = plot(P1, marker=true, ms=2, ylabel="\$P_4\$", size=(400,300), grid=true,
    xlabel="\$\\overline{T}\$", cscale=:log10, legend=:outertopright, 
    title="\$\\bar{z}_s = $zs, \\gamma=$γ, u=v=$u\$")
#savefig("doc/img/estwmigration-si2.pdf") |> texfig

# With assortative mating
# =======================
nrep = 5000
ms   = 10 .^ range(log10(0.1), stop=log10(3), length=10)
zs   = -1.5
ρs   = [fill(r, 3) for r=0:0.2:1.0]
γ    = 0.25
U    = Umat(0.05, 0.05)
nmax = 20000
Nest = 100

sims = pmap(Iterators.product(ms, ρs)) do (m, ρ)
    @info (m, ρ)
    reps = map(1:nrep) do _
        sim = simest2(InfDemeMixEst(m=m, θ=fill(zs, 3), U=U, ρ=ρ), 
            nmax=nmax, Nest=Nest)
    end
    (m, ρ[3], reps)
end

sims = deserialize("data/sims-am1.jls")

P1 = plot(title="(A) Time to establishment")
for col in eachcol(sims)
    t = map(col) do reps
        1/mle_censored_geometric(nmax, getindex.(reps[3], 2))
    end
    s = col[1][2]
    plot!(ms, t, marker=true, ms=2, label="\$\\rho = $s\$")
end
P1 = plot(P1, marker=true, ms=2, xlabel="\$m\$",
    ylabel="\$\\overline{T}\$", xscale=:log10, legend=false, yscale=:log10)

P2 = plot(title="(B) Tetraploid establishment", xscale=:log10)
map(eachcol(sims)) do col
    xx = map(col) do reps
        nsucc = sum(first.(reps[3]))
        sum(map(x->x[3][3] > x[3][1] && x[1], reps[3])), nsucc
    end
    x = first.(xx)
    n = last.(xx)
    @info n
    e = jeffreys_interval.(x, n, 0.05)
    p = x ./ n
    s = col[1][2]
    plot!(ms, p, 
        ribbon=(p .- first.(e), last.(e) .- p), fillalpha=0.2,
        marker=true, ms=2, label="\$\\rho = $s\$")
end
hline!([InfGenetics.cytotype_equilibrium(U)[3]], label="\$\\pi_4\$",
    color=:black, xlabel="\$m\$", ylabel="\$P_4\$", ls=:dash, 
    legend=:outertopright, size=(350,220))

plot(P1, P2, size=(540,200), layout=grid(1,2,widths=[0.415,0.585]),
    margin=3Plots.mm, titlefont=8)

P1 = plot()
for col in eachcol(sims)
    t = map(col) do reps
        1/mle_censored_geometric(nmax, getindex.(reps[3], 2))
    end
    x = map(col) do reps
        nsucc = sum(first.(reps[3]))
        sum(map(x->x[3][3] > x[3][1] && x[1], reps[3])), nsucc
    end
    p = first.(x) ./ last.(x)
    s = col[1][2]
    plot!(ms, t, p, marker=true, ms=2, label="\$\\rho = $s\$")
end
P1 = plot(P1, marker=true, ms=2, zlabel="\$P_4\$", size=(600,600), grid=true,
    ylabel="\$\\overline{T}\$", yscale=:log10, legend=false)


# larger panel of conditions
# ==========================
nrep = 5000
nmax = 10000
Nest = 100
zs   = -1.5
γ    = 0.25
ms   = 10 .^ range(log10(0.1), stop=log10(3), length=10)
ρs   = [fill(r, 3) for r=0:0.25:1.0]
αs   = [zeros(3), [1/2, 1/4, 1/6]]
βs   = [ones(3), [1, √(2/3), √(1/2)], [1, 3/4, 1/2]]
us   = [0.025, 0.05]

params = Iterators.product(ms, ρs, αs, βs, us)
@info length(params)

sims = pmap(params) do (m, ρ, α, β, u)
    @info (m, ρ, α[end], β[end], u)
    reps = map(1:nrep) do _
        M   = InfDemeMixEst(m=m, θ=fill(zs, 3), U=Umat(u,u), α=α, β=β, ρ=ρ)
        sim = simest2(M, nmax=nmax, Nest=Nest)
    end
    (M, reps)
end

sims = deserialize("data/sims-am2.jls")

P1 = plot(title="(A) Time to establishment")
for (r, col) in zip(ρs, eachcol(sims[:,:,1,1,1]))
    t = map(col) do (M, reps)
        1/mle_censored_geometric(nmax, getindex.(reps, 2))
    end
    plot!(ms, t, marker=true, ms=2, label="\$\\rho = $r\$")
end
P1 = plot(P1, marker=true, ms=2, xlabel="\$m\$",
    ylabel="\$\\overline{T}\$", xscale=:log10, 
    legend=false, yscale=:log10)

map(enumerate(αs)) do (i, α)
    P1 = plot(title="\$\\alpha = $α\$", xscale=:log10, yscale=:log10)
    P2 = plot(xscale=:log10)
    P3 = plot(xscale=:log10)
    for (j, (r, col)) in enumerate(zip(ρs, eachcol(sims[:,:,i,1,2])))
        t = map(col) do (M, reps)
            1/mle_censored_geometric(nmax, getindex.(reps, 2))
        end
        plot!(P1, ms, t, marker=true, ms=2, label="\$\\rho = $r\$", legend=false)
        xx = map(col) do (M, reps)
            nsucc = sum(first.(reps))
            sum(map(x->x[3][3] > x[3][1] && x[1], reps)), nsucc
        end
        x = first.(xx)
        n = last.(xx)
        e = jeffreys_interval.(x, n, 0.05)
        p = x ./ n
        @info n, p
        plot!(P2, ms, p, 
            label="\$\\rho=$(@sprintf "%.2f" r[1])\$",
            ribbon=(p .- first.(e), last.(e) .- p), fillalpha=0.2,
            marker=true, ms=2)
        plot!(P3, t, p, marker=true, ms=2, legend=false)
    end
    hline!(P2, [0.000731], color=:black, ls=:dash, label="",
        legend=:topleft, xlabel="\$m\$", ylabel="\$P_4\$")
    plot(P1, P2, P3, layout=(1,3))
end |> x->plot(x..., layout=(length(x),1), size=(750,length(x)*200))

plot(P1, P2, size=(550,200), margin=3Plots.mm)
