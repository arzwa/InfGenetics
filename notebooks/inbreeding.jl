
nrep = 100
m    = 1.0
zs   = 2.0
#ρs   = [fill(r, 3) for r=0:0.2:1.0]
σs   = [[0.0, s, s] for s=0:0.2:1.0]
γ    = 0.25
U    = Umat(0.05, 0.05)
nmax = 10000
Nest = 100

function tetest(x) 
    y = counts(x.c, 2:4)
    y[3] > y[1]
end

sims = pmap(σs) do σ
    @info σ
    reps2 = Vector{Float64}[]
    reps4 = Vector{Float64}[]
    while length(reps2) < nrep || length(reps4) < nrep 
        _, sim = simest(
            InfDemeMixEst(m=m, θ=fill(zs, 3), U=U, σ=σ), nmax=nmax, Nest=Nest)
        ff = map(x->mean(x[2].F), sim)
        if tetest(sim[end][2]) 
            length(reps4) < nrep && push!(reps4, ff) 
        else
            length(reps2) < nrep && push!(reps2, ff)
        end
    end
    reps2, reps4
end

map(sims) do (s2, s4)
    x2 = avg(s2)
    x4 = avg(s4)
    plot(x2, xscale=:log10, ylim=(0,1))
    plot!(x4)
end |> x->plot(x..., legend=false, size=(900,300))
