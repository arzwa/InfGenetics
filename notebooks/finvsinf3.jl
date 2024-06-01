using Distributed, Distributions, StatsBase, Serialization
addprocs(6)
@everywhere using InfGenetics, Distributions

nrep = 10
ngen = 2000
L    = 500
N    = 250
M    = InfDemeMix(U=[1.0 0.0; 0.0 0.0; 0.0 1.0])

@everywhere function cb(x)  # callback function
    _, pop = x
    z2 = pop.z[pop.c .== 2]
    z3 = pop.z[pop.c .== 3]
    z4 = pop.z[pop.c .== 4]
    z2, z3, z4
end

# Tetraploids
Xi = pmap(1:nrep) do _
    P = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(4,N))
    X = InfGenetics.fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end
Xf = pmap(1:nrep) do _
    G = [randn(L,4)/√(2L) for i=1:N]
    P = FinPop(sum.(G), G, fill(4, N))
    X = fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end
vsi = mapreduce(x->map(var, last.(x)), hcat, Xi) 
vsf = mapreduce(x->map(var, last.(x)), hcat, Xf) 

P1 = plot(mean(vsi, dims=2), label="\$L \\rightarrow \\infty\$")
plot!(mean(vsf, dims=2)    , label="\$L = $L\$")
plot!(t->exp(-t/4N), 
    xlabel="generation", 
    ylabel="\$V\$", 
    color=:black, ls=:dash, label="") 

# Tetraploids with double reduction
M   = InfDemeMix(U=[1.0 0.0; 0.0 0.0; 0.0 1.0], α=[1/6, 1/6, 1/6])
Xib = pmap(1:nrep) do _
    P = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(4,N))
    X = InfGenetics.fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end    
Xfb = pmap(1:nrep) do _
    G = [randn(L,4)/√(2L) for i=1:N]
    P = FinPop(sum.(G), G, fill(4, N))
    X = fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end
vsib = mapreduce(x->map(var, last.(x)), hcat, Xib) 
vsfb = mapreduce(x->map(var, last.(x)), hcat, Xfb) 

P2 = plot(mean(vsib, dims=2), label="\$L \\rightarrow \\infty\$")
plot!(mean(vsfb, dims=2)    , label="\$L = $L\$")
plot!(t->exp(-t/4N), 
    xlabel="generation", 
    ylabel="\$V\$", 
    color=:black, ls=:dash, label="") 

# Triploids
M   = InfDemeMix(U=[0.0 0.0; 0.05 0.05; 0.0 0.0], w=[0.0,1.0,0.0], α=[1/6, 1/4, 1/6])
Xic = pmap(1:nrep) do _
    P = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(3,N))
    X = InfGenetics.fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end    
Xfc = pmap(1:nrep) do _
    G = [randn(L,3)/√(2L) for i=1:N]
    P = FinPop(sum.(G), G, fill(3, N))
    X = fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end
vsic = mapreduce(x->map(var, last.(x)), hcat, Xic) 
vsfc = mapreduce(x->map(var, last.(x)), hcat, Xfc) 

P3 = plot(mean(vsic, dims=2), label="\$L \\rightarrow \\infty\$")
plot!(mean(vsfc, dims=2)    , label="\$L = $L\$")
plot!(t->exp(-t/3N), 
    xlabel="generation", 
    ylabel="\$V\$", 
    color=:black, ls=:dash, label="") 

plot!(P1, title="tetraploids")
plot!(P2, title="tetraploids, \$\\alpha=1/6\$")
plot!(P3, title="triploids, \$\\alpha=1/4\$")
plot(P1, P2, P3, layout=(1,3), size=(600,160), margin=3Plots.mm, guidefont=9, titlefont=9)


nrep = 1000
M   = InfDemeMix(U=Umat(0.05,0.05), α=[1/4, 1/4, 1/6])
Xid = pmap(1:nrep) do _
    P = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
    X = InfGenetics.fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end    
Xfd = pmap(1:nrep) do _
    G = [randn(L,3)/√(2L) for i=1:N]
    P = FinPop(sum.(G), G, fill(2, N))
    X = fiterate((x,i)->generation(x, M, z->1), deepcopy(P), ngen, callback=cb)
end

res = map(1:3) do i
    map([Xid, Xfd]) do X 
        xs = map(x->getindex.(x, i), X)
        Z = mapreduce(i->map(x->isempty(x) ? NaN : var(x), xs[i]), hcat, 1:nrep)
        map(x->filter(!isnan, x) |> y->isempty(y) ? NaN : mean(y), eachrow(Z))
    end
end
