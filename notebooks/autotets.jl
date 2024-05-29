
M = InfDemeMix(α=fill(0.0, 3), U=[0. 0.; 0. 0.; 0. 1.0])
N = 200
zinit = rand(Normal(0, √(4M.V)), N)
pop = InfPop(z=zinit, c=fill(4,N))
xs = fiterate((pop,i)->generation(pop, M), pop, 1000, callback=x->(z=x[2].z,F=x[2].F))

plot(map(x->var(x.z), xs))
plot!(t->exp(-t/4N))

α=0.3
M = InfDemeMix(α=fill(α, 3), U=[0. 0.; 0. 0.; 0. 1.0])
N = 200
zinit = rand(Normal(0, √(4M.V)), N)
pop = InfPop(z=zinit, c=fill(4,N))
xs = fiterate((pop,i)->generation(pop, M), pop, 500, callback=x->(z=x[2].z,F=x[2].F))

F = map(x->mean(x.F), xs)
plot(F)
hline!([α/(2+α)])


vz = map(x->var(x.z), xs)
plot(vz)
plot!(t->2*(1+2α)*M.V*exp(-t/(4N)))



α = 0.1
N = 200
M = InfDemeMix(α=fill(α, 3), U=[0. 0.; 0. 0.; 0. 1.0])
pop = InfPop(z=rand(Normal(0, √(4M.V)), N), c=fill(4,N))

for i=1:100
    pop = generation(pop, M, z->1)
    pop.F .= 0.0
    pop.Φ .= 0.0
    @info var(pop.z), 2*(1+2α)*M.V
end
