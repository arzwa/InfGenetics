
using InfGenetics, StatsBase, Distributions, Parameters, Plots

MM = InfDemeMix(U=Umat(0.0, 0.0), V=0.5)
N  = 1000
map([2,4]) do c
    pop = InfPop(z=zeros(N), c=fill(c,N))
    for i=1:20
        pop = generation(pop, MM)
        pop.F .= 0.0
    end
    mean(pop.z), var(pop.z)
end

MM = InfDemeMix(U=Umat(0.1, 0.1), V=0.5)
N  = 1000
pop = InfPop(z=zeros(N), c=fill(2,N))
for i=1:100
    pop = generation(pop, MM)
    pop.F .= 0.0
    i % 10 != 0 && continue
    zs = pop.z[pop.c .== 4]
    @info var(zs), 4*MM.β[3]^2*MM.V
end
