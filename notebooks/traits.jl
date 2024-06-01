using InfGenetics, Distributions, Plots

N = 2000
M = InfDemeMix(U=Umat(1.0, 0.0))
pop = InfPop(z=rand(Normal(-2, √(2M.V)), N), c=fill(2,N))
gen1 = generation(pop, M)

stephist(pop.z, fill=true, fillalpha=0.2)
stephist!(gen1.z, fill=true, fillalpha=0.2)


