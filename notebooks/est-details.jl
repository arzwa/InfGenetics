
M = InfDemeMixEst(
    U=Umat(0.0,0.0),  # no unreduced gametes 
    γ=0.25,           # selection strength
    m=0.0,            # no migration
    θ=fill(1.5, 3),   # maladaptation
    V=0.5)            # segregation variance

Random.seed!(446)
pop = InfPop(z=zeros(10), c=fill(2, 10))
plot()
for i=1:15
    pop = generation(M, pop)
    @info length(pop.z)
    scatter!(fill(i, length(pop.z)), pop.z, color=:black, label="", ms=2)
end
plot!()
