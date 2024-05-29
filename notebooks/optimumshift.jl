# Simulate a shifting optimum in a mixed-ploidy deme

function simshift(M, N, Vs, θ, n1=100, n2=100)
    pop = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
    # equilibrate, assume the population is very large
    ys = [(z=pop.z, c=counts(pop.c, 2:4))]
    for _=1:n1
        pop = generation(pop, M, z->exp(-γ*z^2))
        pop.F .= 0.0
        pop.Φ .= 0.0
        push!(ys, (z=pop.z, c=counts(pop.c, 2:4)))
    end
    for _=1:n2
        pop = generation(pop, M, z->exp(-(z - θ)^2/(2Vs)))
        push!(ys, (z=pop.z, c=counts(pop.c, 2:4)))
    end
    return ys
end

function tetraploid_est(M, N, Vs, θ, n1=100, n2=100, nmax=100)
    k = 0
    while true
        ys = simshift(M, N, Vs, θ, n1, n2)  
        (ys[end].c[3] > ys[end].c[1]) && return (k+1, ys)
        k > nmax && return (k, NaN)
        k += 1
    end
end

# There appears to be a very sharp threshold for establishment that is a
# function of θ, Vz, Vs, β, u and v (I guess). This is probably something one
# could try to predict. Note however that as θ/√Vs becomes large enough, we get
# establishment even when the different ploidy levels are phenotypically
# indistinguishable. There are only two differences between ploidy levels in
# that case: the different rate of inbreeding and the different gametic
# outputs.
N = 200
θ = 9
M = InfDemeMix()
Vs = 2*2M.V   # assume Vs >> Vz

tw, xs2 = tetraploid_est(M, N, Vs, θ, 50, 20)

xs1 = simshift(M, N, Vs, θ, 100, 100)

[exp(-θ/(2Vs) + (1/(2Vs))^2/2*m*M.β[m-1]^2*M.V) for m=2:4]



