using Distributions, Plots, StatsPlots, StatsBase

V = 1.0   # segregation variance
N = 100   # population size
γ = 0.25  # selection strength
θ = 2.0   # maladaptation

# Current population trait values
z = rand(Normal(0, √2V), N)

# Number of replicate offspring populations to simulate
nrep = 2500

# Number of individuals in next generation is NW̄
sims1 = map(1:nrep) do _
    Noff = 10000
    off = map(1:Noff) do _
        z1 = rand(z)
        z2 = rand(z)
        zoff = rand(Normal((z1 + z2)/2, √V))
        woff = exp(γ*(zoff - θ))
        woff, zoff
    end
    w̄ = mean(first.(off))
    𝔼N = N*w̄
    N′ = rand(Poisson(𝔼N))
    z′ = sample(last.(off), Weights(first.(off)), N′, replace=false)
end

# Calculate expected contribution of each family, sample families and sample
# offspring trait values from the Gaussian after selection within each family
sims2 = map(1:nrep) do _
    W = [exp(γ*((z[i]+z[j])/2 - θ) + γ^2/2 * V) for i=1:N, j=1:N]
    𝔼N = sum(W)/N
    N′ = rand(Poisson(𝔼N))
    C = CartesianIndices((1:N, 1:N))
    idx = sample(C, Weights(vec(W)), N′)
    z′ = map(idx) do ix
        i, j = Tuple(ix)
        z̄ = (z[i] + z[j])/2 + γ*V
        rand(Normal(z̄, √V))
    end
end

P1 = density(vcat(sims1...), label="sims1", xlabel="\$z'\$", title="(A)")
density!(vcat(sims2...),label="sims2")
z1 = map(mean, sims1)
z2 = map(mean, sims2)
P2 = density([z1 z2], xlabel="\$\\overline{z}'\$", legend=false, title="(B)") 
v1 = map(std, sims1)
v2 = map(std, sims2)
P3 = density([v1 v2], xlabel="\$\\sqrt{V′}\$", legend=false, title="(C)")
P4 = stephist(map(length, sims1), xlabel="\$N'\$", legend=false, title="(D)")
stephist!(map(length, sims2))
plot(P1,P2,P3,P4,size=(900,200),bottom_margin=5Plots.mm,layout=(1,4))


# With mixed-ploidy, the first way of simulating would have to include all
# failed matings with fitness 0, i.e. if triploid produce 2ν euploid gametes,
# and 1-2ν aneuploid/polyploid ones, these would yield zero-fitness offspring,
# which would have to be included among the `Noff`.

