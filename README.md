# InfGenetics

This repository contains code for conducting simulations of the mixed-ploidy
infinitesimal model.

The `notebooks/` directory contains the code for generating the figures in
Zwaenepoel (2025) (see especially the files `establishment.jl` and
`estwmigration.jl`, associated with figure 2 and figures 3,4 & 5 respectively).

## Basic usage examples

### Finite mixed-ploidy WF population

```julia
# Define the model (use default α, β, V, etc.)
M = InfDemeMix(U=[0.92 0.08; 0.08 0.08; 0.0 0.92]) 

# Initial population of size N=500
N = 500
p = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))

# Simulate the model for 1000 generations, assuming neutrality 
# (fitness function is `z->1`)
xs = fiterate((x,i)->generation(x, M, z->1), p, 1000)
```

### Establishment model

```julia
using InfGenetics

# Define the model
M = InfDemeMixEst(
    U=Umat(0.05,0.05),  # 5% unreduced gametes/euploid gametes in triploids 
    γ=0.25,             # selection strength
    m=0.01,             # migration rate
    θ=fill(2.5, 3),     # maladaptation
    V=0.5)              # segregation variance
    
# Simulate until establishment, storing the population data at each
# generation
sim = simest(M) 

# Simulate until establishment, and don't store the population data at each
# generation
sim = simest2(M) 
```

