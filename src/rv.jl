# approximate the reproductive value of an individual in a given population
# assuming the population composition does not change.
const cs = [2,3,3,4]

# Recursive! Exponential in `n`!
"""
    approx_rv

This gives a crude idea of the reproductive value of an individual dropped in a
population. It calculates the mean relative fitness of the individual, and then
convolves over its different possible families, computing the fitness for the
average offspring within the family... (this is not very clear) i.e. it does
the recursion:
```   
RV = W(z₀)*(Σᵢ Pᵢ W(zᵢ))
```
where Pᵢ is the probability that, conditional on the individual with trait z₀
reproducing, it does so in a family `i` with average offspring trait value zᵢ.
The recursion is broken after `n` iterates. This is not feasible for
populations of decent size and `n` larger than 4 say.
"""
function approx_rv(M, pop, individual, n=3)
    n == 0 && return 1.0
    w, W, Z = mean_relative_fitness(M, pop, individual)
    w == 0.0 && return 0.0
    P = W ./ sum(W)
    w_ = 0.0
    for i=1:size(P,1), j=1:4
        wij = approx_rv(M, pop, (z=Z[i,j], c=cs[j]), n-1)
        w_ += P[i,j]*wij
    end
    return w*w_
end

function mean_relative_fitness(M, pop, individual)
    @unpack z, c = individual
    pop_ = InfGenetics.add_unrelated_individuals(pop, [z], [c], [0.0])
    W, Z, _ = InfGenetics.fitnesses(M, pop_) 
    # ignore the variance (and selection) and go with the mean
    w = mapreduce(i->vec(mean(W[:,:,i], dims=2)), +, 1:4)
    w ./= mean(w)
    w[end], W[end,:,:], Z[end,:,:]
end
