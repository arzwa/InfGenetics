# Simulation methods for a genome of a large but finite number of additive
# unlinked loci with small effect under Mendelian inheritance.

const Genome{T} = Matrix{T}
const Population{T} = Vector{Genome{T}}

ploidy(g::Genome) = size(g)[2]

randgenome(n, k) = Genome(randn(n, k) ./ √(n*k)) 
randgenome(n, k, N) = [Genome(randn(n, k) ./ √(n*k)) for i=1:N]

generate_offspring(pop::Population) =
    map(i->generate_offspring(rand(pop), rand(pop)), 1:length(pop))

function generate_offspring(g::Genome, h::Genome)
    o = similar(g)
    k = ploidy(g) ÷ 2
    for i=1:size(g, 1)
        o[i,1:k] = sample(g[i,:], k, replace=false)
        o[i,k+1:end] = sample(h[i,:], k, replace=false)
    end
    Genome(o)
end

generate_offspring(g::Genome, h::Genome, n) = 
    map(_->generate_offspring(g, h), 1:n)

family(g, h, n) = sum.(generate_offspring(g, h, n))

# evolve a population and get an estimate of the segregation variance in each
# generation
function evolve_pop_track_segvar(pop, t; famsize=100)
    map(1:t) do gen
        pop = generate_offspring(pop)
        var(family(rand(pop), rand(pop), famsize))
    end
end
