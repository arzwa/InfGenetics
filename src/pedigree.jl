# Pedigree data structure
# - each individual has two parents (left, right)
# - we need to be able to prune the pedigree efficiently (remove all irrelevant
#   parts based on who is still alive)
# - we need to be able to efficiently find the generation in which two
#   individuals share an ancestor
# - we might wish to include ploidy level rightaway, because in the mixed
#   ploidy case, a polyploid individual may inherit asymmetrically from its
#   parents (however, this may be a parametric type thing)
# - do we want a recursive data structure? or more a la Graphs.jl?

struct PedigreeNode{T} 
    this::T
    left::Union{PedigreeNode{T},Nothing}
    rght::Union{PedigreeNode{T},Nothing}
end

Base.show(io::IO, n::PedigreeNode) = write(io, "PedigreeNode($(n.this))")

istoplevel(n::PedigreeNode) = isnothing(n.left) && isnothing(n.rght)

# monoidal
wfgeneration!(pop) = [mate(rand(pop), rand(pop)) for i=1:length(pop)]

function mate(a, b)
    k1 = length(a.this)
    k2 = length(b.this)
    y1 = sample(a.this, k1÷2, replace=false)
    y2 = sample(b.this, k2÷2, replace=false)
    PedigreeNode([y1 ; y2], a, b)
end

homozygosity(n::PedigreeNode) = homozygosity(n.this)
homozygosity(x::Vector) = sort(collect(values(countmap(hash.(x)))))

function homozygosity_tetpop(pop)
    hs = Dict([1,1,1,1]=>0., [1,1,2]=>0., [2,2]=>0., [1,3]=>0., [4]=>0.)
    for (k,v) in proportionmap(map(homozygosity, pop))
        hs[k] += v
    end
    [hs[[4]], hs[[1,3]], hs[[2,2]], hs[[1,1,2]], hs[[1,1,1,1]]]
end

function simulate(N, n)
    pop = [PedigreeNode(randn(4), nothing, nothing) for i=1:N]
    hs = map(1:n) do _
        pop = wfgeneration!(pop) 
        homozygosity_tetpop(pop)
    end
    return pop, permutedims(hcat(hs...))
end
