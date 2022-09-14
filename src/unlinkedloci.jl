# Simulation methods for a genome of a large but finite number of additive
# unlinked loci with small effect under Mendelian inheritance.
#
# Consider the model which is also discussed in Barton et al. in an
# example there where we have n loci, each of effect ± 1/√n.
# generalize to k-ploids -> ± 1/√(kn)

const Genome{T} = Matrix{T}
const Population{T} = Vector{Genome{T}}

ploidy(g::Genome) = size(g)[2]

# ϕ denotes a vector of homozygosity states in tetraploids
# ϕ = [abcd, aabc, aabb, aaab, aaaa] in tets

# compute inbreeding coefficient from homozygosity coefficients
F(ϕ) = ϕ[2]/6 + ϕ[3]/3 + ϕ[4]/2 + ϕ[5]

# compute the variance reduction factor from homozygosity coefficients
Φ(ϕ) = first([1 5/6 2/3 1/2 0] * ϕ)  # this is 1 - F(ϕ)...

# generate a random genome with a given probability distribution `ϕ` over
# homozygosity states. `V` is the phenotypic variance of the non-inbred base
# population, i.e. twice the segregation variance.
function randgenome(nloci, ploidy, ϕ, V)
    genotypes = map(i->sample(1:length(ϕ), Weights(ϕ)), 1:nloci)
    xs = if ploidy == 4
        mapreduce(i->randlocustet(genotypes[i], V), hcat, 1:nloci)
    elseif ploidy == 2
        mapreduce(i->randlocusdip(genotypes[i], V), hcat, 1:nloci)
    elseif ploidy == 3
        mapreduce(i->randlocustrip(genotypes[i], V), hcat, 1:nloci)
    else
        @error "should implement"
    end
    return permutedims(xs) ./ √(nloci*ploidy)
end

function randlocustet(genotype, V)
    x = rand(Normal(0, √V), 4)
    (genotype == 2 || genotype == 3) && (x[1] = x[2])
    genotype == 3 && (x[3] = x[4])
    genotype == 4 && (x[1] = x[2] = x[3])
    genotype == 5 && (x[1] = x[2] = x[3] = x[4])
    return x
end

function randlocusdip(genotype, V)
    x = rand(Normal(0, √V), 2)
    genotype == 2 && (x[1] = x[2])
    return x
end

function randlocustrip(genotype, V)
    x = rand(Normal(0, √V), 3)
    genotype == 2 && (x[1] = x[2])
    genotype == 3 && (x[1] = x[2] = x[3])
    return x
end

function gamete(g::Genome, α=0.)
    k = ploidy(g) ÷ 2
    h = similar(g, size(g, 1), k)
    for i=1:size(g, 1)
        h[i,:] .= rand() < α ? rand(g[i,:]) : sample(g[i,:], k, replace=false)
    end
    h
end

function ugamete(g::Genome, ξ=0.5)
    k = ploidy(g)
    @assert k == 2 "currently only for 2n=2"
    h = similar(g, size(g)...)
    for i=1:size(g, 1)
        h[i,:] .= rand() < ξ ? rand(g[i,:]) : sample(g[i,:], k, replace=false)
    end
    h
end

# segregate a triploid in a diploid and haploid gamete
function tripdipgamete(g::Genome, ξ=0.5)
    k = ploidy(g)
    @assert k == 3 "need a triploid"
    h = similar(g, size(g,1), 2)
    for i=1:size(g, 1)    
        h[i,:] .= rand() < ξ ? rand(g[i,:]) : sample(g[i,:], 2, replace=false)
    end
    h
end

function triphapgamete(g::Genome)
    k = ploidy(g)
    @assert k == 3 "need a triploid"
    h = similar(g, size(g,1), 1)
    for i=1:size(g, 1)    
        h[i,:] .= rand(g[i,:])
    end
    h
end

function offspring(g::Genome, h::Genome, α=0.)
    o = similar(g)
    k = ploidy(g) ÷ 2
    for i=1:size(g, 1)
        o[i,1:k]     .= rand() < α ? rand(h[i,:]) : sample(g[i,:], k, replace=false)
        o[i,k+1:end] .= rand() < a ? rand(h[i,:]) : sample(h[i,:], k, replace=false)
    end
    Genome(o)
end

function wf_gen(pop::Population)
    map(i->offspring(rand(pop), rand(pop)), 1:length(pop))
end

function wf(pop, n, callbacks=[mean, var])
    zs = sum.(pop)
    xs = [map(f->f(zs), callbacks)]
    for i=1:n
        i % 10 == 0 && (@info "generation $i")
        pop = wf_gen(pop)
        zs = sum.(pop)
        fs = map(f->f(zs), callbacks)
        push!(xs, fs)
    end
    return zs, xs
end
