
# additive genetic variance in HWLE
function addvarhw(k, a, p)
    Ex  = [i*a*binomial(k,i)*p^i*(1-p)^(k-i) for i=0:k]
    Ex² = [(i*a)^2*binomial(k,i)*p^i*(1-p)^(k-i) for i=0:k]
    sum(Ex²) - sum(Ex)^2, k*a^2*p*(1-p)
end

# segregation variance in autotetraploids
# there are four relevant inbreeding coefficients (5 states of homozygosity)

# Here I estimate Eₓ[Var[Y|X]] by Monte Carlo
using Random, StatsBase
for genotype=[1,2,3,4], α=[0.,0.1,0.2]
    vs = map(1:10000) do _
        x = randn(4) 
        (genotype == 2 || genotype == 3) && (x[1] = x[2])
        genotype == 3 && (x[3] = x[4])
        genotype == 4 && (x[1] = x[2] = x[3])
        genotype == 5 && (x[1] = x[2] = x[3] = x[4])
        ys = map(1:1000) do _
            if rand() < α
                2*rand(x)
            else
                y = shuffle(x)[1:2]
                y[1] + y[2]
            end
        end
        var(ys)
    end
    r = [0., 1/6, 1/3, 1/2, 1.][genotype]
    pred = (1-r)*(1+2α)
    @info (genotype, α) (mean(vs), pred)
end

using Distributions
v = 2.5
map(1:10000) do _
    x = rand(Normal(0,√v), 4)
    var([rand(x) for i=1:10000])
end |> mean
3/4 * v


F(ϕ) = ϕ[2]/6 + ϕ[3]/3 + ϕ[4]/2 + ϕ[5]
Φ(ϕ) = first([1 5/6 2/3 1/2 0] * ϕ)

# for a given probability distribution over homozygosity states, simulate a
# random set of unlinked loci, and generate gametes by Mendelian segregation.
function manyloci(ϕ, nloci, noff, nrep)
    vs = map(1:nrep) do rep1
        # draw the homozygosity states for each locus
        genotypes = map(i->sample(1:5, Weights(ϕ)), 1:nloci)
        # draw the allelic effects, v0/2 = 1
        xs = randn(4, nloci) ./ √nloci
        ys = map(1:noff) do rep2  
            # generate gametes
            z = map(1:nloci) do i
                genotype = genotypes[i]
                x = xs[:,i]
                (genotype == 2 || genotype == 3) && (x[1] = x[2])
                genotype == 3 && (x[3] = x[4])
                genotype == 4 && (x[1] = x[2] = x[3])
                genotype == 5 && (x[1] = x[2] = x[3] = x[4])
                y = shuffle(x)[1:2]
                y[1] + y[2]  # gamete
            end |> sum
        end
        var(ys)
    end
    mean(vs), Φ(ϕ), (1-F(ϕ))
end
        
ϕ = [0.1, 0.54, 0.07, 0.25, 0.04]  # this is from an actual WF simulation
manyloci(ϕ, 100, 100, 1000)

# E[V[Y|X]] = VY - V[E[Y|X]] = VX - V[(X1 + X2 + X3 + X4)/4] = VX - 4VX/16 = 3/4VX

# segregation variance (Var[Y|X]) is half of the total variance.

# the proportional reductions in segregation variance are
# 0, 1/6, 1/3, 1/2 and 1
# for increasing homozygosity

# double reduction leads to an increase in segregation variance (there are more
# distinct gametic genotypes transmitted, i.e. 10 instead of 6).
#
# double reduction leads to faster inbreeding though, so reduction in
# segregation variance in the long term.

using QuadGK, Distributions, Plots

fitfun(θ, Vs) = z->exp(-0.5*(θ-z)^2/Vs)

θ = 1.7
Vs = 1.5
w = fitfun(θ, Vs)
z̄ =  -1.62
V0 = 1.3
plot(w, color=:lightgray, fill=true)
Z, _ = quadgk(z->w(z)*pdf(Normal(z̄, √V0), z), -10., 10.)
plot!(z->pdf(Normal(z̄, √V0), z), color=:black, ls=:dot, lw=2)
plot!(z->(1/Z)*w(z)*pdf(Normal(z̄, √V0), z), color=:black, lw=2)

ϕ = 1/((1/Vs) + (1/V0))
ψ = (z̄/V0 + θ/Vs)*ϕ
plot!(z->pdf(Normal(ψ, √ϕ), z), color=:salmon, fill=true, fillalpha=0.2, lw=0)

# so its a 'posterior' (precision weighted mean)
# Now the expected fitness? should be some function of the normalizing constant
# wolfram alpha -(1/2)*(a^2/V + b^2/W) + (1/2)*(a/V + b/W)^2 * (1/((1/V) + (1/W)))
function theZ(z̄, θ, V0, Vs)
    C = √(Vs/(Vs + V0))
    D = -(θ - z̄)^2/(2*(Vs + V0))
    C*exp(D)
end
theZ(z̄, θ, V0, Vs)  # correct!

