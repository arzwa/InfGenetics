using Pkg; Pkg.activate(@__DIR__); using InfGenetics
using InfGenetics: InfDemeMix, offspring_distribution, 𝔼fitness, generation
using InfGenetics: StabilizingSelection, selection
using Distributions, StatsBase, LaTeXStrings, Plots, PlotThemes
theme(:wong2, titlefont=10, title_loc=:left)

# plot a series of populations
plotpops(pops) = 
    plot(
         plot(map(x->mean(x[2].z), pops), title=L"\bar{z}"),
         plot(map(x->var(x[2].z), pops), title=L"V_z"),
         plot(permutedims(mapreduce(x->counts(x[2].m, 2:4), hcat, pops)), title=L"N"),
         plot(map(x->mean(x[2].F), pops), title=L"\bar{F}", ylim=(0,1)),
         legend=false,
        )

# simulates an adaptive challenge. We should have a branching process
# approximation to get a rough probability of establishment. Then check how
# often tets/dips do establish?
d = InfDemeMix(z = randn(100), 
               m = fill(2, 100),
               U = [0.9 0.1 ; 0.05 0.05; 0.0 0.95], 
               w = [1.0, 1.0, 1.0],
               β = [1.0, 1.5, 2.0], 
               α = 0.2, 
               ξ = 0.5)

context = StabilizingSelection(θ=0., ω=1.0, r=1.1, K=100., ν=2., z̄=-5., V=4d.V)

# comparing dips/tets
pops2 = fiterate((d,_)->InfGenetics.generation(d, context), d(U=[1.0 0.;0. 0.; 0. 0.]), 500)
pops4 = fiterate((d,_)->InfGenetics.generation(d, context), d(m=fill(4,100)), 500)
ps = plot(map(plotpops, [pops2, pops4])..., size=(800,300))

context = StabilizingSelection(θ=0., ω=1.0, r=1.1, K=500., ν=2., z̄=0., V=2d.V)
sims = map(1:100) do i
    pops = fiterate((d,i)->InfGenetics.generation(d, context(θ=i>50 ? 3. : 0.)), d, 100)
    ms = counts(pops[end][2].m, 1:4)
    @info ms
    (argmax(ms), pops)
end;

plotpops(rand(sims)[2])

pops = rand(filter(x->x[1]==4, sims))[2]
plotpops(pops)

pops = fiterate((d,i)->InfGenetics.generation(d, context(θ=3.)),
                pops[end][2], 500)
plotpops(pops)

# ------
# should generate an N × N × 4 array with for each parental pair and each
# possible gametic contribution an offspring distribution object
# with directional selection one just has to increment the offspring mean in
# order to compute the phenotypic distribution after selection.

Y = [offspring_distribution(d, i, j) for i=1:length(d), j=1:length(d)]
Y = directional_selection!.(Y, 0.)
w = 𝔼fitness.(Y)
N = rand(Poisson(mean(w) * length(d)))
P = sample(vec(Y), Weights(vec(w)), N)
O = InfGenetics.randoff.(P, Ref(d))
mean(first.(O))

# after selection
# w(z) * p(z|zi,zj) = exp(b*z) * exp(-(z - m)^2/v)
#                   = exp((b*z*v - z^2 -2zm + m^2)/v)

zfun(b, zi, zj, vij) = z->exp(b*z)*exp(-(z - (zi+zj)/2)^2/(2vij))

b, zi, zj, vij = 0.25, -0.3, 0.2, 1.
zf = zfun(b, zi, zj, vij)
Z, _ = quadgk(zf, -100., 100.)
plot(z->zf(z)/Z)
plot!(z->pdf(Normal((zi + zj)/2 + b*vij, √vij), z))
plot!(z->pdf(Normal((zi + zj)/2, √vij), z))


