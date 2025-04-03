### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ baef8376-1514-11ef-1d35-277242ce9bf3
using Pkg; Pkg.activate("..")

# ╔═╡ 019ac4ed-1679-4683-8d8c-04217475d2f8
using InfGenetics, StatsBase, Distributions, Parameters, Plots

# ╔═╡ 117ed0b3-7d8a-42f7-92aa-94d59d236f39
default(size=(280,220), grid=false, legend=:outertopright, framestyle=:box, 
	fontfamily="Computer modern",  bg_legend=:transparent, titlefont=9, guidefont=9, title_loc=:left, markerstrokewidth=0)

# ╔═╡ cfb8cc81-7801-476a-a2b9-a88286e6c9d7
md"""
## Equilibrium variance
"""

# ╔═╡ 247e5eea-2ba4-4792-80a2-eef8e6d8acdd
let MM = InfDemeMix(U=[1.0 0.0 ; NaN NaN; 0.0 1.0], V=0.5)
	pop = InfPop(z=zeros(500), c=fill(4,500))
	for i=1:20
		pop = generation(pop, MM)
		pop.F .= 0.0
	end
	mean(pop.z), var(pop.z)
end

# ╔═╡ 289ed1fc-c191-412c-b045-42d59c636fef
md"""
## Establishment of a single migrant individual

This is the first thing Barton & Etheridge examined: finding the probability of non-extinction when starting with a single individual randomly drawn from the source population. This sounds simple, but it does not appear analytically tractable.

>While the population size is small, inbreeding rapidly reduces the genetic variance, and so if the population is to survive, individuals with sufficiently large $z$ must be produced in the first few generations. The probability that this will happen depends only on the two dimensionless quantities $\frac{z_0}{\sqrt{2V}}$ and on $\beta\sqrt{2V}$.

They note that $\beta \sqrt{2V}$ is the standard deviation of log fitness of the source population in the new habitat. They set $\sqrt{2V}=1$ and interpret the selectoni gradient and maladaptation of the source as scaled relative to this.
"""

# ╔═╡ feea241c-9248-4d0b-a9af-dcf200a84ce3
md"""
### Strictly diploid model
"""

# ╔═╡ bf0f6791-638b-446b-bbb1-5d37a2e10934
M1 = InfDemeMixEst(
	U=[1.0 0.0 ; NaN NaN ; NaN NaN],
	γ=0.25,  # as in BE18
	m=0.0,   # no migration!
	V=0.25,   # scale relative to √2V = 1
	θ=[-2.0, NaN, NaN],
	c=[1.0, 0.0, 0.0]
)

# ╔═╡ c76e2724-c443-43c6-bea8-fcd416fed746
function simest1(M, x=InfPop([rand(Normal(M.θ[1], √(2M.V)))]))
	xs = [x]
	while 0 < length(x.z) < 100
		x = generation(x, M)
		push!(xs, x)
	end
	return (length(x.z) > 0, xs)
end

# ╔═╡ eb72dab7-b78e-44d5-a397-0c3876a11b56
function simuntil(f) 
	while true
		x, xs = f()
		x && return xs
	end
end

# ╔═╡ f196ce4c-ef7b-4f52-b90e-fca4a2baf98e
function sim(f, n) 
	map(1:n) do _
		x, xs = f()
	end
end

# ╔═╡ ee125d9e-b3a3-4110-8985-1ae7f4bfa8a8
xs = simuntil(()->simest1(reconstruct(M1, γ=0.25), InfPop(z=[-2.0])))

# ╔═╡ 980d4eae-ab5c-4498-94c6-c3de4646d147
plot(); map(enumerate(xs)) do (i,x)
	scatter!(fill(i, length(x.z)), ms=2, x.z, color=:cornflowerblue, legend=false)
end ; plot!(xlabel="generation", ylabel="\$z\$")

# ╔═╡ e3a3e575-442a-45f0-9732-6d86dcc23851
scatter(map(x->length(x.z), xs), map(x->mean(x.z), xs),
	ms=2, color=:cornflowerblue, legend=false, xscale=:log10, 
	xlabel="\$N\$", ylabel="\$\\overline{z}\$")

# ╔═╡ 1e15d33b-6ecb-49b7-bb04-2e63983214c6
md"""
### Single tetraploid migrant from diploid source
"""

# ╔═╡ c316ce4b-589e-4f0a-97a1-1c29405fea0d
M2 = InfDemeMixEst(
	U=[1.0 0.0 ; 0.0 0.0 ; 0.0 1.0],
	γ=0.25,  # as in BE18
	m=0.0,   # no migration!
	V=0.5,   # scale relative to √2V = 1
	θ=[-2.0, NaN, NaN],
	c=[1.0, 0.0, 0.0]
)

# ╔═╡ 19c6212c-c8d3-42d3-b195-fadac24e6835
xs2 = simuntil(()->simest1(M2, InfPop(z=[-2.0], c=[4])))

# ╔═╡ 27911e3d-890a-4ab1-a32a-ea70cdc87671
plot(); map(enumerate(xs2)) do (i,x)
	scatter!(fill(i, length(x.z)), ms=2, x.z, color=:cornflowerblue, legend=false)
end ; plot!(xlabel="generation", ylabel="\$z\$")

# ╔═╡ 36ca5690-c433-420d-ad5b-52fbd8ba5b56
let 
	P1 = plot(map(x->1-mean(x.F), xs)); 
	plot!(map(t->prod(map(x->1 - 1/(2length(x.z)), xs[1:t])), 0:length(xs)-1))
	P2 = plot(map(x->1-mean(x.F), xs2)); 
	plot!(map(t->prod(map(x->1 - 1/(4length(x.z)), xs2[1:t])), 0:length(xs2)-1))
	plot(P1, P2, size=(500,220), legend=false, ylim=(0,1))
end

# ╔═╡ 67e54219-d003-42d4-b934-c3257db56e30
md"""
### Establishment probability from a single migrant
"""

# ╔═╡ 0ea00635-fade-48d3-b971-f910a200b0b9
ps = let n=500_000
	reps1 = sim(()->simest1(reconstruct(M1, V=0.5), InfPop(z=[0.0], c=[2])), n)
    reps2 = sim(()->simest1(M2, InfPop(z=[0.0], c=[4], θ=[2.0,2.0,2.0])), n)
	a = sum(first.(reps1))/n
	b = sum(first.(reps2))/n
	a, b
end

# ╔═╡ 660de9f3-085c-4710-8372-9993559c8dce
md"""
Note that in these simulations, we use the scaling yielding equal equilibrium variances, so the > 3 fold increase in establishment probability is due to reduced inbreeding alone.
"""

# ╔═╡ d0f68cf7-7348-4482-a074-fae8bf2fdce8
reps = map(0.25:0.05:1) do b 
	n = 250000
	M = InfDemeMixEst(
		U=[1.0 0.0 ; 0.0 0.0 ; 0.0 1.0],
		γ=0.25,  # as in BE18
		m=0.0,   # no migration!
		V=0.5,   # scale relative to √2V = 1
		β=[1.0, NaN, √b]
	)
	reps2 = sim(()->simest1(M, InfPop(z=[-2.0], c=[4])), n)
	√b, sum(first.(reps2))/n
end

# ╔═╡ 2433169a-b2c3-447f-b53c-a68e3fe2196f
P1 = begin
	scatter(first.(reps), last.(reps) ./ ps[1], legend=false, color=:black, ms=2.5,
		ylabel="\$P_4/P_2\$", xlabel="\$\\beta\$", ylim=(0,10), size=(260,210))
	hline!([1], color=:firebrick, 
		title="(B) \$z_0=-2, \\gamma=0.25\$")
	vline!([√(0.5)], color=:gray, ls=:dash)
end #; savefig("../doc/img/est1.pdf")

# ╔═╡ f90c5e5d-fdf5-41a6-866c-66b66cce9d5a
let n = 200000
	M = InfDemeMixEst(
		U=[0.95 0.05 ; 0.05 0.05 ; 0.00 0.95],
		γ=0.25,  # as in BE18
		m=0.0,   # no migration!
		V=0.5,   # scale relative to √2V = 1
	)
	reps = sim(()->simest1(M, InfPop(z=[-2.0], c=[2])), n)
	sum(first.(reps))/n
end

# ╔═╡ 9b0907c9-e081-4398-a18f-fe3dec021e1c
X = let gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(2), length=10)]
	map(gs) do g
		n = 500000
		M = InfDemeMixEst(
			U=[1.0 0.0 ; 0.0 0.0 ; 0.0 1.0],
			γ=g, 
			m=0.0,   # no migration!
			V=0.5,   # scale relative to √2V = 1
		)
		reps = map([2,4]) do c
			xs = sim(()->simest1(M, InfPop(z=[-2.0], c=[c])), n)
			sum(first.(xs))/n
		end
		g, reps
	end
end

# ╔═╡ 48fad4ed-1fde-4c4a-b66e-b9a6f4283721
P2 = let i=1
	plot(first.(X)[i:end], first.(last.(X))[i:end], label="diploid", 
		marker=true, ms=2, title="(A) \$z_0=-2, \\beta=\\sqrt{1/2}\$")
	plot!(first.(X)[i:end], last.(last.(X))[i:end], label="tetraploid", 
		legend=:topright, xlabel="\$\\gamma\$", size=(300,240), 
		ylabel="\$P\$", marker=true, ms=2)
	# hline!([X[1][2]])
end

# ╔═╡ fa99b92b-1282-415b-b6dc-3d07b698cbb0
plot(P2, P1, size=(2*270, 210))#; savefig("../doc/img/est2.pdf")

# ╔═╡ e4cb18e9-e05d-440d-b78b-52df8ae7d3cb
let xs=first.(X), ys4=last.(last.(X)), ys2=first.(last.(X))
	rs = ys4 ./ ys2
	i  = 1
	plot(xs[i:end], rs[i:end], color=:black, marker=true, ms=2, ylim=(0,5))
end

# ╔═╡ 045ec5b0-a812-4887-9af7-238dde9866d7
md"""
When the selection gradient is too weak, we get a maladaptive sink? Genetic variation becomes exhausted before the trait mean becomes positive.

Is this also the case with the selection scheme of BE, where an infinite number of juveniles are created upon which selection acts?

Should implement that version as well? It involves more bookkeeping, but since the trait values within the offspring of the same ploidy level within a family is still Gaussian, it is possible.
"""

# ╔═╡ Cell order:
# ╠═baef8376-1514-11ef-1d35-277242ce9bf3
# ╠═019ac4ed-1679-4683-8d8c-04217475d2f8
# ╠═117ed0b3-7d8a-42f7-92aa-94d59d236f39
# ╟─cfb8cc81-7801-476a-a2b9-a88286e6c9d7
# ╠═247e5eea-2ba4-4792-80a2-eef8e6d8acdd
# ╟─289ed1fc-c191-412c-b045-42d59c636fef
# ╟─feea241c-9248-4d0b-a9af-dcf200a84ce3
# ╠═bf0f6791-638b-446b-bbb1-5d37a2e10934
# ╠═c76e2724-c443-43c6-bea8-fcd416fed746
# ╠═eb72dab7-b78e-44d5-a397-0c3876a11b56
# ╠═f196ce4c-ef7b-4f52-b90e-fca4a2baf98e
# ╠═ee125d9e-b3a3-4110-8985-1ae7f4bfa8a8
# ╠═980d4eae-ab5c-4498-94c6-c3de4646d147
# ╠═e3a3e575-442a-45f0-9732-6d86dcc23851
# ╟─1e15d33b-6ecb-49b7-bb04-2e63983214c6
# ╠═c316ce4b-589e-4f0a-97a1-1c29405fea0d
# ╠═19c6212c-c8d3-42d3-b195-fadac24e6835
# ╠═27911e3d-890a-4ab1-a32a-ea70cdc87671
# ╠═36ca5690-c433-420d-ad5b-52fbd8ba5b56
# ╟─67e54219-d003-42d4-b934-c3257db56e30
# ╠═0ea00635-fade-48d3-b971-f910a200b0b9
# ╟─660de9f3-085c-4710-8372-9993559c8dce
# ╠═d0f68cf7-7348-4482-a074-fae8bf2fdce8
# ╠═2433169a-b2c3-447f-b53c-a68e3fe2196f
# ╠═f90c5e5d-fdf5-41a6-866c-66b66cce9d5a
# ╠═9b0907c9-e081-4398-a18f-fe3dec021e1c
# ╠═48fad4ed-1fde-4c4a-b66e-b9a6f4283721
# ╠═fa99b92b-1282-415b-b6dc-3d07b698cbb0
# ╠═e4cb18e9-e05d-440d-b78b-52df8ae7d3cb
# ╟─045ec5b0-a812-4887-9af7-238dde9866d7
