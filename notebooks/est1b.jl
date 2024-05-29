### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 2bf45042-10ad-464a-9329-4012d1acaadb
using Pkg; Pkg.activate("..")

# ╔═╡ 604c8252-4f64-43fb-bab8-07d7abd20c40
using InfGenetics, StatsBase, Distributions, Parameters, Plots

# ╔═╡ 49ab8458-7c71-474e-9077-7800e2d97f44
default(size=(280,220), grid=false, legend=:outertopright, framestyle=:box, 
	fontfamily="Computer modern", fg_legend=:transparent, bg_legend=:transparent, titlefont=9, guidefont=9, title_loc=:left, markerstrokewidth=0)

# ╔═╡ 3dad5ec8-c615-4811-9ae5-a10ed7286151
md"""
## Simulation functions
"""

# ╔═╡ b7b2f2cc-f2ba-4ac9-a16f-797eccfea03e
function simest(M, x=InfPop([rand(Normal(M.θ[1], √(2M.V)))]))
	xs = [x]
	while 0 < length(x.z) < 100
		x = generation(M, x)
		push!(xs, x)
	end
	return (length(x.z) > 0, xs)
end

# ╔═╡ a34aba59-c35a-4d8c-9480-13b73bcdec38
function simuntil(f) 
	while true
		x, xs = f()
		x && return xs
	end
end

# ╔═╡ 8092b6c6-e903-4915-b293-a0114f411535
function sim(f, n) 
	map(1:n) do _
		x, xs = f()
	end
end

# ╔═╡ cf762c4d-c0ac-4667-bd92-327714b6e00c
md"""
## Establishment from a single migrant
"""

# ╔═╡ ac88335c-ef4d-43e5-af0c-33e18095317a
InfDemeMixEst(U=Umat(0.0, 0.0))

# ╔═╡ b1e578d6-a2cb-44e2-81c0-db50d746e086
ps = let n=100000, M=InfDemeMixEst(U=Umat(0.0, 0.0), m=0.0)
	reps1 = sim(()->simest(M, InfPop(z=[-2.0], c=[2])), n)
	reps2 = sim(()->simest(M, InfPop(z=[-2.0], c=[4])), n)
	a = sum(first.(reps1))/n
	b = sum(first.(reps2))/n
	a, b
end

# ╔═╡ 9c76b8e8-4070-4b94-8fed-b5f52704611e
reps = map(0.25:0.05:1) do b 
	n = 250000
	M = InfDemeMixEst(
		U=Umat(0.0, 0.0),
		γ=0.25,  # as in BE18
		m=0.0,   # no migration!
		V=0.5,   # scale relative to √2V = 1
		β=[1.0, 1.0, √b]
	)
	reps2 = sim(()->simest(M, InfPop(z=[-2.0], c=[4])), n)
	√b, sum(first.(reps2))/n
end

# ╔═╡ 11915cc1-6ba3-431d-8797-a75fed6d1b85
P1 = begin
	scatter(first.(reps), last.(reps) ./ ps[1], legend=false, color=:black, ms=2.5,
		ylabel="\$P_4/P_2\$", xlabel="\$\\beta\$", ylim=(0,10), size=(260,210))
	hline!([1], color=:firebrick, 
		title="(B) \$z_0=-2, \\gamma=0.25\$")
	vline!([√(0.5)], color=:gray, ls=:dash)
end #; savefig("../doc/img/est1.pdf")

# ╔═╡ cd0de25f-9068-479c-8999-2a8f266e4374
X = let gs = [0.0 ; 10 .^ range(log10(0.02), stop=log10(4), length=10)]
	map(gs) do g
		n = 5000
		M = InfDemeMixEst(
			U=[1.0 0.0 ; 0.0 0.0 ; 0.0 1.0],
			γ=g, 
			m=0.0,   # no migration!
			V=0.5,   # scale relative to √2V = 1
		)
		reps = map([2,4]) do c
			xs = sim(()->simest(M, InfPop(z=[-1.0], c=[c])), n)
			sum(first.(xs))/n
		end
		g, reps
	end
end

# ╔═╡ 65bf38a3-f91d-42da-bf7b-d414cbf47e1a
P2 = let i=1
	plot(first.(X)[i:end], first.(last.(X))[i:end], label="diploid", 
		marker=true, ms=2, title="(A) \$z_0=-2, \\beta=\\sqrt{1/2}\$")
	plot!(first.(X)[i:end], last.(last.(X))[i:end], label="tetraploid", 
		legend=:topright, xlabel="\$\\gamma\$", size=(300,240), 
		ylabel="\$P\$", marker=true, ms=2)
	# hline!([X[1][2]])
end

# ╔═╡ 7b525da0-45f2-4282-8af3-f7df81546e85
P3 = let i=1, xs=first.(X), ys4=last.(last.(X)), ys2=first.(last.(X))
	rs = ys4 ./ ys2
	plot(first.(X)[i:end], first.(last.(X))[i:end], label="diploid", 
		marker=true, ms=2, title="(A) \$z_0=-2, \\beta=\\sqrt{1/2}\$")
	plot!(first.(X)[i:end], last.(last.(X))[i:end], label="tetraploid", 
		legend=:topright, xlabel="\$\\gamma\$", size=(300,240), 
		ylabel="\$P\$", marker=true, ms=2)
	plot!([],[], color=:black, label="\$P_4/P_2\$", ls=:dash)
	plot!(twinx(), xs[i:end], rs[i:end], label="", legend=false, 
		ylabel="\$P_4/P_2\$", ls=:dash,
		color=:black, alpha=0.4, 
		marker=true, ms=2, ylim=(0,5))
end

# ╔═╡ 827f7fb8-15d8-49dc-9908-c2f7e77923d2
plot(P2, P1, size=(2*270, 210), fg_legend=:transparent, margin=3Plots.mm)#; savefig("../doc/img/est2.pdf")

# ╔═╡ Cell order:
# ╠═2bf45042-10ad-464a-9329-4012d1acaadb
# ╠═604c8252-4f64-43fb-bab8-07d7abd20c40
# ╠═49ab8458-7c71-474e-9077-7800e2d97f44
# ╟─3dad5ec8-c615-4811-9ae5-a10ed7286151
# ╠═b7b2f2cc-f2ba-4ac9-a16f-797eccfea03e
# ╠═a34aba59-c35a-4d8c-9480-13b73bcdec38
# ╠═8092b6c6-e903-4915-b293-a0114f411535
# ╠═cf762c4d-c0ac-4667-bd92-327714b6e00c
# ╠═ac88335c-ef4d-43e5-af0c-33e18095317a
# ╠═b1e578d6-a2cb-44e2-81c0-db50d746e086
# ╠═9c76b8e8-4070-4b94-8fed-b5f52704611e
# ╠═11915cc1-6ba3-431d-8797-a75fed6d1b85
# ╠═cd0de25f-9068-479c-8999-2a8f266e4374
# ╠═65bf38a3-f91d-42da-bf7b-d414cbf47e1a
# ╠═7b525da0-45f2-4282-8af3-f7df81546e85
# ╠═827f7fb8-15d8-49dc-9908-c2f7e77923d2
