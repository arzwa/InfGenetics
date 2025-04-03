### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 99fd9810-18e3-11ef-3ccf-0f13790a9f6e
using Pkg; Pkg.activate("..")

# ╔═╡ 65220b6c-c975-4a05-8377-19cfb88f92d8
using InfGenetics, StatsBase, Distributions, Parameters, Plots

# ╔═╡ 97808c41-877f-454b-a61a-9b3160a176e4
default(size=(280,220), grid=false, 
	fg_legend=:transparent, 
	legend=:outertopright, framestyle=:box, 
	fontfamily="Computer modern", titlefont=9, 
	markerstrokewidth=0,
	guidefont=9, title_loc=:left)

# ╔═╡ a7f8a332-9e11-4025-93ba-d25116994de2
function simest(M, nmax, Nest=100)
	i = 0
	while i < nmax
		xs = [generation(M, InfPop(z=Float64[]))]; i += 1
		while 0 < length(xs[end].z) < Nest 
			push!(xs, generation(M, xs[end])); i += 1
		end
		est = length(xs[end].z) > 0
		if est
			cyt = counts(xs[end].c, 2:4)
			return (est, i, cyt, xs)
		end
	end
	return (false, nmax, [0, 0, 0], InfPop{Float64}[])
end

# ╔═╡ 1e8cd76a-9f61-4da4-8906-9c81beb59924
let M = InfDemeMixEst(U=Umat(0.05, 0.05), γ=0.25, m=1.0, V=0.5)
	simest(M, Inf)
end

# ╔═╡ 424f97c8-7f30-4532-a212-ed353788886a
mss = 10 .^ range(log10(0.1), stop=log10(3), length=10)

# ╔═╡ 585b30e3-5f6b-4294-8399-4d3c52958a0a
Xs = map(mss) do m
	M = InfDemeMixEst(U=Umat(0.05, 0.05), γ=0.25, m=m, V=0.5, θ=fill(-1.0, 3))
	map(1:250) do _
		simest(M, Inf)
	end
end

# ╔═╡ f90e5d0b-6cb7-41fa-bc1b-b46be2ebde31
Ys = map(Xs) do xs
	pest = 1/mean(map(x->x[2], xs))
	ys = filter(x->x[1], xs)
	zs = filter(y->y[3][3] > y[3][1], ys)
	@info pest, length(xs), length(zs), length(zs)/length(ys)
	zs
end

# ╔═╡ a20c9bac-09cc-470f-8b39-8bea86b848e1
function jeffreys_interval(x, n, α)
	if x != 0 && x != n 
		quantile(Beta(x + 0.5, n - x + 0.5), [α/2, 1-α/2])
	elseif x == 0
		[0.0, quantile(Beta(x + 0.5, n - x + 0.5), 1-α/2)]
	else
		[quantile(Beta(x + 0.5, n - x + 0.5), α/2), 1.0]
	end
end

# ╔═╡ 6d12c242-a16e-4f27-b977-008f8968c4da
let tests = map(xs->mean(map(x->x[2], xs)), Xs)
	n  = length(Xs[1])
	xs = length.(Ys) 
	p4 = xs ./ n
	ci = [jeffreys_interval(x, n, 0.05) for x in xs]
	P1 = plot(mss, tests, marker=true, ms=2, legend=false, 
		color=:black, 
		xscale=:log10,
		xlabel="\$m\$", ylabel="\$\\overline{T}\$", 
		title="time to establishment")
	P2 = plot(mss, p4, 
		yerr=(p4 .- first.(ci), last.(ci) .- p4), markerstrokewidth=1,
		legend=false, 
		marker=true, ms=2, 
		xscale=:log10,
		color=:black, title="proportion of tetraploid est.",
		xlabel="\$m\$", ylabel="\$P_4\$")
	hline!([InfDemeMixEst(U=Umat(0.05, 0.05)).c[3]])
	plot(P1, P2, size=(250*2, 200), margin=3Plots.mm)
end#; savefig("../doc/img/est4.pdf")

# ╔═╡ a73a9db1-cbbd-48aa-8cba-85101c6decad
md"""
The red line on the right is the proportion of tetraploid migrants. We do not see an escape from swamping, as the migration rate goes up, the probability of tetraploid est. sinks below the proportion of tetraploid migrants, i.e. MCE is stronger than any escape from swamping effect. Or not? The increase est. probability relative to $\pi_4$ for low $m$, is it only due to lower inbreeding or does it exceed our results for no migration (note, use the right $z_s$ for comparisons...). **A quick simulation suggests the tet. est. probability to be the same for the lowest migration rate here and the single migrant case.**
"""

# ╔═╡ fc443626-0871-4419-8fd3-2e799386dc1b
let n  = length(Xs[1])
	xs = length.(Ys) 
	p4 = xs ./ n
	plot(mss, p4 ./ 0.0035, ms=2, marker=true, color=:black, ylabel="\$P_4/\\pi_4\$")
	hline!([1], legend=false)
end

# ╔═╡ 75aec322-b958-44e4-a63b-5b7946c408e4
let j=1, i=5, Yi=Ys[j][i][4]
	P1 = plot(xlabel="generation", ylabel="\$z\$")
	map(enumerate(Yi)) do (i,x)
		scatter!(P1, fill(i, length(x.z)), x.z, alpha=0.6,
			label="", legend=false, color=:black, ms=2)
	end
	C = hcat(map(x->counts(x.c, 2:4), Yi)...)'
	P2 = plot(log10.(C), marker=true, ms=2, 
		label=["diploid" "triploid" "tetraploid"], legend=:topleft,
		xlabel="generation", ylabel="\$\\log_{10} N\$")
	plot(P1, P2, size=(270*2, 220), margin=3Plots.mm)
end#; savefig("../doc/img/est3.pdf")

# ╔═╡ Cell order:
# ╠═99fd9810-18e3-11ef-3ccf-0f13790a9f6e
# ╠═65220b6c-c975-4a05-8377-19cfb88f92d8
# ╠═97808c41-877f-454b-a61a-9b3160a176e4
# ╠═a7f8a332-9e11-4025-93ba-d25116994de2
# ╠═1e8cd76a-9f61-4da4-8906-9c81beb59924
# ╠═424f97c8-7f30-4532-a212-ed353788886a
# ╠═585b30e3-5f6b-4294-8399-4d3c52958a0a
# ╠═f90e5d0b-6cb7-41fa-bc1b-b46be2ebde31
# ╠═a20c9bac-09cc-470f-8b39-8bea86b848e1
# ╠═6d12c242-a16e-4f27-b977-008f8968c4da
# ╟─a73a9db1-cbbd-48aa-8cba-85101c6decad
# ╠═fc443626-0871-4419-8fd3-2e799386dc1b
# ╠═75aec322-b958-44e4-a63b-5b7946c408e4
