### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ bc36ad7c-18ca-11ef-286f-a523cc8cf5e5
using Pkg; Pkg.activate("..")

# ╔═╡ 18369b1e-76e0-4ba3-9e6d-4fb2a4099307
using InfGenetics, StatsBase, Distributions, Parameters, Plots, SymPy

# ╔═╡ 2c9b1135-5bd0-46bf-ba0b-b2629d4ab511
default(
	size=(280,220), 
	grid=false, 
	fg_legend=:transparent, 
	legend=:outertopright, 
	framestyle=:box, 
	fontfamily="Computer modern", 
	titlefont=9, 
	guidefont=9, 
	title_loc=:left
)

# ╔═╡ d22a3791-2965-4e29-bdcf-6a17b8e6993e
u, v, g, w2, w3, w4 = @syms u v g w2 w3 w4

# ╔═╡ 66822d67-3ca4-4427-8863-0fc3b430b2f9
function nextgen(g1, u, v, w2, w3, w4) 
	g2 = 1 - g1 
	x2 = g1^2
	x3 = 2g1*g2
	x4 = g2^2
	w̄  = x2*w2 + x3*w3 + x4*w4
	x2_ = x2*w2/w̄
	x3_ = x3*w3/w̄
	x4_ = x4*w4/w̄
	g1_ = x2_*(1-u) + x3_*v
	g2_ = x2_*u + x3_*v + x4_*(1-u)
	g1__ = g1_/(g1_ + g2_)	
	g2__ = g2_/(g1_ + g2_)
	g1__
end

# ╔═╡ e4478703-873d-4394-8159-46a42516db22
sols = solve(simplify((nextgen(g, u, v, w2, w3, w4) - g)/g), g)

# ╔═╡ 828aa9c7-038d-4926-a211-83074e4c723a
ex = simplify(subs(sols[2], w3=>w4, w2=>1))

# ╔═╡ ac0f69fb-c68c-48b9-a5ae-f2c75e148b4b
ex2 = simplify(subs((nextgen(g, u, v, w2, w3, w4) - g), w3=>w4, w2=>1))

# ╔═╡ 1eb04a17-5530-4be4-b0fc-568767ffe58d
fun = lambdify(ex2, (g, u, v, w4))

# ╔═╡ ecab70d9-81d6-4ef7-8e3a-67f28cdbfa7c
plot(); map(1:0.3:3) do w4 
	plot!(g->-fun(1-g, 0.05, 0.05, w4), label="\$w=$w4\$")
end; hline!([0], label="")

# ╔═╡ Cell order:
# ╠═bc36ad7c-18ca-11ef-286f-a523cc8cf5e5
# ╠═18369b1e-76e0-4ba3-9e6d-4fb2a4099307
# ╠═2c9b1135-5bd0-46bf-ba0b-b2629d4ab511
# ╠═d22a3791-2965-4e29-bdcf-6a17b8e6993e
# ╠═66822d67-3ca4-4427-8863-0fc3b430b2f9
# ╠═e4478703-873d-4394-8159-46a42516db22
# ╠═828aa9c7-038d-4926-a211-83074e4c723a
# ╠═ac0f69fb-c68c-48b9-a5ae-f2c75e148b4b
# ╠═1eb04a17-5530-4be4-b0fc-568767ffe58d
# ╠═ecab70d9-81d6-4ef7-8e3a-67f28cdbfa7c
