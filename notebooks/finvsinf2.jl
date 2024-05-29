### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 2e34eb4f-d93f-477a-9d6a-a8844dde46f5
using Pkg; Pkg.activate("..")

# ╔═╡ 85e0ef20-c60c-45f8-82b4-eea93e03a7a6
using InfGenetics, StatsBase, Distributions, Parameters, Plots

# ╔═╡ 2b4feedf-f7ff-4b15-a665-2dfc98a1030e
default(size=(280,220), grid=false, fg_legend=:transparent, 
	legend=:outertopright, framestyle=:box, 
	fontfamily="Computer modern", titlefont=9, guidefont=9, title_loc=:left)

# ╔═╡ ec095e7e-e18b-43bb-b510-eb6b9ccc9be0
function smooth(xs::Vector{T}, window, step) where T
    n = 1
    ys = Tuple{Float64,T}[]
    while n+window < length(xs)
        x = n+window/2
        push!(ys, (x, mean(filter(!isnan, xs[n:n+window]))))
        n += step
    end
    xl = (n + length(xs))/2 
    push!(ys, (xl, mean(xs[n:end])))
end

# ╔═╡ d521aceb-0edf-4227-8008-b0695e03423a
function cb(x)  # callback function
    _, pop = x
    z2 = pop.z[pop.c .== 2]
    z3 = pop.z[pop.c .== 3]
    z4 = pop.z[pop.c .== 4]
    z2, z3, z4
end

# ╔═╡ 73ce667b-eb15-4d6f-95a3-6ed70bd6e644
X1 = let N = 250, Ls = [500, 100, 50, 20, 10], ngen=2500
    M = InfDemeMix(U=[1.0 0.0; 0.0 0.0; 0.0 1.0])
    p1 = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(4,N))
    xs1 = fiterate((x,i)->generation(x, M, z->1), deepcopy(p1), ngen, callback=cb)
    xss = map(Ls) do L
         @info L
         Gs = [randn(L,4)/√(2L) for i=1:N]
         p2 = FinPop(sum.(Gs), Gs, fill(4, N))
         fiterate((x,i)->generation(x, M, z->1), deepcopy(p2), ngen, callback=cb)
    end
    xs1, xss
end

# ╔═╡ 0185715e-6d73-445e-9ae6-b6978f642e33
P1 = let Ls = [500, 100, 50, 20, 10], N=250, w=20, s=10
	plot(0:1:2500, t->exp(-t/(4*N)), color=:black, lw=2, alpha=0.8,
		xlabel="generation", ylabel="\$V_z\$",
		label="", size=(280,210), legend=:topright)
	ys = map(xs->getindex.(xs, 3), X1[2][1:end-1])
    for (k,L,y) in reverse(collect(zip(2:length(ys)+1,Ls,ys)))
		plot!(smooth(map(var, y),w,s), lw=1.5, alpha=0.8,
			color=k, label="\$L=$L\$")
    end
	plot!(smooth(map(var, last.(X1[1])), w, s), 
		ylim=(0,Inf), color=1, lw=1.5, alpha=0.8,
		label="\$L \\rightarrow \\infty\$", title="(A) \$\\alpha_4=0\$")
end #; savefig("../doc/img/tetfininf.pdf")

# ╔═╡ 705f985e-8d5c-4188-b90c-b324a40bcc4d
X2 = let N = 250, Ls = [500, 100, 50, 20, 10], ngen=2500
    M = InfDemeMix(U=[1.0 0.0; 0.0 0.0; 0.0 1.0], α=fill(1/6, 3))
    p1 = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(4,N))
    xs1 = fiterate((x,i)->generation(x, M, z->1), deepcopy(p1), ngen, callback=cb)
    xss = map(Ls) do L
         @info L
         Gs = [randn(L,4)/√(2L) for i=1:N]
         p2 = FinPop(sum.(Gs), Gs, fill(4, N))
         fiterate((x,i)->generation(x, M, z->1), deepcopy(p2), ngen, callback=cb)
    end
    xs1, xss
end

# ╔═╡ bc608733-946c-4768-bd3d-d63fde286955
P2 = let X=X2, Ls = [500, 100, 50, 20, 10], N=250, w=20, s=10
	plot(0:1:2500, t->exp(-t/(4*N)), color=:black, lw=2, alpha=0.8,
		xlabel="generation", ylabel="\$V_z\$",
		label="", size=(280,210), legend=:topright)
	ys = map(xs->getindex.(xs, 3), X[2][1:end-1])
    for (k,L,y) in reverse(collect(zip(2:length(ys)+1,Ls,ys)))
		plot!(smooth(map(var, y),w,s), lw=1.5, alpha=0.8,
			color=k, label="\$L=$L\$")
    end
	plot!(smooth(map(var, last.(X[1])), w, s), 
		ylim=(0,1+2/6), color=1, lw=1.5, alpha=0.8,
		label="\$L \\rightarrow \\infty\$",
		title="(B) \$\\alpha_4=1/6\$",
		)
end 

# ╔═╡ 8d712bbb-156c-460c-b1c6-4099f4922403
plot(P1, P2, size=(2*270, 210), margin=3Plots.mm, ylim=(0,1+1.5/6)) #; savefig("../doc/img/tetfininf.pdf")

# ╔═╡ Cell order:
# ╠═2e34eb4f-d93f-477a-9d6a-a8844dde46f5
# ╠═85e0ef20-c60c-45f8-82b4-eea93e03a7a6
# ╠═2b4feedf-f7ff-4b15-a665-2dfc98a1030e
# ╠═ec095e7e-e18b-43bb-b510-eb6b9ccc9be0
# ╠═d521aceb-0edf-4227-8008-b0695e03423a
# ╠═73ce667b-eb15-4d6f-95a3-6ed70bd6e644
# ╠═0185715e-6d73-445e-9ae6-b6978f642e33
# ╠═705f985e-8d5c-4188-b90c-b324a40bcc4d
# ╠═bc608733-946c-4768-bd3d-d63fde286955
# ╠═8d712bbb-156c-460c-b1c6-4099f4922403
