using InfGenetics, Distributions, StatsBase
using Plots, PlotThemes; theme(:hokusai)

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

# look at the trait
function cb(x)
    _, pop = x
    z2 = pop.z[pop.c .== 2]
    z3 = pop.z[pop.c .== 3]
    z4 = pop.z[pop.c .== 4]
    z2, z3, z4
end

N = 500
Ls = [500, 100, 50, 20, 10]
xs = [(α=zeros(3), β=ones(3)), 
      (α=zeros(3), β=[1, √(2/3), √(1/2)]), 
      (α=[1/2, 1/4, 1/6], β=[1, √(2/3), √(1/2)])]
res = map(xs) do (α, β)
    @info α, β
    M = InfDemeMix(α=α, β=β, U=[0.92 0.08; 0.08 0.08; 0.0 0.92])
    p1 = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
    xs1 = fiterate((x,i)->generation(x, M, z->1), deepcopy(p1), 1000, callback=cb)
    xss = map(Ls) do L
        @info L
        Gs = [randn(L,2)/√(2L) for i=1:N]
        p2 = FinPop(sum.(Gs), Gs, fill(2, N))
        fiterate((x,i)->generation(x, M, z->1), deepcopy(p2), 1000, callback=cb)
    end
    xs1, xss
end

map(zip(res, xs)) do ((xs1, xss), (a, b))
    b = round.(b .^ 2, digits=2)
    a = round.(a, digits=2)
    map(1:3) do i
        tt = ["diploid, \$\\beta_2^2=$(b[1]), \\alpha_2=$(a[1])\$", 
              "triploid, \$\\beta_3^2=$(b[2]), \\alpha_3=$(a[2])\$", 
              "tetraploid, \$\\beta_4^2=$(b[3]), \\alpha_4=$(a[3])\$"]
        z = getindex.(xs1, i)
        ys = map(xs->getindex.(xs, i), xss[1:end-1])
        plot(smooth(map(var, z), 20, 10), title=tt[i], label="\$L \\rightarrow \\infty\$")
        plot!(map(var, z), color=1, label="", alpha=0.1)
        for (k,L,y) in zip(2:length(ys)+1,Ls,ys)
            plot!(smooth(map(var, y), 20, 10), color=k, label="\$L=$L\$", title=tt[i])
            plot!(map(var, y), color=k, label="", alpha=0.1)
        end
        plot!(ylim=(0,[1.5, 2.5, 4.2][i]))
    end |> x->plot(x..., size=(750,230), layout=(1,3), legend=:topright,
        xlabel="generation", ylabel="\$V_z\$", bottom_margin=1Plots.mm,
        margin=1Plots.mm, titlefont=9)
end |> x->plot(x..., layout=(3,1), size=(700,650))


ps = [plot() for i=1:3]
map(enumerate(zip(res, xs))) do (j,((xs1, xs2), (a, b)))
    b = map(x->@sprintf("%.2f", x^2), b)
    a = round.(a, digits=2)
    map(1:3) do i
        tt = ["\$\\beta^2=$(b[1]), \\alpha=$(a[1])\$", 
              "\$\\beta^2=$(b[2]), \\alpha=$(a[2])\$", 
              "\$\\beta^2=$(b[3]), \\alpha=$(a[3])\$"]
        tp = ["diploid", "triploid", "tetraploid"]
        z = getindex.(xs1, i)
        plot!(ps[i], smooth(map(var, z), 20, 10), 
            label=tt[i], title="$(tp[i])", alpha=1.7, legend=:topright)
    end 
end
plot(ps..., layout=(1,3), size=(700,220), xlabel="generation",
    ylabel="\$V_z\$", bottom_margin=4Plots.mm, margin=2Plots.mm)


# look at inbreeding
function cb2(x)
    _, pop = x
    z2 = pop.F[pop.c .== 2]
    z3 = pop.F[pop.c .== 3]
    z4 = pop.F[pop.c .== 4]
    z2, z3, z4
end

xsb = [
    (α=zeros(3), β=[1, √(2/3), √(1/2)]), 
    (α=fill(1/6, 3), β=[1, √(2/3), √(1/2)]),
    (α=[1/2, 1/4, 1/6], β=[1, √(2/3), √(1/2)]),
]
res = map(xsb) do (α, β)
    @info α, β
    M = InfDemeMix(α=α, β=β, U=[0.92 0.08; 0.08 0.08; 0.0 0.92])
    p1 = InfPop(z=rand(Normal(0, √(2M.V)), N), c=fill(2,N))
    xs1 = fiterate((x,i)->generation(x, M, z->1), deepcopy(p1), 2000, callback=cb2)
end

ps = [plot() for i=1:3]
map(enumerate(zip(res, xsb))) do (j,(xs1, (a, b)))
    b = map(x->@sprintf("%.2f", x^2), b)
    a = map(x->@sprintf("%.2f", x  ), a)
    map(1:3) do i
        tt = ["\$\\alpha=$(a[1])\$", 
              "\$\\alpha=$(a[2])\$", 
              "\$\\alpha=$(a[3])\$"]
        tp = ["diploid", "triploid", "tetraploid"]
        F = getindex.(xs1, i)
        plot!(ps[i], map(mean, F), 
            label=tt[i], title="$(tp[i])", alpha=1.7, legend=:bottomright)
    end 
end
plot(ps..., layout=(1,3), size=(650,200), xlabel="generation", ylim=(0,1),
    ylabel="\$\\overline{F}\$", bottom_margin=4Plots.mm, margin=2Plots.mm)

# Unreduced gametes
us = [0.01, 0.025, 0.05, 0.1]
N = 200
res = map(us) do u
    U = [1-u u; u u; 0.0 1-u]
    NN = ceil(Int, N/(1-2u))
    @info u, (1-2u), Nerv(U, N)/N
    ne = Nerv(U, N)/N
    NN = ceil(Int, N/ne)
#    @info u, NN, (1-2u)
# Ne = (1-2u)N seems to work better than this eigenvalue thing?
    M = InfDemeMix(U=U, α=zeros(3))
    p1 = InfPop(z=rand(Normal(0, √(2M.V)), NN), c=fill(2,NN))
    xs1 = fiterate((x,i)->generation(x, M, z->1), deepcopy(p1), 2000, callback=cb2)
end;

P1 = plot()
P2 = plot()
map(zip(us, res)) do (u, x)
    F2 = getindex.(x, 1)
    F  = map(y->vcat(y...), x)
    U = [1-u u; u u; 0.0 1-u]
    #NN = ceil(Int, N/(1-2u))
    ne = Nerv(U, N)/N
    NN = ceil(Int, N/ne)
    us = @sprintf "%.3f" u
    plot!(P1, map(mean, F), title="mixed-ploidy population", 
        label="\$N=$NN, u=$us\$",alpha=1.7, legend=:bottomright)
    plot!(P2, map(mean, F2), title="diploids", label="\$N=$NN, u=$us\$",
        alpha=1.7, legend=:bottomright)
end
plot!(P1, t->1-exp(-t/2N), color=:orange, ls=:dash, label="\$1-e^{-t/2N_e}\$")
plot!(P2, t->1-exp(-t/2N), color=:orange, ls=:dash, label="\$1-e^{-t/2N_e}\$")
plot!(P1, P2, size=(510,220), xlabel="generation", ylim=(0,1),
    ylabel="\$\\overline{F}\$", bottom_margin=4Plots.mm, margin=2Plots.mm)

InfGenetics.mpeq(0.05, 0.05)

