using InfGenetics, Statistics, Test, StatsBase
using InfGenetics: randgen_noinbreeding

@testset "No inbreeding" begin
    d = [InfDeme(randn(100), 1., k=k) for k=[1,2,4]]
    for deme in d
        @test var(evolve(x->randgen_noinbreeding(x), deme, 100).z) ≈ 2. atol=0.2
    end
end

# verify segregation variance
function construct_tet(n, β=exp(randn()), ngen=4)
    x = randn(n, 4) ./ √(β*n*4)
    y = deepcopy(x)
    for i=1:ngen
        g1 = hcat([sample(x[k,:], 2, replace=false) for k=1:n]...)
        g2 = hcat([sample(x[k,:], 2, replace=false) for k=1:n]...)
        x = permutedims([g1 ; g2])
    end
    cs = combinations(1:4, 2)
    F = 0.
    for i=1:n
        for c in cs
            F += isequal(x[i,c]...) ? 1 : 0
        end
    end
    F /= (6n)
    return y, x, F
end

function gametes(x, d=0., n=1000)
    map(1:n) do i
        gamete(x, d)
    end
end

function gamete(x, d)
    z = 0.
    for k=1:size(x,1)
        z += rand() < d ? 2rand(x[k,:]) : sum(sample(x[k,:], 2, replace=false))
    end
    z
end

test = map(1:10) do rep
    d = 0.2
    y, x, F = construct_tet(100, 1, rand(1:4))
    v0 = gametes(y, 0.) |> var  # segvar no inbreeding, no dr
    v1 = gametes(x, 0.) |> var  # segvar with inbreeding, no dr
    v2 = gametes(y, d)  |> var  # segvar no inbreeding, with dr
    v3 = gametes(x, d)  |> var  # segvar inbreeding and dr
    # why 2?
    (F=F, v0=v0, v0d=v2, v0d_=v0*(1+2d), vf=v1, vf_=v0*(1-F), vfd=v3, vfd_=v0*(1+2d)*(1-F))
end |> DataFrame

# replace = false leads to v = v0(1+2d)
# replace = true  leads to v = v0(1+ d)
function v0(x, d=0.)
   map(1:1_000_000) do i
       rand() < d ? 2rand(x) : sum(sample(x, 2, replace=true))
   end |> var
end

x = randn(4)
v = v0(x)
@info "" v v0(x, 0.3) v*1.3


model = InfDemeMixEst(m=0.1, θ=fill(2.5, 3), U=Umat(0.05,0.05))
pop = InfPop(z=randn(20), c=fill(2, 20))
@btime InfGenetics.generation(model, pop)
