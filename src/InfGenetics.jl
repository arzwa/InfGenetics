module InfGenetics

using LinearAlgebra, Distributions, StatsBase, Random, Parameters

export jeffreys_interval, Umat, simest, simest2, simuntil, fiterate, maxα, strs
export InfDemeMixEstPair, InfDemeMixEstSelf, InfDemeMixEst, InfDemeMix, InfDemeHS, InfPop, FinPop, generation

include("infmix.jl")
include("finmix.jl")
include("infmixest.jl")
include("simulation.jl")
include("rv.jl")
#include("infmixhs.jl")
#include("infmixestpair.jl")
#include("infmixestself.jl")

# iterate a function
function fiterate(fun, x, n; callback=identity)
    map(1:n) do i
        x = fun(x, i)
        callback((i, x))
    end
end

function jeffreys_interval(x, n, α)
	if x != 0 && x != n 
		quantile(Beta(x + 0.5, n - x + 0.5), [α/2, 1-α/2])
	elseif x == 0
		[0.0, quantile(Beta(x + 0.5, n - x + 0.5), 1-α/2)]
	else
		[quantile(Beta(x + 0.5, n - x + 0.5), α/2), 1.0]
	end
end

mpeq(u, v) = (3 - 3u - 6v + √((u + 2v - 1)*(5u + 2v - 1)))/(4 - 2u -8v)
Umat(u, v) = [1-u u ; v v ; 0.0 1-u]
const maxα = [1/2, 1/4, 1/6]
const strs = ["", "diploid", "triploid", "tetraploid"]

end # module
