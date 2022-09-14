module InfGenetics

using LinearAlgebra, Distributions, StatsBase, Random, Parameters

export fiterate

include("infdememix.jl")

# iterate a function (monoidal)
function fiterate(fun, x, n; callback=identity)
    map(1:n) do i
        x = fun(x, i)
        callback((i, x))
    end
end

end # module
