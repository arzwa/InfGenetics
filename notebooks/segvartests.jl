using Random, StatsBase, Distributions

function segtetraploid(X, δs, α, n)
    @assert sum(δs) == 1
    hom = sample(1:5, Weights(δs))
    shuffle!(X)
    if hom == 2
        X[1] = X[2]
    elseif hom == 3
        X[1] = X[2]
        X[3] = X[4]
    elseif hom == 4
        X[1] = X[2] = X[3]
    elseif hom ==5
        X[1] = X[2] = X[3] = X[4]
    end
    map(1:n) do _
        if rand() < α 
            # double reduction
            2*rand(X)
        else
            X1, X2 = sample(X, 2, replace=false)
            X1 + X2
        end
    end
end

δs = [0.6, 0.2, 0.1, 0.07, 0.03]
F = δs[2]/6 + δs[3]/3 + δs[4]/2 + δs[5]
Vx = 1
α = 0.2
V = map(1:10000) do _
    # Var[Y|X]
    var(segtetraploid(rand(Normal(0,√Vx), 4), δs, α, 1000))
end |> mean  # E[Var[Y|X]]
@info V, (1+2α)*(1-F)*Vx

function segtriploid_diploid(X, δs, ξ, n)
    @assert sum(δs) == 1
    hom = sample(1:3, Weights(δs))
    shuffle!(X)
    if hom == 2
        X[1] = X[2]
    elseif hom == 3
        X[1] = X[2] = X[3]
    end
    map(1:n) do _
        if rand() < ξ
            2*rand(X)
        else
            X1, X2 = sample(X, 2, replace=false)
            X1 + X2
        end
    end
end

δs = [0.6, 0.25, 0.15]
F = δs[2]/3 + δs[3]
Vx = 1
ξ = 0.4
V = map(1:10000) do _
    # Var[Y|X]
    var(segtriploid_diploid(rand(Normal(0,√Vx), 3), δs, ξ, 1000))
end |> mean  # E[Var[Y|X]]
@info V, (2/3)*(1+3ξ)*(1-F)*Vx


function segtriploid_haploid(X, δs, n)
    @assert sum(δs) == 1
    hom = sample(1:3, Weights(δs))
    shuffle!(X)
    if hom == 2
        X[1] = X[2]
    elseif hom == 3
        X[1] = X[2] = X[3]
    end
    map(1:n) do _
        rand(X)
    end
end

δs = [0.6, 0.25, 0.15]
F = δs[2]/3 + δs[3]
Vx = 1
V = map(1:10000) do _
    # Var[Y|X]
    var(segtriploid_haploid(rand(Normal(0,√Vx), 3), δs, 5000))
end |> mean  # E[Var[Y|X]]
@info V, (2/3)*(1-F)*Vx



function segdiploid_diploid(X, F, ξ, n)
    @assert sum(δs) == 1
    if rand() < F
        X[1] = X[2]
    end
    map(1:n) do _
        rand() < ξ ? 2rand(X) : X[1]+X[2]
    end
end

F = 0.3
ξ = 0.6
Vx = 1
V = map(1:10000) do _
    # Var[Y|X]
    var(segdiploid_diploid(rand(Normal(0,√Vx), 2), F, ξ, 5000))
end |> mean  # E[Var[Y|X]]
@info V, 2ξ*(1-F)*Vx
