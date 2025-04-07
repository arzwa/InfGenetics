
function condition(n, pop, M, nmax, Nest) 
    N = popsize(pop)
    (M.m == 0.0 && N == 0)   || (n >  nmax) || (N >= Nest) 
    # extinct & no migration || too long    || established
end

# Simulate until establishment
function simest(M; x0=InfPop(z=Float64[]), nmax=Inf, Nest=100)
    xs = [(1, generation(M, x0))]
    while !condition(last(xs)..., M, nmax, Nest)
        (n, x) = last(xs) 
        x_ = generation(M, x, Nmax=Nest) # no need to simulate > Nest
        if popsize(x_) == 0 
            # extinction -- prune the stored generations
            xs = [(n + 1, x_)]
        else 
            push!(xs, (n + 1, x_))
        end
	end
    return (popsize(last(xs)[2]) > 0, xs)
end

# don't store the pop's
function simest2(M; x0=InfPop(z=Float64[]), nmax=Inf, Nest=100)
    x = generation(M, x0); n = 1
    while !condition(n, x, M, nmax, Nest)
        x = generation(M, x, Nmax=Nest); n += 1
	end
    N = popsize(x)
    c = N > 0 ? counts(x.c, 2:4) : [0, 0, 0]
    return (N > 0, n, c)
end

# weak migration
function _simest3(M; x=InfPop(z=Float64[]), nmax=Inf, Nest=100)
    p0 = pdf(Poisson(M.m), 0)
    M0 = reconstruct(M, m=0)
    n = 0
    while !condition(n, x, M, nmax, Nest)
        if popsize(x) == 0
            n += rand(Geometric(1-p0))
            x = add_n_migrants(M, 1, x)  
            # XXX this should be the distribution conditional on being nonzero,
            # for small m this is very likely 1
        end
        x = generation(M0, x, Nmax=Nest); n+= 1
	end
    N = popsize(x)
    c = N > 0 ? counts(x.c, 2:4) : [0, 0, 0]
    return (N > 0, n, c)
end

function simuntil(f) 
	while true
		x, xs = f()
		x && return xs
	end
end

