
M = InfDemeMixEst(m=1.0, γ=0.0, σ=ones(3), θ=fill(0.0, 3), U=Umat(0.0,0.0))
x = InfPop(z=zeros(4))
W, Z, V = InfGenetics.families(M, x)
InfGenetics.selfing!(W, fill(Inf, 3), x.c)
sum(W)

