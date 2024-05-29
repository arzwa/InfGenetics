
function cytotype_equilibrium(U)
    u21, u22 = U[1,:]
    u31, u32 = U[2,:]
      _, u42 = U[3,:]
    g1 = (u21 - 4*u31 - 2*u32 + 2*u42)/(
        2*(u21 + u22 - 2*u31 - 2*u32 + u42)) - sqrt(
            u21^2 - 4*u21*u32 + 8*u22*u31 - 4*u22*u42 + 4*u32^2)/(
                2*(-u21 - u22 + 2*u31 + 2*u32 - u42))
    [g1^2, 2g1*(1-g1), (1-g1)^2] 
end

function rvmatrix(U)
    p = cytotype_equilibrium(U)
    Z1 = sum(U[:,1] .* p)
    Z2 = sum(U[:,2] .* p)
    M = [
        p[1]*U[1,1]/Z1                          p[2]*U[2,1]/Z1                             0.0              ;
        p[1]*U[1,1]/(3Z1) + 2p[1]*U[1,2]/(3Z2)  p[2]*U[2,1]/(3Z1) + 2p[2]*U[2,2]/(3Z2)   2*p[3]*U[3,2]/(3Z2);
        p[1]*U[1,2]/Z2                          p[2]*U[2,2]/Z2                             p[3]*U[3,2]/Z2
    ]
    return M
end

function rvs(U)
    M = rvmatrix(U)
    _, evecs = eigen(M')
    v = evecs[:,3]
    v ./ sum(v)
end

function Nerv(U, N)
    alpha = rvs(U)
    pis = cytotype_equilibrium(U)
    Nis = N .* pis
    1 / sum((alpha .^ 2) ./ Nis)
end

