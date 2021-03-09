using InfGenetics, Statistics

@testset "No inbreeding" begin
    d = [InfDeme(randn(100), 1., k=k) for k=[1,2,4]]
    for deme in d
        @test var(evolve(x->randgen(x), d, 100)) ≈ 2.
    end
end
