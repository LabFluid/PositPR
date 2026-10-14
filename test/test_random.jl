# Random number generation: rand (uniform) and randn (normal, Box-Muller).

using Statistics: mean, std

@testset "Random numbers" begin

    @testset "rand – Posit{$ps,$es}" for (ps, es) in [(8, 1), (8, 2), (16, 2)]
        T = Posit{ps, es}
        @test rand(T) isa T
        v = rand(T, 3, 4)
        @test size(v) == (3, 4) && eltype(v) == PositType(T)
        sample = Posit2x.(rand(T, 4000))
        @test all(x -> 0 <= x <= 1, sample)
        @test all(x -> 0 <= x < 1, sample)            # half-open interval [0, 1)
        @test abs(mean(sample) - 0.5) < 0.03
        @test abs(std(sample) - sqrt(1 / 12)) < 0.03
        # Reproducibility with an explicit generator.
        @test rand(Xoshiro(1), T, 10) == rand(Xoshiro(1), T, 10)
    end

    # randn uses log, sqrt, cos and sin in Posit arithmetic (while loops): via gcall.
    @testset "randn – Posit{$ps,$es}" for (ps, es) in [(8, 2), (16, 2)]
        T = Posit{ps, es}
        @test gcall(randn, T) isa T
        v = gcall(randn, T, 2, 3)
        @test size(v) == (2, 3) && eltype(v) == PositType(T)
        @test length(gcall(randn, T, 5)) == 5                 # odd count
        sample = Posit2x.(gcall(randn, Xoshiro(3), T, 400; timeout = 300.0))
        @test abs(mean(sample)) < 0.2
        @test abs(std(sample) - 1) < 0.2
        @test gcall(randn, Xoshiro(9), T, 6) == gcall(randn, Xoshiro(9), T, 6)
    end
end
