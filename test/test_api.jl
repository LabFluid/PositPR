# Public interface: exports, parametric type, constants and constructors.

using PositPR: NaR, isNaR, smallest

@testset "API and constants" begin

    @testset "Exports declared in PositPR.jl" begin
        # All these functions are in the module's `export` list and must be
        # available after `using PositPR`.
        exported = names(PositPR)
        for f in (:Posit, :useed, :x2Posit, :Posit2x, :Posit2Posit, :PositConverter,
                  :PositDecoder, :PositEncoder, :nextPosit, :prevPosit, :isNaR, :NaR)
            @test f in exported
        end
    end

    @testset "Parametric type Posit{ps,es}" begin
        @test Posit{8, 2} <: Real
        @test PositPR.ps(Posit{16, 1}(1.0)) == 16
        @test PositPR.es(Posit{16, 1}(1.0)) == 1
        @test typeof(Posit{8, 2}(1.0)) === Posit{8, 2, UInt8}
        @test Posit{8, 2}(1.0) isa Posit{8, 2}
        @test PositPR.useed(0) == 2
        @test PositPR.useed(1) == 4
        @test PositPR.useed(2) == 16
        @test PositPR.useed(3) == 256
        for (ps, es) in [(6, 1), (8, 2), (10, 3), (12, 0), (16, 2), (24, 2), (32, 2)]
            # Arbitrary configurations must work.
            @test Posit2x(Posit{ps, es}(1.0)) == 1.0
            @test Posit2x(Posit{ps, es}(-2.0)) == -2.0
        end
    end

    @testset "Storage type chosen from ps" begin
        # The bit pattern is stored in the smallest unsigned integer that fits.
        for (ps, U) in [(2, UInt8), (8, UInt8), (9, UInt16), (12, UInt16), (16, UInt16),
                        (17, UInt32), (32, UInt32), (33, UInt64), (64, UInt64), (65, UInt128), (128, UInt128)]
            T = PositType(ps, 2)
            @test T === Posit{ps, 2, U}
            @test isconcretetype(T) && isbitstype(T)
            @test sizeof(T) == sizeof(U)
            @test typeof(Posit{ps, 2}(1.5)) === T
        end
        @test Posit8{2} === Posit{8, 2, UInt8}
        @test Posit16{2} === Posit{16, 2, UInt16}
        @test Posit32{2} === Posit{32, 2, UInt32}
        @test Posit64{2} === Posit{64, 2, UInt64}
        @test Posit16{2}(2.5) == Posit{16, 2}(2.5)
        @test PositType(Posit{16, 2}) === Posit16{2}
        @test_throws ArgumentError PositType(129, 2)(1.0)
        # Arrays of posits use the concrete type and are compact.
        v = Posit16{2}[1, 2.5, 3]
        @test eltype(v) === Posit16{2}
        @test sizeof(v) == 3 * 2
        @test eltype(zeros(Posit{8, 2}, 4)) === Posit8{2}
        # Values hash like the equal Float64, as Julia requires for == across types.
        @test hash(Posit{8, 2}(1.5)) == hash(1.5)
    end

    @testset "Constants of $T" for T in (Posit{8, 0}, Posit{8, 2}, Posit{16, 2}, Posit{32, 2})
        ps, es = PositPR.ps(one(T)), PositPR.es(one(T))
        @test pattern(zero(T)) == 0
        @test pattern(one(T)) == UInt64(1) << (ps - 2)
        @test pattern(NaR(T)) == nar_pattern(ps)
        @test Posit2x(one(T)) == 1.0
        @test Posit2x(zero(T)) == 0.0
        @test iszero(zero(T)) && !iszero(one(T))
        @test isone(one(T)) && !isone(zero(T))
        @test isNaR(NaR(T)) && !isNaR(zero(T)) && !isNaR(one(T))
        @test one(one(T)) == one(T)
        @test zero(one(T)) == zero(T)

        # maxpos = useed^(ps-2), minpos = useed^-(ps-2)
        @test Posit2x(typemax(T)) == 2.0^(2^es * (ps - 2))
        @test pattern(typemax(T)) == nar_pattern(ps) - 1
        @test Posit2x(smallest(T)) == 2.0^(-2^es * (ps - 2))
        @test pattern(smallest(T)) == 1
        @test Posit2x(typemin(T)) == -2.0^(2^es * (ps - 2))
        @test pattern(typemin(T)) == nar_pattern(ps) + 1
    end

    @testset "Constructors and copy" begin
        T = Posit{16, 2}
        p = T(3.25)
        @test T(p) == p
        @test copy(p) == p
        @test copy(p) !== nothing
        @test Posit2x(T(3)) == 3.0
        @test Posit2x(T(3 // 4)) == 0.75
        @test Posit2x(T(Float32(0.5))) == 0.5
        @test Posit2x(T(big"2.5")) == 2.5
        @test T(true) == one(T)
        @test T(false) == zero(T)
    end

    @testset "zeros and ones" begin
        T = Posit{8, 2}
        z = zeros(T, 2, 3)
        o = ones(T, 4)
        @test size(z) == (2, 3) && eltype(z) == PositType(T)
        @test all(iszero, z)
        @test size(o) == (4,) && eltype(o) == PositType(T)
        @test all(isone, o)
    end
end
