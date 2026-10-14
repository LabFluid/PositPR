# Conversions: x2Posit, Posit2x, Posit2Posit, PositConverter, Int, floor,
# Float64, convert / promote.

using PositPR: Posit2Posit, PositConverter

"Random real values spread on a logarithmic scale over the whole range of the type."
function log_uniform_samples(ps, es, n; rng = Xoshiro(2024))
    span = 2^es * (ps - 2) + 2               # slightly beyond maxpos/minpos
    return [rand(rng, (-1, 1)) * 2.0^(span * (2rand(rng) - 1)) for _ in 1:n]
end

"Midpoints between consecutive posits and their immediate neighbours."
function midpoints(ps, es)
    sorted = sort(collect(real_patterns(ps)); by = u -> signed_pattern(u, ps))
    vals = [ref_value(u, ps, es) for u in sorted]
    xs = Float64[]
    for i in 1:length(vals)-1
        m = (vals[i] + vals[i+1]) / 2
        push!(xs, m, prevfloat(m), nextfloat(m))
    end
    return xs
end

@testset "Conversions" begin

    @testset "x2Posit: exactly representable values, Posit{$ps,$es}" for (ps, es) in CONFIGS_8
        T = Posit{ps, es}
        @test report("x2Posit(exact value)", check_all(real_patterns(ps)) do u
            x = ref_value(u, ps, es)
            got = pattern(x2Posit(x, T))
            got == u ? nothing : "x = $x: got $(bitstr(got, ps)), expected $(bitstr(u, ps))"
        end)
        # The internal representation (struct) must also match the one built by the
        # decoder, otherwise `==` (which compares field by field) fails.
        @test report("x2Posit(value) == PositDecoder(bits)", check_all(real_patterns(ps)) do u
            x2Posit(ref_value(u, ps, es), T) == from_pattern(u, ps, es) ? nothing :
                "x = $(ref_value(u, ps, es)): $(x2Posit(ref_value(u, ps, es), T)) ≠ $(from_pattern(u, ps, es))"
        end)
    end

    @testset "x2Posit: rounding at midpoints, Posit{$ps,$es}" for (ps, es) in CONFIGS_8
        T = Posit{ps, es}
        xs = midpoints(ps, es)
        res = check_all(xs) do x
            got, want = pattern(x2Posit(x, T)), ref_encode(x, ps, es)
            got == want ? nothing : "x = $x: got $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
        end
        @test report("correct rounding (RNE)", res)
        @test report("error ≤ 1 ULP", check_all(xs) do x
            ulp_distance(pattern(x2Posit(x, T)), ref_encode(x, ps, es), ps) <= 1 ? nothing : "x = $x"
        end)
    end

    @testset "x2Posit: random values over the whole range, Posit{$ps,$es}" for (ps, es) in [(8, 1), (8, 2), (16, 1), (16, 2), (32, 2)]
        T = Posit{ps, es}
        xs = log_uniform_samples(ps, es, 2000)
        @test report("correct rounding (RNE)", check_all(xs) do x
            got, want = pattern(x2Posit(x, T)), ref_encode(x, ps, es)
            got == want ? nothing : "x = $x: got $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
        end)
        @test report("error ≤ 1 ULP", check_all(xs) do x
            ulp_distance(pattern(x2Posit(x, T)), ref_encode(x, ps, es), ps) <= 1 ? nothing : "x = $x"
        end)
    end

    @testset "x2Posit: saturation and special values" for T in (Posit{8, 2}, Posit{16, 2})
        maxpos, minpos = Posit2x(typemax(T)), Posit2x(PositPR.smallest(T))
        # Posits never overflow to infinity nor underflow to zero.
        @test Posit2x(x2Posit(1e300, T)) == maxpos
        @test Posit2x(x2Posit(-1e300, T)) == -maxpos
        @test Posit2x(x2Posit(maxpos * 64, T)) == maxpos
        @test Posit2x(x2Posit(1e-300, T)) == minpos
        @test Posit2x(x2Posit(-1e-300, T)) == -minpos
        @test Posit2x(x2Posit(minpos / 64, T)) == minpos
        @test iszero(x2Posit(0.0, T))
        @test iszero(x2Posit(-0.0, T))
        @test iszero(x2Posit(0, T))
        @test isNaR(x2Posit(Inf, T))
        @test isNaR(x2Posit(-Inf, T))
        @test isNaR(x2Posit(NaN, T))
    end

    @testset "x2Posit: input types" begin
        T = Posit{16, 2}
        @test Posit2x(x2Posit(3, T)) == 3.0
        @test Posit2x(x2Posit(-7, T)) == -7.0
        @test Posit2x(x2Posit(1000, T)) == 1000.0
        @test Posit2x(x2Posit(3 // 8, T)) == 0.375
        @test Posit2x(x2Posit(Float32(1.25), T)) == 1.25
        @test Posit2x(x2Posit(big"0.125", T)) == 0.125
        @test Posit2x(x2Posit(π, T)) == ref_value(ref_encode(BigFloat(π), 16, 2), 16, 2)
        @test x2Posit(1.5, Val(16), Val(2)) == x2Posit(1.5, T)
    end

    @testset "Posit2x, Float64 and float" begin
        for (ps, es) in CONFIGS_8
            @test all(real_patterns(ps)) do u
                p = from_pattern(u, ps, es)
                v = ref_value(u, ps, es)
                Posit2x(p) == v && Float64(p) == v && float(p) == v
            end
        end
        @test isnan(Posit2x(NaR(Posit{8, 2})))           # NaR -> NaN
    end

    @testset "Posit2Posit: widening is exact" begin
        for ((ps0, es0), (ps1, es1)) in [((8, 0), (16, 0)), ((8, 1), (16, 1)), ((8, 2), (16, 2)),
                                         ((8, 1), (16, 2)), ((8, 2), (16, 3)), ((8, 2), (32, 2)),
                                         ((16, 2), (32, 2)), ((8, 2), (12, 2))]
            pats = ps0 == 8 ? collect(real_patterns(ps0)) :
                   filter(!=(nar_pattern(ps0)), sample_patterns(ps0; step = 13))
            @test report("Posit{$ps0,$es0} → Posit{$ps1,$es1}", check_all(pats) do u
                p = from_pattern(u, ps0, es0)
                q = Posit2Posit(p, ps1, es1)
                v = ref_value(u, ps0, es0)
                typeof(q) == PositType(ps1, es1) && Posit2x(q) == v && pattern(q) == ref_encode(v, ps1, es1) ?
                    nothing : "$v -> $(Posit2x(q))"
            end)
        end
    end

    @testset "Posit2Posit: narrowing rounds correctly" begin
        for ((ps0, es0), (ps1, es1)) in [((16, 2), (8, 2)), ((16, 1), (8, 1)), ((16, 2), (8, 1)),
                                         ((16, 1), (8, 2)), ((32, 2), (16, 2))]
            pats = filter(!=(nar_pattern(ps0)), sample_patterns(ps0; step = ps0 == 32 ? 999_983 : 13))
            @test report("RNE Posit{$ps0,$es0} → Posit{$ps1,$es1}", check_all(pats) do u
                v = ref_value(u, ps0, es0)
                got = pattern(Posit2Posit(from_pattern(u, ps0, es0), ps1, es1))
                want = ref_encode(v, ps1, es1)
                got == want ? nothing : "$v: got $(ref_value(got, ps1, es1)), expected $(ref_value(want, ps1, es1))"
            end)
            @test report("≤ 1 ULP Posit{$ps0,$es0} → Posit{$ps1,$es1}", check_all(pats) do u
                v = ref_value(u, ps0, es0)
                got = pattern(Posit2Posit(from_pattern(u, ps0, es0), ps1, es1))
                ulp_distance(got, ref_encode(v, ps1, es1), ps1) <= 1 ? nothing : "$v"
            end)
        end
        @test iszero(Posit2Posit(zero(Posit{16, 2}), 8, 2))
        @test isNaR(Posit2Posit(NaR(Posit{16, 2}), 8, 2))
    end

    @testset "PositConverter converts to the wider format" begin
        a, b = Posit{8, 2}(1.5), Posit{16, 2}(2.25)
        a2, b2 = PositConverter(a, b)
        @test typeof(a2) == typeof(b2) == PositType(16, 2)
        @test Posit2x(a2) == 1.5 && Posit2x(b2) == 2.25
        c2, d2 = PositConverter(b, a)
        @test typeof(c2) == typeof(d2) == PositType(16, 2)
    end

    @testset "Int, Integer and convert" begin
        T = Posit{16, 2}
        for n in -200:200
            @test Int(T(n)) == n
        end
        @test Integer(T(42)) == 42
        @test convert(Int, T(-17)) == -17
        @test Int(zero(T)) == 0
        # For non-integer values the library rounds to the nearest integer.
        for x in (0.25, 0.75, 1.25, 2.75, -0.75, -2.25, 100.5 + 0.25)
            @test abs(Int(T(x)) - Posit2x(T(x))) <= 0.5
        end
    end

    @testset "floor" begin
        for (ps, es) in [(8, 1), (8, 2)]
            @test report("floor Posit{$ps,$es}", check_all(real_patterns(ps)) do u
                v = ref_value(u, ps, es)
                abs(v) > 1e6 && return nothing
                got = Posit2x(floor(from_pattern(u, ps, es)))
                got == floor(v) ? nothing : "floor($v): got $got"
            end)
        end
    end

    @testset "convert and promotion" begin
        T = Posit{16, 2}
        @test convert(T, 2.5) == T(2.5)
        @test convert(T, 3) == T(3)
        @test convert(T, true) == one(T)
        @test convert(T, T(1.5)) == T(1.5)
        @test promote_type(PositType(T), Float64) == PositType(T)
        @test promote_type(Int, PositType(T)) == PositType(T)
        @test promote_type(PositType(T), Bool) == PositType(T)
        @test typeof(T(1) + 1.5) == PositType(T)
        @test Posit2x(T(1) + 1.5) == 2.5
        @test Posit2x(2.0 * T(1.5)) == 3.0
        @test T(1.5) == 1.5
        @test Base.decompose(T(1.5)) |> ((n, p, d),) -> n * 2.0^p / d == 1.5
        @test Base.decompose(T(-0.375)) |> ((n, p, d),) -> n * 2.0^p / d == -0.375
    end
end
