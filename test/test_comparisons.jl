# Comparisons, abs, eps, nextPosit / prevPosit.

using PositPR: nextPosit, prevPosit

@testset "Comparisons and neighbours" begin

    @testset "==, <, >, <=, >=, isless on all pairs, Posit{$ps,$es}" for (ps, es) in [(8, 0), (8, 2)]
        pats = collect(real_patterns(ps))
        posits = Dict(u => from_pattern(u, ps, es) for u in pats)
        pairs = [(a, b) for a in pats for b in pats]
        for (name, op) in [("==", ==), ("<", <), (">", >), ("<=", <=), (">=", >=), ("isless", isless)]
            @test report(name, check_all(pairs) do (a, b)
                got = op(posits[a], posits[b])
                want = op(ref_value(a, ps, es), ref_value(b, ps, es))
                got == want ? nothing : "$(ref_value(a, ps, es)) $name $(ref_value(b, ps, es)): got $got"
            end)
        end
    end

    @testset "NaR" begin
        T = Posit{8, 2}
        nar = NaR(T)
        @test nar == nar                         # Posit standard: NaR == NaR
        @test isNaR(nar)
        for x in (typemin(T), -one(T), zero(T), one(T), typemax(T))
            @test nar < x                        # NaR is smaller than every real
            @test !(nar > x)
            @test nar != x
        end
    end

    @testset "Comparison with other Real types" begin
        T = Posit{16, 2}
        @test T(1) < 2.0
        @test T(3) > 2
        @test T(2) <= 2.0
        @test T(2) >= 2
        @test 0.5 < T(1)
        @test 3 > T(2.5)
        @test 1.0 <= T(1)
        @test 4 >= T(3.5)
        @test isless(T(1), 1.5) && isless(1.0, T(1.5))
        @test T(0.5) >= Posit{8, 2}(0.25)        # different Posit types
        @test sort([T(3), T(-1), T(0.5), T(-2.25), zero(T)]) ==
              [T(-2.25), T(-1), zero(T), T(0.5), T(3)]
        @test maximum([T(1), T(7), T(-9)]) == T(7)
        @test minimum([T(1), T(7), T(-9)]) == T(-9)
    end

    @testset "abs and negation, Posit{$ps,$es} (exhaustive)" for (ps, es) in CONFIGS_8
        @test report("abs", check_all(real_patterns(ps)) do u
            got = Posit2x(abs(from_pattern(u, ps, es)))
            got == abs(ref_value(u, ps, es)) ? nothing : "abs($(ref_value(u, ps, es))) = $got"
        end)
        @test report("-x", check_all(real_patterns(ps)) do u
            got = pattern(-from_pattern(u, ps, es))
            got == (~u + 1) & mask(ps) ? nothing : "-($(ref_value(u, ps, es))) = $(ref_value(got, ps, es))"
        end)
    end

    @testset "nextPosit / prevPosit walk along the posit line, Posit{$ps,$es}" for (ps, es) in CONFIGS_8
        # Every pattern from -maxpos to maxpos, except the end points.
        inner = [u for u in real_patterns(ps) if u != nar_pattern(ps) - 1 && u != nar_pattern(ps) + 1]
        @test report("nextPosit", check_all(inner) do u
            u == nar_pattern(ps) - 1 && return nothing
            got = pattern(gcall(nextPosit, from_pattern(u, ps, es)))
            want = (u + 1) & mask(ps)
            got == want ? nothing : "next($(ref_value(u, ps, es))) = $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
        end)
        @test report("prevPosit", check_all(inner) do u
            got = pattern(gcall(prevPosit, from_pattern(u, ps, es)))
            want = (u - 1) & mask(ps)
            got == want ? nothing : "prev($(ref_value(u, ps, es))) = $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
        end)
    end

    @testset "eps" begin
        for (ps, es) in [(8, 0), (8, 2), (16, 1), (16, 2)]
            T = Posit{ps, es}
            fs = ps - 3 - es                     # fraction bits in the binade [1,2)
            @test Posit2x(gcall(eps, T)) == 2.0^-fs
            @test Posit2x(gcall(eps, one(T))) == 2.0^-fs
            @test gcall(eps, zero(T)) == PositPR.smallest(T)
        end
        # eps(x) = distance to the next posit of larger magnitude (rounded).
        for (ps, es) in [(8, 1), (8, 2)]
            pos = [u for u in real_patterns(ps) if 0 < u < nar_pattern(ps) - 1]
            @test report("eps(x) Posit{$ps,$es}", check_all(pos) do u
                got = Posit2x(gcall(eps, from_pattern(u, ps, es)))
                # the spacing is not always a posit: we expect the rounded value
                spacing = ref_value(u + 1, ps, es) - ref_value(u, ps, es)
                want = ref_value(ref_encode(spacing, ps, es), ps, es)
                got == want ? nothing : "eps($(ref_value(u, ps, es))) = $got, expected $want"
            end)
        end
    end
end
