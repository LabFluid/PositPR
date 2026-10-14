# Quire: exact accumulation, fused dot product, fused sum and fused multiply-add.
#
# The reference is the exact result computed in BigFloat (1024 bits, exact for
# these sums) and rounded once by the oracle's ref_encode.

"Exact Σ xᵢ·yᵢ (BigFloat) from bit patterns."
exact_dot(xs, ys, ps, es) = sum(big_value(x, ps, es) * big_value(y, ps, es) for (x, y) in zip(xs, ys); init = big(0.0))

"Random real (non-NaR) patterns."
function random_real_patterns(ps, n, rng)
    out = UInt64[]
    while length(out) < n
        u = rand(rng, UInt64) & mask(ps)
        u != nar_pattern(ps) && push!(out, u)
    end
    return out
end

@testset "Quire" begin

    @testset "Quire(p) holds p exactly, Posit{$ps,$es} (exhaustive)" for (ps, es) in CONFIGS_8
        @test report("Posit(Quire(p)) == p", check_all(real_patterns(ps)) do u
            p = from_pattern(u, ps, es)
            Posit(Quire(p)) == p ? nothing : "p = $(ref_value(u, ps, es))"
        end)
    end

    @testset "fma(a, b, c) rounds once, Posit{$ps,$es} (all pairs a, b)" for (ps, es) in [(8, 0), (8, 2)]
        pats = collect(real_patterns(ps))
        cs = UInt64[0, 0x41, 0xb3, mask(ps) >> 1, 1]
        for c in cs
            @test report("fma(a, b, $(ref_value(c, ps, es)))", check_all(((a, b) for a in pats for b in pats)) do (a, b)
                got = pattern(fma(from_pattern(a, ps, es), from_pattern(b, ps, es), from_pattern(c, ps, es)))
                exact = big_value(a, ps, es) * big_value(b, ps, es) + big_value(c, ps, es)
                want = ref_encode(exact, ps, es)
                got == want ? nothing :
                    "fma($(ref_value(a, ps, es)), $(ref_value(b, ps, es)), $(ref_value(c, ps, es))): got $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
            end)
        end
    end

    @testset "fdp and fsum round once, Posit{$ps,$es} (random vectors)" for (ps, es) in [(8, 2), (16, 1), (16, 2), (32, 2)]
        rng = Xoshiro(ps * 100 + es)
        cases = [(random_real_patterns(ps, n, rng), random_real_patterns(ps, n, rng)) for n in (1, 2, 5, 50, 500) for _ in 1:40]
        @test report("fdp(xs, ys)", check_all(cases) do (xs, ys)
            got = pattern(fdp([from_pattern(u, ps, es) for u in xs], [from_pattern(u, ps, es) for u in ys]))
            want = ref_encode(exact_dot(xs, ys, ps, es), ps, es)
            got == want ? nothing : "n = $(length(xs)): got $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
        end)
        @test report("fsum(xs)", check_all(cases) do (xs, _)
            got = pattern(fsum([from_pattern(u, ps, es) for u in xs]))
            want = ref_encode(sum(big_value(u, ps, es) for u in xs; init = big(0.0)), ps, es)
            got == want ? nothing : "n = $(length(xs)): got $(ref_value(got, ps, es)), expected $(ref_value(want, ps, es))"
        end)
    end

    @testset "Accumulation API" begin
        T = Posit{16, 2}
        q = Quire(T)
        @test iszero(q) && !isNaR(q)
        @test iszero(Posit(q))
        qmuladd!(q, T(3), T(4))                    # 12
        qadd!(q, T(0.5))                           # 12.5
        qsub!(q, T(2))                             # 10.5
        qmulsub!(q, T(1.5), T(3))                  # 6
        @test Posit(q) == T(6)
        @test Posit{16, 2}(q) == T(6)
        @test Posit16{2}(q) == T(6)
        @test BigFloat(q) == 6
        qneg!(q)
        @test Posit(q) == T(-6)
        qabs!(q)
        @test Posit(q) == T(6)
        r = Quire(T(1))
        @test Posit(q + r) == T(7)
        @test Posit(q - r) == T(5)
        @test Posit(-q) == T(-6)
        @test Posit(q) == T(6)                     # + and - do not modify q
        qadd!(q, r)
        @test Posit(q) == T(7)
        qclear!(q)
        @test iszero(q)
        @test sprint(show, Quire(T(2.5))) == "Quire{16,2}(2.5)"
    end

    @testset "Exactness where ordinary arithmetic loses information" begin
        T = Posit{16, 2}
        # Catastrophic cancellation: (2^20 + 1) - 2^20
        xs = T[2^20, 1, -2^20]
        @test iszero(sum(xs))                      # ordinary sum loses the 1
        @test fsum(xs) == one(T)                   # the quire keeps it
        # Many small terms: Σ 0.001 over 1000 terms
        a = fill(T(0.001), 1000)
        b = fill(one(T), 1000)
        exact = 1000 * BigFloat(T(0.001))
        @test pattern(fdp(a, b)) == ref_encode(exact, 16, 2)
        # Products below minpos and above maxpos are kept exactly
        tiny, huge = PositPR.smallest(T), typemax(T)
        q = Quire(T)
        qmuladd!(q, tiny, tiny)                    # minpos² (not representable as a posit)
        qmuladd!(q, huge, tiny)                    # 1
        @test Posit(q) == one(T)                   # 1 + minpos² rounds to 1
        @test BigFloat(q) == 1 + BigFloat(tiny)^2  # but the quire holds it exactly
        qmuladd!(q, huge, huge)
        qmulsub!(q, huge, huge)
        @test BigFloat(q) == 1 + BigFloat(tiny)^2  # maxpos² - maxpos² cancels exactly
    end

    @testset "NaR and edge cases" begin
        T = Posit{8, 2}
        @test isNaR(Posit(qadd!(Quire(T), NaR(T))))
        @test isNaR(Posit(qmuladd!(Quire(T), NaR(T), one(T))))
        @test isNaR(fsum(T[1, NaR(T), 2]))
        @test isNaR(fdp(T[1, 2], T[NaR(T), 1]))
        @test isNaR(fma(NaR(T), one(T), one(T)))
        @test iszero(fsum(T[]))
        @test iszero(fdp(T[], T[]))
        @test_throws DimensionMismatch fdp(T[1, 2], T[1])
        # Every posit appears with its opposite, so the exact sum of all posits is 0
        all8 = [from_pattern(u, 8, 2) for u in real_patterns(8)]
        @test iszero(fsum(all8))
    end
end
