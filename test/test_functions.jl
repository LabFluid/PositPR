# Elementary functions: pi, exp, log, sin, cos, tan, tanh, a^b (b Posit).
#
# For transcendental functions we require an error of at most 1 ULP with respect
# to the correctly rounded value (computed in BigFloat with 1024 bits).
# They all contain `while` loops, so they are called through `gcall` (watchdog.jl).

"Sample of up to `maxn` patterns with |x| in [lo, hi] (exhaustive for 8 bits)."
function patterns_in(ps, es, lo, hi; positive_only = false, maxn = 120)
    pool = ps == 8 ? collect(real_patterns(ps)) : sample_patterns(ps; step = 5)
    pats = [u for u in pool if u != 0 && u != nar_pattern(ps) &&
            lo <= abs(ref_value(u, ps, es)) <= hi && (!positive_only || u < nar_pattern(ps))]
    step = max(1, length(pats) ÷ maxn)
    return pats[1:step:end]
end

function test_function(name, lib_f, big_f, ps, es, pats)
    cases = run_cases(p -> gcall(lib_f, p), big_f, ps, es, [(u,) for u in pats], x -> "$name($x)")
    @test report("$name in Posit{$ps,$es}: no exceptions", check_all(cases) do c
        c.got === nothing ? "$(c.desc): $(c.error)" : nothing
    end)
    @test report("$name in Posit{$ps,$es}: error ≤ 1 ULP", check_all(cases) do c
        c.got !== nothing && ulp_distance(c.got, c.want, ps) <= 1 ? nothing : show_case(c, ps, es)
    end)
end

"`true` if |f(args...) - expected| < tol, guarded against infinite loops."
approx_call(f, args, expected, tol) = abs(Posit2x(gcall(f, args...)) - expected) < tol

@testset "Elementary functions" begin

    @testset "pi" begin
        for (ps, es) in [(8, 2), (16, 1), (16, 2), (32, 2)]
            T = Posit{ps, es}
            got = pattern(gcall(PositPR.pi, T))
            want = ref_encode(BigFloat(π), ps, es)
            @test ulp_distance(got, want, ps) <= 1
        end
    end

    @testset "Notable values" begin
        T = Posit{16, 2}
        @test isone(gcall(exp, zero(T)))
        @test iszero(gcall(log, one(T)))
        @test iszero(gcall(sin, zero(T)))
        @test isone(gcall(cos, zero(T)))
        @test iszero(gcall(tan, zero(T)))
        @test iszero(gcall(tanh, zero(T)))
        @test approx_call(log, (T(2),), log(2), 1e-3)
        @test approx_call(exp, (one(T),), ℯ, 1e-2)
        @test approx_call(sin, (T(π / 2),), 1.0, 1e-3)
        @test approx_call(cos, (T(π),), -1.0, 1e-3)
        @test gcall(sincos, T(0.5)) == (gcall(sin, T(0.5)), gcall(cos, T(0.5)))
        @test_throws Exception gcall(log, zero(T))
        @test_throws Exception gcall(log, T(-1))
    end

    # Domains chosen so that results stay within the range of the type.
    functions = [
        ("exp",  exp,  exp,  0.0,    4.0,    false),
        ("log",  log,  log,  2.0^-8, 2.0^8,  true),
        ("sin",  sin,  sin,  0.0,    8.0,    false),
        ("cos",  cos,  cos,  0.0,    8.0,    false),
        ("tan",  tan,  tan,  0.0,    1.5,    false),
        ("tanh", tanh, tanh, 0.0,    6.0,    false),
    ]

    @testset "$name" for (name, f, g, lo, hi, posonly) in functions
        for (ps, es) in [(8, 2), (16, 2)]
            pats = patterns_in(ps, es, lo, hi; positive_only = posonly, maxn = ps == 8 ? 200 : 60)
            test_function(name, f, g, ps, es, pats)
        end
    end

    @testset "Real power a^b" begin
        T = Posit{16, 2}
        @test gcall(^, T(4), T(0.5)) == T(2)
        @test gcall(^, T(4), T(-0.5)) == T(0.5)
        @test isone(gcall(^, T(3), zero(T)))
        @test iszero(gcall(^, zero(T), T(2)))
        @test approx_call(^, (T(2), T(1.5)), 2^1.5, 1e-2)
        @test approx_call(^, (T(10), T(0.3)), 10^0.3, 1e-2)
        @test approx_call(^, (T(-2), T(3)), -8.0, 1e-2)
        @test_throws Exception gcall(^, T(-2), T(0.5))
    end
end
