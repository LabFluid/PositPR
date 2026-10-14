# Arithmetic: +, -, *, /, sqrt, ^, operations with integers, shifts, mod, max, round.

"Result of a library operation compared with the reference value."
struct OpCase
    desc::String          # "a op b"
    got::Union{UInt64, Nothing}
    canonical::Bool       # result == PositDecoder(PositEncoder(result))
    error::String
    want::UInt64          # correctly rounded pattern
    representable::Bool   # the exact result is a posit
end

function run_cases(lib_op, exact_op, ps, es, inputs, desc; maxhangs = 2)
    cases = OpCase[]
    nhangs = 0
    for args in inputs
        exact = exact_op(map(u -> big_value(u, ps, es), args)...)
        want = ref_encode(exact, ps, es)
        representable = isfinite(exact) && big_value(want, ps, es) == exact
        got, canonical, err = nothing, false, ""
        if nhangs >= maxhangs
            err = "not run (after $nhangs infinite loops)"
        else
            try
                r = lib_op(map(u -> from_pattern(u, ps, es), args)...)
                got = pattern(r)
                canonical = r == from_pattern(got, ps, es)
            catch e
                e isa HangError && (nhangs += 1)
                err = e isa HangError ? "infinite loop (> $(WATCHDOG_TIMEOUT) s)" :
                      first(sprint(showerror, e), 100)
            end
        end
        d = desc(map(u -> ref_value(u, ps, es), args)...)
        push!(cases, OpCase(d, got, canonical, err, want, representable))
    end
    return cases
end

show_case(c, ps, es) = "$(c.desc): got $(c.got === nothing ? c.error : ref_value(c.got, ps, es)), " *
                       "expected $(ref_value(c.want, ps, es))"

"The four levels of strictness described in runtests.jl."
function test_levels(cases, ps, es)
    @test report("no exceptions", check_all(cases) do c
        c.got === nothing ? "$(c.desc): $(c.error)" : nothing
    end)
    @test report("result in canonical form", check_all(cases) do c
        c.got === nothing || c.canonical ? nothing : "$(c.desc): non-canonical internal fields"
    end)
    @test report("exact when representable", check_all(cases) do c
        !c.representable || c.got == c.want ? nothing : show_case(c, ps, es)
    end)
    @test report("correct rounding (RNE)", check_all(cases) do c
        c.got == c.want ? nothing : show_case(c, ps, es)
    end)
    @test report("error ≤ 1 ULP", check_all(cases) do c
        c.got !== nothing && ulp_distance(c.got, c.want, ps) <= 1 ? nothing : show_case(c, ps, es)
    end)
end

const BINARY_OPS = [("+", +), ("-", -), ("*", *), ("/", /)]

"Random pattern pairs (no NaR) for types too large for exhaustive testing."
function random_pairs(ps, n; rng = Xoshiro(7))
    draw() = (u = rand(rng, UInt64) & mask(ps); u == nar_pattern(ps) ? UInt64(0) : u)
    return [(draw(), draw()) for _ in 1:n]
end

@testset "Arithmetic operations" begin

    @testset "Posit{$ps,$es} – all pairs" for (ps, es) in CONFIGS_8
        pats = collect(real_patterns(ps))
        for (name, op) in BINARY_OPS
            @testset "a $name b" begin
                pairs = [(a, b) for a in pats for b in pats if !(name == "/" && b == 0)]
                cases = run_cases(op, op, ps, es, pairs, (x, y) -> "$x $name $y")
                test_levels(cases, ps, es)
            end
        end
    end

    @testset "Posit{$ps,$es} – random pairs" for (ps, es, n) in [(16, 1, 3000), (16, 2, 3000), (32, 2, 500)]
        for (name, op) in BINARY_OPS
            @testset "a $name b" begin
                pairs = [(a, b) for (a, b) in random_pairs(ps, n) if !(name == "/" && b == 0)]
                cases = run_cases(op, op, ps, es, pairs, (x, y) -> "$x $name $y")
                test_levels(cases, ps, es)
            end
        end
    end

    @testset "Simple everyday cases" begin
        T = Posit{16, 2}
        x(v) = T(v)
        @test x(1) + x(1) == x(2)
        @test x(1.5) + x(2.25) == x(3.75)
        @test x(3) - x(1) == x(2)
        @test x(1) - x(3) == x(-2)
        @test zero(T) - x(3) == x(-3)
        @test x(3) - zero(T) == x(3)
        @test x(1.5) * x(2.5) == x(3.75)
        @test x(-1.5) * x(2.5) == x(-3.75)
        @test x(3) / x(2) == x(1.5)
        @test x(1) / x(4) == x(0.25)
        @test x(-6) / x(3) == x(-2)
        @test x(0.1) + x(0.2) == T(0.1 + 0.2) || abs(Posit2x(x(0.1) + x(0.2)) - 0.3) < 1e-3
    end

    @testset "Algebraic properties, Posit{$ps,$es} (exhaustive)" for (ps, es) in [(8, 1), (8, 2)]
        pats = collect(real_patterns(ps))
        P(u) = from_pattern(u, ps, es)
        T = Posit{ps, es}
        pairs = [(a, b) for a in pats for b in pats]
        @test report("a + b == b + a", check_all(pairs) do (a, b)
            P(a) + P(b) == P(b) + P(a) ? nothing : "$(ref_value(a, ps, es)), $(ref_value(b, ps, es))"
        end)
        @test report("a * b == b * a", check_all(pairs) do (a, b)
            P(a) * P(b) == P(b) * P(a) ? nothing : "$(ref_value(a, ps, es)), $(ref_value(b, ps, es))"
        end)
        @test report("a - b == a + (-b)", check_all(pairs) do (a, b)
            P(a) - P(b) == P(a) + (-P(b)) ? nothing : "$(ref_value(a, ps, es)), $(ref_value(b, ps, es))"
        end)
        @test report("a + 0 == a, a * 1 == a, a / 1 == a, a - 0 == a, 0 - a == -a", check_all(pats) do a
            p = P(a)
            (p + zero(T) == p && p * one(T) == p && p / one(T) == p && p - zero(T) == p &&
             zero(T) - p == -p) ? nothing : "a = $(ref_value(a, ps, es))"
        end)
        @test report("a - a == 0 e a + (-a) == 0", check_all(pats) do a
            iszero(P(a) - P(a)) && iszero(P(a) + (-P(a))) ? nothing : "a = $(ref_value(a, ps, es))"
        end)
        @test report("a * 0 == 0", check_all(pats) do a
            iszero(P(a) * zero(T)) ? nothing : "a = $(ref_value(a, ps, es))"
        end)
        @test report("a / a == 1 (a ≠ 0)", check_all(filter(!iszero, pats)) do a
            isone(P(a) / P(a)) ? nothing : "a = $(ref_value(a, ps, es)): $(Posit2x(P(a) / P(a)))"
        end)
        @test report("-(-a) == a", check_all(pats) do a
            -(-P(a)) == P(a) ? nothing : "a = $(ref_value(a, ps, es))"
        end)
    end

    @testset "NaR and division by zero" begin
        T = Posit{8, 2}
        nar = NaR(T)
        # The Posit standard defines NaR as the result of x/0 and as a "contagious"
        # value. The library may throw an exception, but must never silently
        # return a real number.
        nar_or_throws(f) = try isNaR(f()) catch; true end
        @test nar_or_throws(() -> one(T) / zero(T))
        @test nar_or_throws(() -> zero(T) / zero(T))
        @test nar_or_throws(() -> one(T) / 0)
        for x in (zero(T), one(T), T(-3.5))
            @test nar_or_throws(() -> nar + x)
            @test nar_or_throws(() -> x + nar)
            @test nar_or_throws(() -> nar - x)
            @test nar_or_throws(() -> nar * x)
            @test nar_or_throws(() -> x * nar)
            @test nar_or_throws(() -> nar / x)
        end
        @test nar_or_throws(() -> sqrt(nar))
    end

    @testset "sqrt, Posit{$ps,$es} (all positives)" for (ps, es) in vcat(CONFIGS_8, [(16, 2)])
        pos = ps == 8 ? [u for u in real_patterns(ps) if 0 <= u < nar_pattern(ps)] :
              [u for u in sample_patterns(ps; step = 11) if 0 <= u < nar_pattern(ps)]
        cases = run_cases(sqrt, sqrt, ps, es, [(u,) for u in pos], x -> "sqrt($x)")
        test_levels(cases, ps, es)
    end

    @testset "sqrt: special cases" begin
        T = Posit{16, 2}
        @test sqrt(T(4)) == T(2)
        @test sqrt(T(0.25)) == T(0.5)
        @test sqrt(T(16)) == T(4)
        @test sqrt(T(2.25)) == T(1.5)
        @test isone(sqrt(one(T)))
        @test iszero(sqrt(zero(T)))
        @test_throws Exception sqrt(T(-1))
    end

    @testset "Integer power a^n" begin
        T = Posit{16, 2}
        @test T(2)^0 == one(T)
        @test T(2)^1 == T(2)
        @test T(2)^3 == T(8)
        @test T(-2)^3 == T(-8)
        @test T(-2)^2 == T(4)
        @test T(0.5)^2 == T(0.25)
        @test T(2)^-1 == T(0.5)
        @test T(2)^-3 == T(0.125)
        @test T(1.5)^2 == T(1.5) * T(1.5)
        @test T(3)^4 == T(81)
        for (ps, es) in [(8, 2), (16, 2)]
            pats = [u for u in real_patterns(ps) if u != 0 && abs(ref_value(u, ps, es)) <= 8][1:3:end]
            for n in (2, 3)
                cases = run_cases(p -> p^n, x -> x^n, ps, es, [(u,) for u in pats], x -> "($x)^$n")
                @test report("Posit{$ps,$es} x^$n with error ≤ 1 ULP", check_all(cases) do c
                    c.got !== nothing && ulp_distance(c.got, c.want, ps) <= 1 ? nothing : show_case(c, ps, es)
                end)
            end
        end
    end

    @testset "Operations with integers" begin
        for (ps, es) in [(8, 2), (16, 2)]
            T = Posit{ps, es}
            pats = [u for u in real_patterns(ps) if abs(ref_value(u, ps, es)) <= 64][1:7:end]
            for n in (-5, -3, -1, 0, 1, 2, 3, 6, 10), (name, f, g) in [
                    ("a + n", (p, n) -> p + n, (x, n) -> x + n),
                    ("n + a", (p, n) -> n + p, (x, n) -> n + x),
                    ("a - n", (p, n) -> p - n, (x, n) -> x - n),
                    ("a * n", (p, n) -> p * n, (x, n) -> x * n),
                    ("n * a", (p, n) -> n * p, (x, n) -> n * x),
                    ("a / n", (p, n) -> p / n, (x, n) -> x / n),
                    ("n / a", (p, n) -> n / p, (x, n) -> n / x)]
                (n == 0 && name == "a / n") && continue
                inputs = name == "n / a" ? filter(!=(0), pats) : pats
                cases = run_cases(p -> f(p, n), x -> g(x, n), ps, es, [(u,) for u in inputs],
                                  x -> replace(replace(name, "a" => "($x)"), "n" => "$n"))
                @test report("Posit{$ps,$es} $name (n = $n), error ≤ 1 ULP", check_all(cases) do c
                    c.got !== nothing && ulp_distance(c.got, c.want, ps) <= 1 ? nothing : show_case(c, ps, es)
                end)
            end
        end
    end

    @testset "Shifts << and >> (multiplication/division by 2^n)" begin
        for (ps, es) in [(8, 1), (8, 2), (16, 2)]
            pats = [u for u in real_patterns(ps) if u != 0 && 2.0^-6 <= abs(ref_value(u, ps, es)) <= 2.0^6]
            for n in 1:3
                for (name, f, g) in [("<<", p -> p << n, x -> x * 2^n), (">>", p -> p >> n, x -> x / 2^n)]
                    cases = run_cases(f, g, ps, es, [(u,) for u in pats], x -> "($x) $name $n")
                    @test report("Posit{$ps,$es} a $name $n == a·2^(±$n)", check_all(cases) do c
                        c.got == c.want ? nothing : show_case(c, ps, es)
                    end)
                end
            end
        end
    end

    @testset "mod, max, round" begin
        T = Posit{16, 2}
        @test mod(T(5), T(3)) == T(2)
        @test mod(T(5.5), T(2)) == T(1.5)
        @test mod(T(-5), T(3)) == T(1)              # same convention as Base.mod
        @test mod(T(7), T(-3)) == T(-2)
        @test_throws Exception mod(T(1), zero(T))
        @test max(T(1), T(2)) == T(2)
        @test max(T(-1), T(-2)) == T(-1)
        @test max(T(1), 3.0) == T(3)
        @test max(-1.0, T(-2)) == T(-1)
        # Base.round must round to an integer (ties to even), as for Float64.
        @test round(T(2.5)) == T(2)
        @test round(T(3.5)) == T(4)
        @test round(T(-1.25)) == T(-1)
        @test round(T(7)) == T(7)
    end
end
