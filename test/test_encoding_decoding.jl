# PositDecoder (bits -> fields) and PositEncoder (fields -> bits).

@testset "Encoding and decoding" begin

    @testset "Fields of known patterns – Posit{8,1}" begin
        # (pattern, s, k, rs, e, fi)
        cases = [
            ("01000000", 0,  0, 2, 0, "0000"),   #  1
            ("01010000", 0,  0, 2, 1, "0000"),   #  2
            ("01100000", 0,  1, 3, 0, "000"),    #  4
            ("01001000", 0,  0, 2, 0, "1000"),   #  1.5
            ("00100000", 0, -1, 2, 0, "0000"),   #  0.25
            ("00110000", 0, -1, 2, 1, "0000"),   #  0.5
            ("01111111", 0,  6, 7, 0, ""),       #  maxpos = 4^6
            ("00000001", 0, -6, 7, 0, ""),       #  minpos = 4^-6
        ]
        for (bp, s, k, rs, e, fi) in cases
            p = PositDecoder(bp, 8, 1)
            @test p.s == s
            @test p.k == k
            @test p.rs == rs
            @test p.e == e
            @test p.fi == fi || (isempty(fi) && p.fi == "0")
            @test p.fs == length(fi)
            @test p.sn == 0
        end
    end

    @testset "Special cases (zero and NaR)" for (ps, es) in CONFIGS_8
        z = PositDecoder("0"^ps, ps, es)
        n = PositDecoder("1" * "0"^(ps - 1), ps, es)
        @test iszero(z)
        @test isNaR(n)
        @test z.sn == 1 && z.s == 0
        @test n.sn == 1 && n.s == 1
        @test PositEncoder(z) == "0"^ps
        @test PositEncoder(n) == "1" * "0"^(ps - 1)
    end

    @testset "Decoded value = reference, Posit{$ps,$es} (exhaustive)" for (ps, es) in CONFIGS_8
        @test report("value of each pattern", check_all(real_patterns(ps)) do u
            got = Posit2x(from_pattern(u, ps, es))
            want = ref_value(u, ps, es)
            got == want ? nothing : "$(bitstr(u, ps)): got $got, expected $want"
        end)
    end

    @testset "Decoded value = reference, Posit{$ps,$es} (sample)" for (ps, es) in [(16, 1), (16, 2), (32, 2)]
        pats = filter(!=(nar_pattern(ps)), sample_patterns(ps; step = ps == 32 ? 999_983 : 7))
        @test report("value of each pattern", check_all(pats) do u
            got = Posit2x(from_pattern(u, ps, es))
            want = ref_value(u, ps, es)
            got == want ? nothing : "$(bitstr(u, ps)): got $got, expected $want"
        end)
    end

    @testset "Round trip encode(decode(bits)) == bits, Posit{$ps,$es}" for (ps, es) in vcat(CONFIGS_8, [(16, 1), (16, 2), (32, 2)])
        pats = ps == 8 ? collect(all_patterns(ps)) : sample_patterns(ps; step = ps == 32 ? 999_983 : 7)
        @test report("round trip", check_all(pats) do u
            got = PositEncoder(from_pattern(u, ps, es))
            got == bitstr(u, ps) ? nothing : "$(bitstr(u, ps)) -> $got"
        end)
    end

    @testset "Symmetry: decode(-bits) == -decode(bits), Posit{$ps,$es}" for (ps, es) in CONFIGS_8
        @test report("sign symmetry", check_all(real_patterns(ps)) do u
            neg = (~u + 1) & mask(ps)
            a = Posit2x(from_pattern(u, ps, es))
            b = Posit2x(from_pattern(neg, ps, es))
            a == -b ? nothing : "$(bitstr(u, ps)): $a vs $(bitstr(neg, ps)): $b"
        end)
    end

    @testset "Order: increasing patterns (as integers) give increasing values" begin
        for (ps, es) in CONFIGS_8
            sorted = sort(collect(real_patterns(ps)); by = u -> signed_pattern(u, ps))
            vals = [Posit2x(from_pattern(u, ps, es)) for u in sorted]
            @test issorted(vals; lt = <=) && allunique(vals)
        end
    end
end
