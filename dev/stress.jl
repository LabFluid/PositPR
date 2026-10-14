# Extra check: random operations against the reference oracle (not part of the suite).
using Test, PositPR, Random
include(joinpath(@__DIR__, "..", "test", "reference.jl"))
println("ambiguities: ", length(Test.detect_ambiguities(PositPR)))
for (ps, es, N) in [(8, 2, 20000), (12, 1, 20000), (16, 2, 50000), (32, 2, 50000), (32, 3, 20000), (64, 2, 20000)]
    T = PositType(ps, es); rng = Xoshiro(ps * 10 + es); bad = 0; cnt = 0
    for _ in 1:N
        a = PositPR.from_bits(T, rand(rng, UInt128)); b = PositPR.from_bits(T, rand(rng, UInt128))
        (isNaR(a) || isNaR(b) || iszero(b)) && continue
        for op in (+, -, *, /)
            cnt += 1
            PositPR.bits(op(a, b)) == UInt128(ref_encode(op(BigFloat(a), BigFloat(b)), ps, es)) || (bad += 1)
        end
        if PositPR.sign_bit(a) == 0
            cnt += 1
            PositPR.bits(sqrt(a)) == UInt128(ref_encode(sqrt(BigFloat(a)), ps, es)) || (bad += 1)
        end
    end
    println("$T: $cnt random operations, RNE mismatches = $bad")
end
stop_watchdog()
