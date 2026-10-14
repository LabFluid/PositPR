using PositPR
T = Posit{16,2}
# Catastrophic cancellation: (big + 1) - big
xs = T[2^20, 1, -2^20]
println("sum with + : ", sum(xs), "   fsum (quire): ", fsum(xs))
# Dot product with many small terms
a = [T(0.001) for _ in 1:1000]; b = [T(1) for _ in 1:1000]
naive = foldl((acc, i) -> acc + a[i] * b[i], 1:1000; init = zero(T))
println("dot with + and * : ", naive, "   fdp (quire): ", fdp(a, b), "   exact: ", 1000 * Float64(T(0.001)))
q = Quire(T); qmuladd!(q, T(3), T(4)); qsub!(q, T(2)); println(q, " -> ", Posit(q))
println(fma(T(1.5), T(2.5), T(-3)), "  ", isNaR(Posit(qadd!(Quire(T), NaR(T)))))
