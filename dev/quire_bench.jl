using Test, PositPR
println("ambiguities: ", length(Test.detect_ambiguities(PositPR)))
for T in (Posit8{2}, Posit16{2}, Posit32{2})
    xs = [T(randn()) for _ in 1:10_000]; ys = [T(randn()) for _ in 1:10_000]
    fdp(xs, ys); sum(xs .* ys)
    t_q = minimum(@elapsed(fdp(xs, ys)) for _ in 1:5) / length(xs) * 1e9
    a_q = @allocated fdp(xs, ys)
    naive(xs, ys) = (acc = zero(eltype(xs)); for i in eachindex(xs); acc += xs[i] * ys[i]; end; acc)
    naive(xs, ys)
    t_n = minimum(@elapsed(naive(xs, ys)) for _ in 1:5) / length(xs) * 1e9
    println(T, ": fdp $(round(t_q, digits=1)) ns/termo ($(round(a_q / length(xs), digits=1)) B/termo)   ",
            "acc += x*y $(round(t_n, digits=1)) ns/termo   quire = $(PositPR.standard_quire_bits(PositPR.ps(xs[1]))) bits no padrão")
end
