# Validação da infraestrutura e do Experimento I contra BigFloat.
#   julia --project=benchmarks benchmarks/test/runtests.jl

using Test
isdefined(@__MODULE__, :ns_run) || include(joinpath(@__DIR__, "..", "near_singular.jl"))

const P16 = PositType(16, 2)
const P32 = PositType(32, 2)
bf(x) = BigFloat(x)
two(k) = BigFloat(2)^k
I_big(n) = [i == j ? bf(1) : bf(0) for i in 1:n, j in 1:n]

@testset "Benchmarks: validação" begin

@testset "conversão com um único arredondamento" begin
    # 1 + 2^-11 + 2^-30 está acima do ponto médio entre 1 e 1+2^-10 (Float16);
    # passando por Float32 o 2^-30 se perde e o empate vai para 1 (par).
    x = 1 + two(-11) + two(-30)
    @test exact(tofmt(Float16, x)) == 1 + two(-10)
    @test exact(Float16(Float32(x))) == 1          # o arredondamento duplo que evitamos
    # o mesmo para Posit{16,2}: 11 bits de fração em [1, 2)
    y = 1 + two(-12) + two(-40)
    @test exact(tofmt(P16, y)) == 1 + two(-11)
    @test exact(tofmt(P16, 1 + two(-12))) == 1      # empate → par
    @test representable(P16, 1 + two(-11)) && !representable(Float16, 1 + two(-11))
    @test representable(Float16, 1 + two(-10)) && representable(P16, 1 + two(-10))
    # Float32 e Posit{32,2}
    @test exact(tofmt(Float32, 1 + two(-24) + two(-60))) == 1 + two(-23)
    @test representable(P32, 1 + two(-27)) && !representable(P32, 1 + two(-28))
end

@testset "CSV: ida e volta" begin
    path = joinpath(mktempdir(), "t.csv")
    rows = [(formato = "Posit{16,2}", j = 3, erro = 0.125, ok = true, s = "a\"b"),
            (formato = "Float16", j = 4, erro = NaN, ok = false, s = "")]
    write_csv(path, rows)
    back = read_csv(path)
    @test back[1][:formato] == "Posit{16,2}" && back[1][:j] == 3 && back[1][:erro] == 0.125
    @test back[1][:ok] == "true" && back[1][:s] == "a\"b"
    @test back[2][:formato] == "Float16" && isnan(back[2][:erro]) && back[2][:s] == ""
end

@testset "dados comuns a todos os formatos" begin
    rng = MersenneTwister(2024)
    Ts = (Float16, P16)
    for _ in 1:300
        x = bf(rand(rng)) * two(rand(rng, -6:6)) * (rand(rng, Bool) ? 1 : -1)
        v, ok = common_value(x, Ts)
        @test ok
        @test all(T -> representable(T, v), Ts)
        @test abs(v - x) <= two(-8) * abs(x)          # perto do original
    end
end

@testset "métricas em BigFloat" begin
    for δ in (two(-1), two(-7), two(-20), two(-100))
        A = BigFloat[1 1; 1 1+δ]
        λ = ((2 + δ) + sqrt(4 + δ^2)) / 2, ((2 + δ) - sqrt(4 + δ^2)) / 2   # simétrica: σ = |λ|
        smax, smin = svals2(A)
        @test smax ≈ λ[1] rtol = 1e-70
        @test smin ≈ λ[2] rtol = 1e-60
        @test cond2(A) ≈ λ[1] / λ[2] rtol = 1e-60
        @test abs(cond2(A) * δ / 4 - 1) < 2δ + 1e-60  # κ ≈ 4/δ
        b = BigFloat[2, 2 + δ]
        x = solve2_exact(A, b)
        @test x == BigFloat[1, 1]
        @test forward_error(x, x) == 0 && residual(A, x, b) == 0 && backward_error(A, x, b) == 0
    end
    @test solve2_exact(BigFloat[1 1; 1 1], BigFloat[1, 1]) === nothing
    @test cond2(BigFloat[1 1; 1 1]) == Inf
    # erro direto e retroativo de uma perturbação conhecida
    A, b, x = BigFloat[2 0; 0 1], BigFloat[2, 1], BigFloat[1, 1]
    @test forward_error(BigFloat[1, 1 + two(-10)], x) ≈ two(-10) / sqrt(bf(2)) rtol = 1e-70
    @test backward_error(A, BigFloat[1, 1 + two(-10)], b) ≈ two(-10) / (2 * sqrt(1 + (1 + two(-10))^2) + sqrt(bf(5))) rtol = 1e-70
end

@testset "álgebra linear genérica" begin
    rng = MersenneTwister(7)
    # em BigFloat o algoritmo reproduz a solução exata (até a precisão de 256 bits)
    for n in (2, 3, 6), pivot in (false, true)
        A = bf.(randn(rng, n, n)) + n * I_big(n)
        xstar = bf.(randn(rng, n))
        b = [sum(A[i, k] * xstar[k] for k in 1:n) for i in 1:n]
        x, st, _ = gauss_solve(A, b; pivot)
        @test st === :ok
        @test maximum(abs.(x .- xstar)) < 1e-70
    end
    # pivô nulo sem pivoteamento; com pivoteamento resolve
    A = BigFloat[0 1; 1 0]; b = BigFloat[3, 5]
    @test gauss_solve(A, b)[2] === :pivo_zero
    @test gauss_solve(A, b; pivot = true)[1] == BigFloat[5, 3]
    # pivoteamento escolhe o maior |aᵢₖ|; empate mantém a ordem
    @test lu_generic(BigFloat[1 2; 3 4]; pivot = true).perm == [2, 1]
    @test lu_generic(BigFloat[1 2; -1 4]; pivot = true).perm == [1, 2]
    # o tipo dos dados é preservado (nenhuma promoção escondida)
    for T in (Float16, P16, Float32, P32)
        x, st, F = gauss_solve(tofmt(T, BigFloat[2 1; 1 3]), tofmt(T, BigFloat[3, 4]))
        @test eltype(x) === T && eltype(F.LU) === T && st === :ok
    end
end

@testset "subdot: fusão = um arredondamento do valor exato" begin
    rng = MersenneTwister(11)
    for _ in 1:300
        n = rand(rng, 1:8)
        a = tofmt(P16, bf.(randn(rng, n)))
        b = tofmt(P16, bf.(randn(rng, n)))
        c = tofmt(P16, bf(randn(rng)))
        ex = exact(c) - sum(exact(a[k]) * exact(b[k]) for k in 1:n)
        @test subdot(c, a, b, true) == tofmt(P16, ex)
        # sem fusão: dois arredondamentos por termo, da esquerda para a direita
        acc = c
        for k in 1:n
            acc = acc - a[k] * b[k]
        end
        @test subdot(c, a, b, false) == acc
    end
    # Float16 ignora `fused` (não há quire para float)
    a, b = Float16[1, 2], Float16[3, 4]
    @test subdot(Float16(1), a, b, true) == subdot(Float16(1), a, b, false)
end

@testset "Experimento I: previsões analíticas" begin
    # Família diádica, x* = [1, 1]: os dados são exatos enquanto 1+δ e 2+δ forem
    # representáveis; a eliminação é exata (l₂₁ = 1, u₂₂ = δ por Sterbenz).
    lastexact = Dict("Float16" => 9, "Posit{16,2}" => 10, "Posit{16,2}+quire" => 10)
    singular = Dict("Float16" => 11, "Posit{16,2}" => 12, "Posit{16,2}+quire" => 12)
    for fmt in FORMATS_16, j in 1:16
        P = ns_problem(:diadica, j, :uns)
        rows, _ = ns_run(fmt, P, (j = j,))
        lu = only(filter(r -> r.algoritmo == "lu", rows))
        pp = only(filter(r -> r.algoritmo == "lu_pp", rows))
        if j <= lastexact[fmt.name]
            @test lu.status == "ok" && lu.erro_orig == 0
        elseif j < singular[fmt.name]
            @test lu.status == "ok" && lu.erro_orig > 0 && lu.erro_arm == 0
        else
            @test lu.status == "pivo_zero" && !lu.distingue
        end
        # o pivoteamento não muda nada aqui (empate |a₂₁| = |a₁₁|)
        @test isequal(lu.erro_orig, pp.erro_orig)
    end
    # nos dados diádicos o quire não tem o que corrigir: resultados idênticos
    for j in 1:16, xk in NS_XEXACT
        P = ns_problem(:diadica, j, xk)
        r1, _ = ns_run(FORMATS_16[2], P, (j = j,))
        r2, _ = ns_run(FORMATS_16[3], P, (j = j,))
        for (a, b) in zip(r1, r2)
            @test isequal(a.erro_orig, b.erro_orig) && a.status == b.status
        end
    end
end

@testset "Experimento I: estrutura rotacionada" begin
    for fam in NS_FAMILIES, j in (1, 5, 10, 20), xk in NS_XEXACT
        Pc = ns_problem(fam, j, xk)
        Pr = ns_problem(fam, j, xk; estrutura = :rotacionada)
        @test cond2(Pr.A) ≈ cond2(Pc.A) rtol = 1e-60     # rotação não muda κ₂
        @test det2(Pr.A) ≈ Pc.delta rtol = 1e-60           # nem o determinante
        @test solve2_exact(Pr.A, Pr.b) ≈ Pr.x rtol = 1e-60  # nem a solução
        x, st, _ = gauss_solve(Pr.A, Pr.b)                  # em BigFloat: exato
        @test st === :ok && forward_error(x, Pr.x) < 1e-60
    end
    # dados genéricos: nenhum formato de 16 bits representa G·A exatamente
    Pr = ns_problem(:diadica, 4, :uns; estrutura = :rotacionada)
    for T in (Float16, P16)
        @test exact(tofmt(T, Pr.A)) != Pr.A
    end
end

@testset "Experimento I: refinamento e modo comum" begin
    # Em BigFloat (sem erro de arredondamento) o refinamento termina na 1ª iteração.
    P = ns_problem(:generica, 8, :alternado)
    F = lu_generic(P.A)
    x, st, it, _ = refine(F, P.b, P.A, P.b)
    @test st in (:residuo_zero, :estacionario) && it <= 1
    @test forward_error(x, P.x) < 1e-70
    # Refinamento contra os dados armazenados converge para a solução armazenada
    # (κ(A_T)·u pequeno em j = 4).
    for fmt in FORMATS_16
        P = ns_problem(:generica, 4, :uns)
        rows, _ = ns_run(fmt, P, (j = 4,))
        r = only(filter(r -> r.algoritmo == "ri_arm_esc", rows))
        @test r.erro_arm < 2.0^-9
    end
    # Modo comum: os mesmos dados (exatos) em todos os formatos
    Ts = unique(f.T for f in FORMATS_16)
    for fam in NS_FAMILIES, j in 1:10, xk in NS_XEXACT
        Pc = ns_common(ns_problem(fam, j, xk), Ts)
        @test Pc.ok == (abs(Pc.delta - delta_value(fam, j)) <= delta_value(fam, j) / 2)
        j <= 8 && @test Pc.ok
        for T in Ts
            @test exact(tofmt(T, Pc.A)) == Pc.A && exact(tofmt(T, Pc.b)) == Pc.b
        end
        rows = reduce(vcat, [ns_run(f, Pc, (j = j,))[1] for f in FORMATS_16])
        @test all(r -> r.dados_exatos, rows)
    end
end

end
