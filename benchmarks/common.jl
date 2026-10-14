# Infraestrutura comum dos benchmarks de estabilidade numérica.
#
# Princípios (ver FUNDAMENTACAO.md):
#   * os problemas são construídos em BigFloat (precisão REF_PREC) e convertidos
#     para cada formato uma única vez, diretamente de BigFloat (nunca passando por
#     Float32/Float64, o que causaria arredondamento duplo);
#   * os algoritmos são genéricos no tipo T e só usam operações de T (laços
#     escritos aqui, sem BLAS); variantes de precisão mista são nomeadas como tal;
#   * as métricas são calculadas em BigFloat a partir do valor EXATO dos operandos
#     armazenados, fora da região cronometrada.

using PositPR
using Printf
using Random
using BenchmarkTools
using Statistics

const REF_PREC = 256
setprecision(BigFloat, REF_PREC)

# ----------------------------------------------------------------------------
# Formatos comparados

"""
    Fmt(name, T, fused)

Formato numérico de um experimento. `fused = true` indica que os produtos internos
e as operações `c - a*b` usam acumulação exata (quire / fma) com um único
arredondamento.
"""
struct Fmt
    name::String
    T::DataType
    fused::Bool
end

const FORMATS_16 = [
    Fmt("Float16", Float16, false),
    Fmt("Posit{16,2}", PositType(16, 2), false),
    Fmt("Posit{16,2}+quire", PositType(16, 2), true),
]

const FORMATS_32 = [
    Fmt("Float32", Float32, false),
    Fmt("Posit{32,2}", PositType(32, 2), false),
    Fmt("Posit{32,2}+quire", PositType(32, 2), true),
]

# ----------------------------------------------------------------------------
# Conversões exatas

"Valor exato de um número armazenado (Float16/32/64, Posit ou BigFloat)."
exact(x::Union{Float16, Float32, Float64}) = BigFloat(x)
exact(x::Posit) = BigFloat(x)
exact(x::BigFloat) = x
exact(x::AbstractArray) = map(exact, x)

"""
    tofmt(T, x::BigFloat)

Arredonda `x` para o formato `T` com UM arredondamento (RNE). Float16(::BigFloat)
e Float32(::BigFloat) arredondam corretamente em Julia; Posit{ps,es}(::BigFloat)
usa o arredondamento correto da PositPR.
"""
tofmt(::Type{T}, x::BigFloat) where {T} = T(x)
tofmt(::Type{T}, x::AbstractArray) where {T} = map(v -> tofmt(T, v), x)

"`true` se `x` é exatamente representável em `T`."
representable(::Type{T}, x::BigFloat) where {T} = exact(tofmt(T, x)) == x

"""
    common_value(x, Ts)

Valor representável em TODOS os formatos `Ts`, obtido arredondando `x` para cada
formato até chegar a um ponto fixo. Devolve `(valor, ok)`; `ok = false` se não
houver ponto fixo em poucas iterações.
"""
function common_value(x::BigFloat, Ts)
    v = x
    for _ in 1:8
        w = v
        for T in Ts
            w = exact(tofmt(T, w))
        end
        all(T -> representable(T, w), Ts) && return (w, true)
        v = w
    end
    return (v, false)
end

common_value(x::AbstractArray, Ts) = (r = map(v -> common_value(v, Ts), x); (first.(r), all(last.(r))))

# ----------------------------------------------------------------------------
# Operações que podem ser fundidas

"`c - a*b` com dois arredondamentos (`fused = false`) ou um só (`fused = true`)."
@inline fms(a, b, c, fused::Bool) = fused ? fma(-a, b, c) : c - a * b

"""
    dot_T(x, y, fused)

Produto interno no formato de `x`: laço com `acc += x[i]*y[i]` (dois arredondamentos
por termo) ou, se `fused`, `fdp` (quire: um único arredondamento no fim).
"""
function dot_T(x::AbstractVector{T}, y::AbstractVector{T}, fused::Bool) where {T}
    fused && return fdp(x, y)
    acc = zero(T)
    @inbounds for i in eachindex(x, y)
        acc += x[i] * y[i]
    end
    return acc
end

# ----------------------------------------------------------------------------
# Rastreamento das magnitudes intermediárias (só na execução de diagnóstico)

"""
    MagTracker()

Registra `log2|v|` dos valores intermediários não nulos de um algoritmo, para
saber se o cálculo fica perto de 1 ou passa por regiões em que a precisão cônica
do posit diminui. Nas execuções cronometradas passa-se `nothing`.
"""
mutable struct MagTracker
    lo::Float64
    hi::Float64
    n::Int
    outside::Int      # valores fora de [2^-4, 2^4) (onde Posit{16,2} tem 11 bits)
end
MagTracker() = MagTracker(Inf, -Inf, 0, 0)

@inline track!(::Nothing, v) = nothing
function track!(t::MagTracker, v)
    (iszero(v) || isnan(Float64(exact(v)))) && return nothing
    l = Float64(log2(abs(exact(v))))
    t.lo = min(t.lo, l)
    t.hi = max(t.hi, l)
    t.n += 1
    (l < -4 || l >= 4) && (t.outside += 1)
    return nothing
end

# ----------------------------------------------------------------------------
# Métricas em BigFloat

norm2(v::AbstractVector{BigFloat}) = sqrt(sum(abs2, v))
matvec(A::AbstractMatrix{BigFloat}, x::AbstractVector{BigFloat}) =
    [sum(A[i, j] * x[j] for j in axes(A, 2)) for i in axes(A, 1)]

"Valores singulares de uma matriz 2×2 (fórmula fechada, em BigFloat)."
function svals2(A::AbstractMatrix{BigFloat})
    a, b, c, d = A[1, 1], A[1, 2], A[2, 1], A[2, 2]
    s1 = a^2 + b^2 + c^2 + d^2
    det = a * d - b * c
    disc = sqrt(max(s1^2 - 4det^2, big(0)))
    smax = sqrt((s1 + disc) / 2)
    smin = abs(det) / smax                   # σmax·σmin = |det|, sem cancelamento
    return smax, smin
end

opnorm2(A::AbstractMatrix{BigFloat}) = size(A) == (2, 2) ? svals2(A)[1] :
    error("opnorm2 só implementada para 2×2 neste experimento")
cond2(A::AbstractMatrix{BigFloat}) = size(A) == (2, 2) ? (s = svals2(A); s[2] == 0 ? big(Inf) : s[1] / s[2]) :
    error("cond2 só implementada para 2×2 neste experimento")

"Determinante exato (em BigFloat) de uma matriz 2×2."
det2(A::AbstractMatrix{BigFloat}) = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]

"Solução exata (em BigFloat) de um sistema 2×2, pela regra de Cramer."
function solve2_exact(A::AbstractMatrix{BigFloat}, b::AbstractVector{BigFloat})
    det = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    iszero(det) && return nothing
    return [(b[1] * A[2, 2] - A[1, 2] * b[2]) / det, (A[1, 1] * b[2] - A[2, 1] * b[1]) / det]
end

"Erro direto relativo ‖x - x*‖₂ / ‖x*‖₂."
forward_error(x, xstar) = norm2(exact(x) .- xstar) / norm2(xstar)

"Resíduo relativo ‖b - Ax‖₂ / ‖b‖₂ (b ≠ 0) ou ‖b - Ax‖₂ (b = 0)."
function residual(A, x, b)
    r = norm2(b .- matvec(A, exact(x)))
    nb = norm2(b)
    return iszero(nb) ? r : r / nb
end

"Erro retroativo normwise η = ‖b - Ax‖ / (‖A‖‖x‖ + ‖b‖)."
backward_error(A, x, b) = norm2(b .- matvec(A, exact(x))) / (opnorm2(A) * norm2(exact(x)) + norm2(b))

"Converte BigFloat em Float64 para salvar (Inf/NaN preservados)."
f64(x) = x === nothing ? NaN : Float64(x)

# ----------------------------------------------------------------------------
# Temporização (fora dos diagnósticos)

"""
    time_ns_per_call(f, args...)

Tempo mediano de `f(args...)` em ns, com aquecimento (BenchmarkTools compila e
descarta as primeiras execuções). Nenhum diagnóstico em BigFloat roda aqui.
"""
function time_ns_per_call(f, args...)
    f(args...)                                   # aquecimento / compilação
    b = @benchmark $f($(args)...) samples = 200 evals = 20 seconds = 1
    return median(b).time, median(b).allocs
end

# ----------------------------------------------------------------------------
# CSV simples (sem dependências)

csvfield(x::AbstractString) = occursin(r"[,\"]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
csvfield(x::AbstractFloat) = isnan(x) ? "NaN" : isinf(x) ? (x > 0 ? "Inf" : "-Inf") : @sprintf("%.17g", x)
csvfield(x) = string(x)

"Grava um vetor de NamedTuples (mesmas chaves) em CSV."
function write_csv(path, rows)
    isempty(rows) && return path
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), ","))
        for r in rows
            println(io, join(csvfield.(values(r)), ","))
        end
    end
    return path
end

"Lê um CSV gravado por `write_csv` como vetor de Dicts (números convertidos)."
function read_csv(path)
    lines = readlines(path)
    hdr = Symbol.(csvsplit(lines[1]))
    parse_field(s) = (v = tryparse(Float64, s); v === nothing ? String(s) : v)
    return [Dict(zip(hdr, parse_field.(csvsplit(l)))) for l in lines[2:end]]
end

"Divide uma linha CSV respeitando campos entre aspas (\"Posit{16,2}\")."
function csvsplit(line::AbstractString)
    fields, buf, inq, i = String[], IOBuffer(), false, firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inq && c == '"' && i < lastindex(line) && line[nextind(line, i)] == '"'
            write(buf, '"'); i = nextind(line, i)           # aspas escapadas ""
        elseif c == '"'
            inq = !inq
        elseif c == ',' && !inq
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end
