module PositPR

    include("posit/type.jl")
    include("posit/auxiliar.jl")
    include("posit/encoding_decoding.jl")
    include("posit/rounding.jl")
    include("posit/operations.jl")
    include("posit/conversions.jl")
    include("posit/comparisons.jl")
    include("posit/trig.jl")
    include("posit/random.jl")
    include("posit/quire.jl")

    # `pi(T)` is not exported to avoid clashing with `Base.pi`: use `PositPR.pi(T)`.
    # sin, cos, exp, log, ... extend the Base functions and need no export.
    export
    Posit, Posit8, Posit16, Posit32, Posit64, PositType, useed, #type
    x2Posit, Posit2x, Posit2Posit, PositConverter, #conversions
    PositDecoder, PositEncoder, #encoding_decoding
    nextPosit, prevPosit, isNaR, NaR, #comparisons
    Quire, qadd!, qsub!, qmuladd!, qmulsub!, qneg!, qabs!, qclear!, fdp, fsum, #quire
    shift_left, to_decimal, exponent_update #auxiliar
end
