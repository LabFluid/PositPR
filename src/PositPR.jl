module PositPR

    include("posit/type.jl")
    include("posit/auxiliar.jl")
    include("posit/encoding_decoding.jl")
    include("posit/operations.jl")
    include("posit/conversions.jl")
    include("posit/comparisons.jl")
    include("posit/trig.jl")
    include("posit/random.jl")

    export 
    Posit, useed #type
    x2Posit, Posit2x, Posit2Posit, PositConverter, #conversions
    PositDecoder, PositEncoder, #encoding_decoding
    nextPosit, prevPosit, isNaR, NaR, #comparisons
    pi, sin, cos, tan, tanh #trig
    shift_left, to_decimal, exponent_update #auxiliar
end
 
