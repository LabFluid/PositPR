"""
    PositDecoder(BP, ps, es)

Posit whose bit pattern is the string `BP` of `ps` characters '0'/'1'.
"""
function PositDecoder(BP::AbstractString, ::Val{ps}, ::Val{es}) where {ps, es}
    length(BP) == ps || throw(ArgumentError("expected $ps bits, got \"$BP\""))
    return from_bits(Posit{ps, es}, parse(UInt128, BP; base = 2))
end

PositDecoder(BP::AbstractString, ps::Integer, es::Integer) = PositDecoder(BP, Val(ps), Val(es))

"""
    PositEncoder(p)

Bit pattern of the posit as a string of `ps` characters '0'/'1'.
"""
PositEncoder(p::Posit{ps, es}) where {ps, es} = string(bits(p); base = 2, pad = ps)

"""
    Posit{ps,es}(s, sn, k, rs, ers, e, fi, fs, f, useed)

Build a posit from the decoded fields of the original struct (kept for
compatibility). Only `s`, `sn`, `k`, `e` and the fraction bits `fi[1:fs]` are
used: the value is `(-1)^s · 2^(2^es·k + e) · (1.fi)₂`.
"""
function Posit{ps, es}(s::Integer, sn::Integer, k::Integer, rs::Integer, ers::Integer,
                       e::Integer, fi::AbstractString, fs::Integer, f::Real,
                       useed::Real) where {ps, es}
    T = PositType(ps, es)
    if sn == 1
        return s == 1 ? NaR(T) : zero(T)
    end
    frac = fs == 0 ? UInt128(0) : parse(UInt128, fi[1:fs]; base = 2)
    M = (UInt128(1) << fs) | frac
    return posit_round(T, s, M, k * (1 << es) + e - fs)
end

"Expanded exponent E = 2^es·k + e."
expand_exponent(a::Posit) = decode_fields(a).E
