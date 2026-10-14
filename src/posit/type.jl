"""
    Posit{ps, es, U} <: Real

A posit with `ps` bits in total and at most `es` exponent bits (2 ≤ ps ≤ 128).

The bit pattern is stored in the smallest unsigned integer that holds `ps` bits:
`U = UInt8` up to 8 bits, `UInt16` up to 16, `UInt32` up to 32, `UInt64` up to 64
and `UInt128` up to 128. `U` is chosen automatically, so in practice one writes

    Posit{16,2}(1.5)        # a Posit{16, 2, UInt16}, 2 bytes
    zero(Posit{8,2})        # a Posit{8, 2, UInt8},   1 byte

`Posit{ps,es}` is the family of all posits with those parameters; the concrete
type is `Posit{ps,es,U}`. For type annotations (array element types, struct
fields) use the concrete type: `PositType(ps, es)` or the aliases `Posit8{es}`,
`Posit16{es}`, `Posit32{es}` and `Posit64{es}`.

A posit is a plain `isbits` value: arrays of posits are compact, `==` compares
integers and `hash` works. The decoded fields of the original implementation
(`p.s`, `p.k`, `p.e`, `p.fi`, ...) are still available as properties and are
computed on demand (see `Base.getproperty` below).

Build posits from real numbers with `Posit{ps,es}(x)` / `x2Posit`, or from bit
patterns with `PositDecoder` / `from_bits`.
"""
struct Posit{ps, es, U <: Unsigned} <: Real
    bits::U

    function Posit{ps, es, U}(::Val{:bits}, bits::Unsigned) where {ps, es, U}
        (ps isa Int && 2 <= ps <= 128) ||
            throw(ArgumentError("Posit{ps,es} requires 2 ≤ ps ≤ 128, got ps = $ps"))
        (es isa Int && es >= 0) ||
            throw(ArgumentError("Posit{ps,es} requires es ≥ 0, got es = $es"))
        U === storage_type(ps) ||
            throw(ArgumentError("Posit{$ps,$es} is stored in $(storage_type(ps)), not $U"))
        return new{ps, es, U}((UInt128(bits) & posit_mask(ps)) % U)
    end
end

"Smallest unsigned integer type that holds `ps` bits."
@inline storage_type(ps::Integer) =
    ps <= 8 ? UInt8 : ps <= 16 ? UInt16 : ps <= 32 ? UInt32 : ps <= 64 ? UInt64 : UInt128

"""
    PositType(ps, es)

Concrete posit type with `ps` bits and `es` exponent bits, e.g.
`PositType(16, 2) == Posit{16, 2, UInt16}`. Use it in type annotations.
"""
@inline PositType(ps::Integer, es::Integer) = Posit{ps, es, storage_type(ps)}
PositType(::Type{<:Posit{ps, es}}) where {ps, es} = PositType(ps, es)

# Concrete aliases for the usual sizes: Posit16{2} == Posit{16, 2, UInt16}
const Posit8{es}  = Posit{8, es, UInt8}
const Posit16{es} = Posit{16, es, UInt16}
const Posit32{es} = Posit{32, es, UInt32}
const Posit64{es} = Posit{64, es, UInt64}

ps(::Posit{p, e}) where {p, e} = p
es(::Posit{p, e}) where {p, e} = e

"Mask with the low `ps` bits set."
posit_mask(ps::Integer) = ps >= 128 ? typemax(UInt128) : (UInt128(1) << ps) - UInt128(1)

"Posit of type `T` (family or concrete type) with the given bit pattern (modulo 2^ps)."
@inline from_bits(::Type{PT}, b::Unsigned) where {ps, es, PT <: Posit{ps, es}} =
    PositType(ps, es)(Val(:bits), b)
from_bits(::Type{PT}, b::Integer) where {ps, es, PT <: Posit{ps, es}} =
    PositType(ps, es)(Val(:bits), UInt128(mod(big(b), big(1) << ps)))

"Bit pattern of the posit, widened to a UInt128 (all internal arithmetic uses UInt128)."
@inline bits(p::Posit) = UInt128(getfield(p, :bits))

# Bit patterns of the special values
@inline nar_bits(ps) = UInt128(1) << (ps - 1)
@inline maxpos_bits(ps) = nar_bits(ps) - UInt128(1)

@inline sign_bit(p::Posit{ps, es}) where {ps, es} = Int(bits(p) >> (ps - 1))
@inline is_special(p::Posit{ps, es}) where {ps, es} = (bits(p) & (nar_bits(ps) - 1)) == 0

"Two's complement of the pattern: the bit pattern of -p (0 and NaR are fixed points)."
@inline negate_bits(u::UInt128, ps) = (~u + UInt128(1)) & posit_mask(ps)

"""
    decode_fields(p)

Decode the bit pattern with integer operations. Returns a NamedTuple with the
sign `s`, regime value `k`, regime size `rs`, number of exponent bits present
`ers`, exponent `e`, fraction size `fs`, fraction bits `frac` (as an integer) and
the expanded exponent `E = 2^es·k + e`. For zero and NaR all fields but `s` are 0.
"""
@inline function decode_fields(p::Posit{ps, es}) where {ps, es}
    u = bits(p)
    s = Int(u >> (ps - 1))
    if is_special(p)
        return (s = s, k = 0, rs = 0, ers = 0, e = 0, fs = 0, frac = UInt128(0), E = 0)
    end
    m = s == 1 ? negate_bits(u, ps) : u        # magnitude
    x = m << (129 - ps)                        # bits after the sign, left-aligned
    if (x >> 127) == 1
        run = leading_ones(x)
        k = run - 1
    else
        run = leading_zeros(x)
        k = -run
    end
    # run equal bits + terminating bit (absent when the regime fills the word)
    rs = min(run + 1, ps - 1)
    rest = x << rs
    nrem = ps - 1 - rs
    ers = min(es, nrem)
    # Exponent bits cut off by the regime are the least significant ones and are
    # zero: reading the top es bits (padded with zeros) gives the right value.
    e = es == 0 ? 0 : Int(rest >> (128 - es))
    fs = nrem - ers
    frac = fs == 0 ? UInt128(0) : (rest << es) >> (128 - fs)
    return (s = s, k = k, rs = rs, ers = ers, e = e, fs = fs, frac = frac, E = k * (1 << es) + e)
end

# Compatibility with the original field-based struct: p.s, p.sn, p.k, p.rs, p.ers,
# p.e, p.fi, p.fs, p.f and p.useed are computed from the bit pattern.
@inline function Base.getproperty(p::Posit{ps, es}, name::Symbol) where {ps, es}
    name === :bits && return getfield(p, :bits)
    name === :s && return sign_bit(p)
    name === :sn && return Int(is_special(p))
    name === :useed && return 2.0^(2^es)
    d = decode_fields(p)
    name === :k && return d.k
    name === :rs && return d.rs
    name === :ers && return d.ers
    name === :e && return d.e
    name === :fs && return d.fs
    name === :fi && return d.fs == 0 ? "0" : string(d.frac; base = 2, pad = d.fs)
    name === :f && return d.fs == 0 ? 0.0 : Float64(d.frac) / 2.0^d.fs
    throw(ArgumentError("Posit has no property $name"))
end

Base.propertynames(::Posit) = (:bits, :s, :sn, :k, :rs, :ers, :e, :fi, :fs, :f, :useed)

Base.copy(p::Posit) = p

function Base.show(io::IO, p::Posit{ps, es}) where {ps, es}
    print(io, "Posit{", ps, ",", es, "}(")
    if isNaR(p)
        print(io, "NaR")
    elseif ps <= 32
        print(io, Posit2x(p))
    else
        print(io, BigFloat(p))
    end
    print(io, ")")
end
