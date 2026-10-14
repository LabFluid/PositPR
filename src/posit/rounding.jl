# The library's single rounding point.
#
# Every operation computes its result EXACTLY (or with spare bits plus a
# "sticky" bit telling whether anything was discarded) in the form
#
#     x = (-1)^s · M · 2^P        (integer M ≥ 0)
#
# and calls `posit_round`, which builds the bit string (regime, exponent,
# fraction), rounds it to nearest with ties to even (RNE), saturates at
# ±minpos/±maxpos (never rounds to 0 or to NaR) and returns the posit.
#
# Integers are kept in a UInt128 whenever the intermediate value fits in 127 bits
# (always the case for ps ≤ 32); otherwise BigInt is used.

"Number of bits of a nonnegative integer."
bitlength(M::UInt128) = 128 - leading_zeros(M)
bitlength(M::BigInt) = iszero(M) ? 0 : ndigits(M; base = 2)
bitlength(M::Integer) = bitlength(big(M))

"`M` as a UInt128 if it (plus `extra` bits of headroom) fits in 127 bits, else as a BigInt."
function fit_int(M::Integer, extra::Integer = 0)
    if bitlength(M) + extra <= 127
        return M isa UInt128 ? M : UInt128(M)
    end
    return big(M)
end

"Integer mantissa with the hidden bit: `|p| = mantissa(p) · 2^(E - fs)`."
function mantissa(p::Posit)
    d = decode_fields(p)
    return (UInt128(1) << d.fs) | d.frac
end

"Exact decomposition `(s, M, P)` with `p = (-1)^s · M · 2^P` (M is a UInt128)."
function exact_parts(p::Posit)
    is_special(p) && return (0, UInt128(0), 0)
    d = decode_fields(p)
    return (d.s, (UInt128(1) << d.fs) | d.frac, d.E - d.fs)
end

function exact_parts(n::Integer)
    s = n < 0 ? 1 : 0
    return (s, fit_int(abs(big(n))), 0)
end

"""
    posit_round(T, s, M, P, sticky = false)

Posit of type `T` nearest to `(-1)^s · M · 2^P`, rounding the bit string with
RNE. `sticky = true` means the true value has nonzero bits below `M` (inexact
result, as in a division).
"""
function posit_round(::Type{PT}, s::Integer, M::Integer, P::Integer,
                     sticky::Bool = false) where {ps, es, PT <: Posit{ps, es}}
    T = PositType(ps, es)
    if iszero(M)
        sticky || return zero(T)
        return from_bits(T, s == 1 ? negate_bits(UInt128(1), ps) : UInt128(1))
    end
    n = bitlength(M)
    E = P + n - 1                              # |x| = 1.f · 2^E
    k = fld(E, 1 << es)
    e = E - k * (1 << es)

    if k > ps - 2                              # above maxpos: saturate
        v = maxpos_bits(ps)
    elseif k < -(ps - 2)                       # below minpos: saturate
        v = UInt128(1)
    else
        rlen = k >= 0 ? k + 2 : -k + 1         # regime length, terminator included
        L = rlen + es + n - 1                  # length of regime | exponent | fraction
        v = L <= 127 ? round_body(UInt128, ps, es, k, e, M, n, L, sticky) :
                       round_body(BigInt, ps, es, k, e, M, n, L, sticky)
    end
    return from_bits(T, s == 1 ? negate_bits(v, ps) : v)
end

"Assemble regime | exponent | fraction, keep ps-1 bits and apply RNE."
function round_body(::Type{I}, ps, es, k, e, M, n, L, sticky) where {I <: Integer}
    one_ = one(I)
    regime = k >= 0 ? ((one_ << (k + 1)) - one_) << 1 : one_
    frac = I(M) & ((one_ << (n - 1)) - one_)   # M without the hidden bit
    body = (regime << (es + n - 1)) | (I(e) << (n - 1)) | frac

    drop = L - (ps - 1)                        # bits that do not fit after the sign
    if drop <= 0
        v = body << (-drop)
        guard = false
    else
        v = body >> drop
        guard = isodd(body >> (drop - 1))      # first discarded bit
        sticky = sticky || !iszero(body & ((one_ << (drop - 1)) - one_))
    end
    if guard && (sticky || isodd(v))
        v += one_                              # RNE: ties go to the even pattern
    end
    v = UInt128(v)
    iszero(v) && (v = UInt128(1))                          # never 0
    v > maxpos_bits(ps) && (v = maxpos_bits(ps))            # never NaR
    return v
end

"""
    posit_round_div(T, s, Mn, Pn, Md, Pd)

Round the quotient `(-1)^s · (Mn·2^Pn) / (Md·2^Pd)`: integer division produces
spare quotient bits and the remainder becomes the sticky bit.
"""
function posit_round_div(::Type{PT}, s, Mn::Integer, Pn::Integer,
                         Md::Integer, Pd::Integer) where {ps, es, PT <: Posit{ps, es}}
    bn, bd = bitlength(Mn), bitlength(Md)
    K = max(0, bd + ps + 3 - bn)               # quotient gets at least ps + 2 bits
    if bn + K <= 127 && bd <= 127
        q, r = divrem(UInt128(Mn) << K, UInt128(Md))
    else
        q, r = divrem(big(Mn) << K, big(Md))
    end
    return posit_round(Posit{ps, es}, s, q, Pn - Pd - K, !iszero(r))
end

# ------------------------------------------------------------------ extended precision
# Elementary functions are evaluated in a posit with GUARD_BITS extra bits and the
# result is rounded back, so that errors accumulated in the series do not show up
# in the final result.

const GUARD_BITS = 16

wide_type(::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = PositType(min(ps + GUARD_BITS, 128), es)
widen_posit(a::Posit{ps, es}) where {ps, es} = Posit2Posit(a, Val(min(ps + GUARD_BITS, 128)), Val(es))
narrow_posit(w::Posit, ::Type{PT}) where {ps, es, PT <: Posit{ps, es}} = Posit2Posit(w, Val(ps), Val(es))
