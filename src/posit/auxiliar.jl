# Small helpers on bit strings, kept for compatibility with the original API.
# The arithmetic itself works on integers (see rounding.jl and operations.jl).

"Integer value of a bit string (\"\" is 0)."
function to_decimal(bitstr::String)
    isempty(bitstr) && return 0
    return parse(Int, bitstr; base = 2)
end

"Bit string of `decimal_value` with at least `num_bits` characters."
to_binary(decimal_value::Integer, num_bits::Integer) = string(decimal_value; base = 2, pad = num_bits)

"Append `shift` zeros to the right of a bit string."
shift_left(str::String, shift::Int) = str * repeat("0", shift)

"Two's complement of a bit string (invert and add 1)."
function twos_complement(str::String)
    n = length(str)
    v = (~parse(UInt128, str; base = 2) + UInt128(1)) & posit_mask(n)
    return string(v; base = 2, pad = n)
end

"""
    exponent_update(E, es)

Split an expanded exponent `E = 2^es·k + e` into the regime value `k` and the
exponent `e` (Euclidean division: `0 ≤ e < 2^es`).
"""
function exponent_update(e_expanded::Integer, es::Integer)
    k = fld(e_expanded, 1 << es)
    return k, e_expanded - k * (1 << es)
end
