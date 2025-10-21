
import Base: >>
function >>(a::Posit{ps, es}, b::Int64) where {ps, es}
    if a.e >= b
        return Posit{ps, es}(a.s, a.sn, a.k, a.rs, a.ers, a.e - b, a.fi, a.fs, a.f, a.useed)
    else
        two_to_es = 1<<es
        reduction_in_k = Int(ceil((b-a.e)/two_to_es))
        return Posit{ps, es}(a.s, a.sn, a.k - reduction_in_k, a.rs, a.ers, a.e + reduction_in_k*two_to_es - b, a.fi, a.fs, a.f, a.useed)
    end
    y
end

import Base: <<
function <<(a::Posit{ps, es}, b::Int64) where {ps, es}
    two_to_es = 1<<es
    addition_in_k = div(b+a.e, two_to_es)
    return Posit{ps, es}(a.s, a.sn, a.k + addition_in_k, a.rs, a.ers, (b+a.e) % two_to_es, a.fi, a.fs, a.f, a.useed)
end

import Base: round
function round(p::Posit{ps, es}) where {ps, es}
    # If there are no fraction bits, nothing to round
    if p.fs == 0
        return Posit{ps, es}(p.s, p.sn, p.k, p.rs, p.ers, p.e, "", 0, 0.0, p.useed)
    end

    fi = p.fi
    e = p.e
    k = p.k

    # Check if we need to round (there are discarded bits)
    if length(fi) > p.fs
        first_discarded_bit = fi[p.fs + 1]
        fi = fi[1:p.fs]  # truncate

        if first_discarded_bit == '1'
            # Add 1 to the last bit of the current fi
            rounded, carry = sum_frac(fi, "0"^(p.fs - 1) * "1", 0)
            fi = rounded
            e += carry
            k, e = exponent_update(k * 2^es + e, es)
        end
    else
        fi = join_frac(fi, p.fs)
    end

    f = to_decimal(fi) / 2.0^p.fs

    return Posit{ps, es}(p.s, p.sn, k, p.rs, p.ers, e, fi, p.fs, f, p.useed)
end

import Base: +
function +(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    # Special cases
    if iszero(a)
        return b
    elseif iszero(b)
        return a
    end

    if a == twos_complement(b)
        return zero(Posit{ps, es})
    end

    # Ensure a is the larger one (in magnitude)
    if abs(a) < abs(b)
        b, a = a, b
    end

    # Convert to common format if needed
    a, b = PositConverter(a, b)

    # If signs differ, use subtraction logic
    if a.s != b.s
        if a > b
            return a - twos_complement(b)
        else
            return b - twos_complement(a)
        end
    end

    # Aligned addition
    temp_a = 2^es * a.k + a.e
    temp_b = 2^es * b.k + b.e
    shifts_in_f = temp_a - temp_b
    fi_sum, new_e = sum_frac_plus(shift_right("1" * b.fi, shifts_in_f), "1" * a.fi, temp_a)
    new_k, new_e = exponent_update(new_e, es)

    new_rs = new_k ≥ 0 ? new_k + 2 : -new_k + 1
    if new_rs >= ps - 1
        return new_k ≥ 0 ? typemax(Posit{ps, es}) : zero(Posit{ps, es})
    end

    new_ers = max(0, min(es, ps - new_rs - 1))
    new_fs = max(0, ps - new_rs - new_ers - 1)
    fi_sum = join_frac(fi_sum, new_fs)

    f_val = to_decimal(fi_sum) / 2^new_fs

    # Construct new Posit
    p = Posit{ps, es}(a.s, 0, new_k, new_rs, new_ers, new_e, fi_sum, new_fs, f_val, a.useed)
    return round(p)
end

import Base: -
function -(a::Posit{ps, es}) where {ps, es}
    return twos_complement(a)
end
	
function -(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    # Special cases
    if iszero(a)
        return b
    elseif iszero(b)
        return a
    end

    if a == b
        return zero(Posit{ps, es})
    end

    # Sign difference → convert to addition
    if a.s != b.s
        p = abs(a) + abs(b)
        return Posit{ps, es}(p.s == a.s ? p.s : a.s, 0, p.k, p.rs, p.ers, p.e, p.fi, p.fs, p.f, a.useed)
    end

    # Ensure a > b in magnitude
    a, b = PositConverter(a, b)

    s = a.s
    if a < b
        if a.s == 1
            a, b = abs(a), abs(b)
        else
            s = 1
            b, a = a, b
        end
    elseif a.s == 1
        s = 0
        a, b = abs(b), abs(a)
    end

    # Aligned subtraction
    temp_a = 2^es * a.k + a.e
    temp_b = 2^es * b.k + b.e
    shifts_in_f = temp_a - temp_b

    fi, e = sub_frac("1" * a.fi, shift_right("1" * b.fi * "0", shifts_in_f), temp_a)
    if fi isa Nothing
        return zero(Posit{ps, es})
    end
    k, e = exponent_update(e, es)

    rs = k ≥ 0 ? k + 2 : -k + 1
    if rs >= ps - 1
        return k ≥ 0 ? typemax(Posit{ps, es}) : zero(Posit{ps, es})
    end

    ers = max(0, min(es, ps - rs - 1))
    fs = ps - rs - ers - 1
    fi = rpad(fi, fs, '0')[1:fs]
    f = to_decimal(fi) / 2.0^fs

    p = Posit{ps, es}(s, 0, k, rs, ers, e, fi, fs, f, a.useed)
    return round(p)
end

import Base: *
function *(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    if iszero(a) || iszero(b)
        return zero(Posit{ps, es})
    elseif isone(a)
        return b
    elseif isone(b)
        return a
    end

    if a.s != b.s
		s = 1
	else
		s = 0
	end

    fi = "0"
    e = 0

    # Fast path: one of the fractions is a power of two
    if ispower_oftwo(a.fi)
        fi = b.fi
    elseif ispower_oftwo(b.fi)
        fi = a.fi
    else
        # General case: multiply and add hidden bits
        fi_prod = mul_frac(a.fi, b.fi)
        sum, temp_e = sum_frac(a.fi, b.fi, 0)

        # Adjust using hidden bit representation (e.g., 10..., 1..., 0...)
        fi, e = sum_frac(string(temp_e) * sum, "0" * fi_prod, 0)

        if e == 1
            fi = '1' * fi[2:end]
        elseif fi[1] == '0'
            fi = fi[2:end]
        else
            fi = '0' * fi[2:end]
            e += 1
        end
    end

    # Compute new exponent: k + e
    e_expanded = expand_exponent(a) + expand_exponent(b) + e
    k, e = exponent_update(e_expanded, es)

    # Regime size
    rs = k ≥ 0 ? k + 2 : -k + 1
    if rs >= ps - 1
        return k ≥ 0 ? typemax(Posit{ps, es}) : zero(Posit{ps, es})
    end

    ers = max(0, min(es, ps - rs - 1))
    fs = ps - rs - ers - 1

    fi = join_frac(fi, fs)
    f = to_decimal(fi) / 2.0^fs

    p = Posit{ps, es}(s, 0, k, rs, ers, e, fi, fs, f, a.useed)
    return round(p)
end

function *(a::Posit{ps, es}, b::Integer) where {ps, es}
    # Zero cases
    if iszero(a)
        return a
    end

    if b == 0
        return zero(Posit{ps, es})
    end

    # Handle negative integer: a * (-b) = -(a * b)
    if b < 0
        return -(a * (-b))
    end

    # Fast path: a * 1 or a == 1
    if isone(a)
        return x2Posit(b, Val(ps), Val(es))
    end

    if b == 1
        return a
    end

    # If b is a power of two, use left shift
    if ispower_oftwo(b)
        n = floor(Int, log2(b))
        return a << n
    end

    # accumulate shifted versions
    b_bin = lstrip(bitstring(b), '0')
    p = zero(Posit{ps, es})

    for (i, char) in enumerate(reverse(b_bin))
        enable_sum = parse(Int, char)
        p += a * (enable_sum * 2^(i - 1))
    end

    rs = p.k ≥ 0 ? p.k + 2 : -p.k + 1
    if rs >= ps - 1
        return p.k ≥ 0 ? typemax(Posit{ps, es}) : zero(Posit{ps, es})
    end
    ers = max(0, min(es, ps - rs - 1))
    fs = ps - rs - ers - 1

    fi = p.fi[1:fs]
    f = to_decimal(fi) / 2.0^fs

    return Posit{ps, es}(p.s, 0, p.k, rs, ers, p.e, fi, fs, f, p.useed)
end

import Base: /
function /(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    if iszero(b)
        throw("can't divide by 0")
    end

    if iszero(a)
        return a
    end

    if isone(b)
        return a
    end

    if a.s != b.s
		s = 1
	else
		s = 0
	end

    fi = "0"
    e = 0
    k = 0

    # If one of the fractions is power of two, shortcut
    if ispower_oftwo(a.fi) || ispower_oftwo(b.fi)
        fi = ispower_oftwo(a.fi) ? b.fi : a.fi
    end

    # Exponent calculation (k and e)
    e_expanded = expand_exponent(a) - expand_exponent(b)
    k, e = exponent_update(e_expanded, es)

    # Regime size
    rs = k ≥ 0 ? k + 2 : -k + 1
    if rs >= ps - 1
        return k ≥ 0 ? typemax(Posit{ps, es}) : zero(Posit{ps, es})
    end

    # Exponent and fraction sizes
    ers = max(0, min(es, ps - rs - 1))
    fs = max(0, ps - rs - ers - 1)

    # Fraction division logic
    if a.fi == b.fi
        fi = "0"
    else
        fi = div_frac("1" * a.fi, "1" * b.fi, fs + 2)

        # Normalize: ensure MSB is 1
        if fi[1] == '0'
            e -= 1
            k, e = exponent_update(2^es * k + e, es)
            fi = fi[3:end]
        else
            fi = fi[2:end-1]
        end
    end

    f = to_decimal(fi) / 2.0^fs
    p = Posit{ps, es}(s, 0, k, rs, ers, e, fi, fs, f, a.useed)
    return round(p)
end

function /(a::Posit{ps, es}, b::Integer) where {ps, es}
    if b == 0
        error("can't divide by 0")
    end

    if iszero(a)
        return zero(Posit{ps, es})
    end

    if b == 1
        return a
    end

    # If b is a power of two, use right shift
    if ispower_oftwo(b)
        n = floor(Int, log2(b))
        return a >> n
    end

    # General case: factorize b and divide iteratively
    factors = factorize(b)
    p = copy(a)

    for (factor, count) in factors
        if ispower_oftwo(factor)
            n = floor(Int, log2(factor))
            p = p >> n
        else
            p = p / x2Posit(factor, Val(ps), Val(es))
        end
    end

    return p
end

/(b::Integer, a::Posit{ps, es}) where {ps, es} = x2Posit(b, Val(ps), Val(es)) / a

import Base: sqrt
function sqrt(a::Posit{ps, es}) where {ps, es}
    if a.s == 1
        throw(ArgumentError("square root not defined for negative numbers."))
    end

    if isone(a)
        return a
    end
 
    temp_exp = a.e - a.fs
    needs_shift = false
    shift_amount = 0

    # Perfect square: no fraction
    exp = 2^es * a.k + a.e
    if all(c -> c == '0', a.fi) && exp ≥ 2 && iseven(exp)
        return a >> (exp ÷ 2)
    end
    
    # Adjust exponent and calculate square root of the fraction
    if iseven(a.e)
        e = a.e ÷ 2
        fs = a.fs
        if isodd(temp_exp)
            fi = sqrt_frac(shift_left("1" * a.fi, 1), fs + 1)[2:end]
        else
            fi = sqrt_frac("1" * a.fi, fs + 1)[2:end]
        end
        ers = a.ers
    else
        e = (a.e - 1) ÷ 2
        ers = max(0, min(es, ps - a.rs - 1))
        fs = ps - a.rs - ers - 1
        if isodd(temp_exp)
            fi = sqrt_frac(shift_left("1" * a.fi, 1), fs + 1)[2:end]
        else
            fi = sqrt_frac("1" * a.fi, fs + 1)[2:end]
        end
    end

    # Regime adjustment
    if iseven(a.k)
        k = div(a.k, 2)
    else
        k = div(a.k - 1, 2)
        needs_shift = true
        shift_amount = 1 << (es - 1)
    end

    rs = k ≥ 0 ? k + 2 : -k + 1
    if rs >= ps - 1
        return k ≥ 0 ? typemax(Posit{ps, es}) : zero(Posit{ps, es})
    end

    f = to_decimal(fi) / 2.0^fs
    p = Posit{ps, es}(a.s, a.sn, k, rs, ers, e, fi, fs, f, a.useed)
    return needs_shift ? p << shift_amount : p
end

import Base: exp
function exp(a::Posit{ps, es}) where {ps, es}
    if iszero(a)
        return one(Posit{ps, es})
    end

    # Negative exponent: use exp(-x) = 1 / exp(x)
    if a < zero(Posit{ps, es})
        return one(Posit{ps, es}) / exp(-a)
    end

    exp_a = one(Posit{ps, es})
    term = one(Posit{ps, es})

    n = 1
    while abs(term) > eps(exp_a)
        # term = term * a / n
        term *= a / x2Posit(n, Val(ps), Val(es))
        exp_a += term
        n += 1
    end

    return exp_a
end

import Base: log
function log(a::Posit{ps, es}) where {ps, es}
    if a.s == 1 || iszero(a)
        throw(DomainError(a, "log is only defined for positive real inputs."))
    end
    # Extract exponent as integer: k * 2^es + e
    exponent = 2^es * a.k + a.e

    # Rescale to remove the exponent, bringing y closer to 1
    k = 0
    e = 0
    y = Posit{ps, es}(a.s, a.sn, k, a.rs, a.ers, e, a.fi, a.fs, a.f, a.useed)
    # Apply transformation: log(x) = 2 * atanh((x - 1)/(x + 1))
    y = (y - one(Posit{ps, es})) / (y + one(Posit{ps, es}))

    sum = zero(Posit{ps, es})
    last_sum = copy(sum)

    # Series expansion of atanh(y)
    for k in 0:50
        term = (y ^ (2k + 1)) / x2Posit(2k + 1, Posit{ps, es})
        sum += term
        if last_sum == sum
            break
        end
        last_sum = copy(sum)
    end

    # Use log(2) in high precision
    log2 = big"0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687542001481020570"

    return x2Posit(exponent * log2, Posit{ps, es}) + (x2Posit(2.0, Posit{ps, es}) * sum)
end
import Base: log
log(a::Posit{ps, es}) where {ps, es} = log(a)  # log == log

import Base: ^
function ^(a::Posit{ps, es}, b::Integer) where {ps, es}
    if b == 0
        return one(Posit{ps, es})
    
    # a^(-b) = 1 / a^b
    elseif b < 0
        return one(Posit{ps, es}) / (a ^ (-b))
    end


    # Exponentiation by squaring
    result = one(Posit{ps, es})
    base = a
    exp = b

    while exp > 0
        if isodd(exp)
            result *= base
        end
        base *= base
        exp ÷= 2
    end

    return result
end

function ^(a::Posit{ps, es}, b::Posit{ps, es}) where {ps, es}
    # a^0 = 1
    T = Posit{ps,es}

    if iszero(a)
        return iszero(b) ? one(T) : zero(T)

    elseif iszero(b)
        return one(T)

    # a^(±1/2) = ±sqrt(a)
    elseif abs(b) == x2Posit(0.5, T)
        return b > zero(T) ? sqrt(a) : one(T) / sqrt(a)

    # a^(-b) = 1 / a^b
    elseif b < zero(Posit{ps, es})
        return one(Posit{ps, es}) / (a ^ (-b))

    else
        # Se a < 0 e b não é inteiro → erro de domínio
        if a.s == 1
            b_int = int_floor(b)
            b_is_integer = b == x2Posit(b_int, T)
            if !b_is_integer
                throw(DomainError((a, b), "Exponentiation a^b is undefined for a < 0 when b is not an integer."))
            end

            result = exp(b * log(abs(a)))
            return iseven(b_int) ? result : -result
        end

        return exp(b * log(a))
    end
end

function +(a::Posit{ps, es}, b::Integer) where {ps, es}
    b_posit = x2Posit(b, Posit{ps, es})
    return a + b_posit
end

function -(a::Posit{ps, es}, b::Integer) where {ps, es}
    b_posit = x2Posit(b, Posit{ps, es})
    return a * b_posit
end

Base.:+(a::Integer, b::Posit{ps, es}) where {ps, es} = b + a
Base.:*(a::Integer, b::Posit{ps, es}) where {ps, es} = b * a

import Base: floor
import Base: mod

function mod(a::Posit{ps, es}, b::Posit{ps, es})::Posit{ps, es} where {ps, es}
    if iszero(b)
        error("Division by zero in mod function")
    end
    # Compute a % b as: a - b * floor(a / b)
    return a - (b * floor(a / b))
end

import Base: max
function max(a::Posit{ps, es}, b::Posit{ps2, es2}) where {ps, es, ps2, es2}
    return a >= b ? a : b
end

function max(a::Posit{ps, es}, b::T) where {ps, es, T <: Real}
    b_posit = x2Posit(b, Posit{ps, es})
    return max(a, b_posit)
end

function max(a::T, b::Posit{ps, es}) where {ps, es, T <: Real}
    return max(b, a)
end
