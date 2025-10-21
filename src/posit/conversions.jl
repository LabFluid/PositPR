function x2Posit(x::Real, ::Val{ps}, ::Val{es}) where {ps, es}
	y = abs(x)
	useed = 2.0^(2^es)

	# valores padrão
	s, sn, k, rs, ers, e, fi, fs, f = 0, 0, 0, 0, 0, 0, "0", 0, 0.0

	# zero
	if x == 0.0
		sn, f = 1, 0.0
		return Posit{ps, es}(s, sn, k, rs, ers, e, fi, fs, f, useed)
	end

	# infinito
	if y == Inf
		sn, s = 1, 1
		return Posit{ps, es}(s, sn, k, rs, ers, e, fi, fs, f, useed)
	end

	s = ifelse(x < 0, 1, 0)

	# cálculo do regime
	if y >= 1
		ri = 1
		while y ≥ useed && ri ≤ ps
			k += 1
			y /= useed
			ri += 1
		end
		rs = ri + 1
	else
		ri = 0
		while y < 1 && ri ≤ ps
			k -= 1
			y *= useed
			ri += 1
		end
		rs = ri + 1
	end

	rs = min(rs, ps - 1)

	# cálculo do expoente
	if ri ≤ ps
		e = 0
		while y ≥ 2
			y /= 2
			e += 1
		end
	end

	ers = max(0, min(es, ps - rs - 1))

	if ers == 0
		e = 0
	end

	fs = ps - rs - ers - 1

	if fs != 0
		f_raw = y - 1
		frac = Int(round(f_raw * 2. ^ fs))
        if frac == 2^fs
            fi = "0" ^ fs
            f = 0.0
            e += 1
            if e == 2^ers
                k += 1
                e = 0
            end
        else
            fi = join_frac(to_binary(frac, fs), fs)
            f = to_decimal(fi) / 2^fs
        end
	else
		fi = ""
		f = 0.0
	end

	return Posit{ps, es}(s, sn, k, rs, ers, e, fi, fs, f, useed)
end

x2Posit(x::Real, ::Type{Posit{ps, es}}) where {ps, es} = x2Posit(x, Val(ps), Val(es))

function Posit2x(p::Posit) # float representation of a Posit
	if (p.sn != 0)
		return p.s == 0. ? 0. : Inf
	end
	s,k, e, f, useed = p.s, p.k, p.e, p.f, p.useed
    return (-1)^s*Float64(useed)^k * 2.0^e * (1+f)
end


import Base: convert
import Base: promote_rule
promote_rule(::Type{Posit{ps, es}}, ::Type{Bool}) where {ps, es} = Posit{ps, es}
promote_rule(::Type{Bool}, ::Type{Posit{ps, es}}) where {ps, es} = Posit{ps, es}
promote_rule(::Type{Posit{ps, es}}, ::Type{T}) where {ps, es, T<:AbstractFloat} = Posit{ps, es}
promote_rule(::Type{T}, ::Type{Posit{ps, es}}) where {ps, es, T<:AbstractFloat} = Posit{ps, es}
promote_rule(::Type{Posit{ps, es}}, ::Type{T}) where {ps, es, T<:Integer} = Posit{ps, es}
promote_rule(::Type{T}, ::Type{Posit{ps, es}}) where {ps, es, T<:Integer} = Posit{ps, es}
promote_rule(::Type{T}, ::Type{PositPR.Posit{ps, es}}) where {ps, es, T<:AbstractIrrational} = Posit{ps, es}
promote_rule(::Type{PositPR.Posit{ps, es}}, ::Type{T}) where {ps, es, T<:AbstractIrrational} = Posit{ps, es}

convert(::Type{Posit{ps, es}}, x::Real) where {ps, es} = x2Posit(x,Posit{ps,es})
convert(::Type{Posit{ps, es}}, x::Posit{ps, es}) where {ps, es} = x
convert(::Type{Posit{ps, es}}, x::Bool) where {ps, es} = x ? one(Posit{ps, es}) : zero(Posit{ps, es})

Posit{ps,es}(x::Real) where {ps,es} = x2Posit(x,Posit{ps,es})

Posit{ps,es}(x::Posit{ps,es}) where {ps,es} = copy(x)

function Posit2Posit(p0::Posit{ps0, es0}, ::Val{ps_new}, ::Val{es_new}) where {ps0, es0, ps_new, es_new}
    useed_new = 2.0^(2^es_new)

    # Special case: zero or infinity (sn = 1)
    if p0.sn == 1
        return Posit{ps_new, es_new}(p0.s, 1, 0, 0, 0, 0, "0", 0, 0.0, useed_new)
    end

    # Initialize fields for the new Posit
    k = 0
    rs = 2
    e = 0
    ers = 0
    fs = 0
    fi = "0"
    f = 0.0

    if es_new == es0
        # Same exponent size, only adjusting precision
        k, rs, ers, e = p0.k, p0.rs, p0.ers, p0.e

    elseif es_new > es0
        # More exponent bits: fewer regime bits needed

        temp_regime = 2^(es_new - es0)  # how many old regimes fit in a new regime

        k = p0.k ÷ temp_regime  # adjust regime index

        if p0.k < 0  # for numbers in [0,1] interval
            k = (p0.k + 1) ÷ temp_regime - 1
        end

        # recover exponent bits "folded" into previous regime
        e = p0.e + mod(p0.k, temp_regime) * 2^es0

    else
        # Fewer exponent bits: more bits absorbed by regime

        temp_regime = 2^(es0 - es_new)

        k = p0.k * temp_regime + p0.e ÷ (2^es_new)

        # remove exponent bits that overflowed into the regime
        e = p0.e - mod(k, temp_regime) * 2^es_new
    end

    # Regime size
    rs = k ≥ 0 ? k + 2 : -k + 1
    rs = min(rs, ps_new - 1)  # clamp to max allowed size

    # Exponent size
    ers = max(0, min(es_new, ps_new - rs - 1))

    # Fraction size
    fs = ps_new - rs - ers - 1

    # Fraction value
    if fs > p0.fs
        fi = p0.fi * leading_zeros_inverse(fs - length(p0.fi))[1:end-1]
    else
        fi = p0.fi[1:fs]
    end

    if fs != 0
        f = to_decimal(fi) / 2.0^fs
    else
        fi = "0"
        f = 0.0
    end

    return Posit{ps_new, es_new}(p0.s, 0, k, rs, ers, e, fi, fs, f, useed_new)
end

Posit2Posit(p::Posit{ps0, es0}, ps_new::Integer, es_new::Integer) where {ps0, es0} = Posit2Posit(p, Val(ps_new), Val(es_new))

function PositConverter(p1::Posit{ps1, es1}, p2::Posit{ps2, es2}) where {ps1, es1, ps2, es2}
    if ps1 > ps2
        # Convert p2 to the type of p1
        p2 = Posit2Posit(p2, Val(ps1), Val(es1))
    elseif ps1 < ps2 || es1 != es2
        # Convert p1 to the type of p2
        p1 = Posit2Posit(p1, Val(ps2), Val(es2))
    end

    return p1, p2
end

import Base: Int

function Int(a::Posit{ps, es}) where {ps, es}
    if iszero(a)
        return 0
    end

    e_expanded = expand_exponent(a)

    if e_expanded < -1
        return 0
    end

    n_bits = min(e_expanded, length(a.fi))

    if n_bits > 0
        bits = rpad(a.fi[1:n_bits], e_expanded, '0')
        abs_value = parse(Int, bits, base=2) + 2^e_expanded
    else 
        abs_value = 1
    end

    if e_expanded > -1 && length(a.fi) > e_expanded && a.fi[e_expanded + 1] == '1'
        abs_value += 1
    end

    return ifelse(a.s == 1, -abs_value, abs_value)
end

import Base: Integer
Base.Integer(a::Posit{ps,es}) where {ps,es} = Int(a)

convert(::Type{Int}, a::Posit{ps, es}) where {ps, es} = Int(a)

function int_floor(a::Posit{ps, es}) where {ps, es}
    if iszero(a)
        return 0
    end

    e_expanded = expand_exponent(a)

    if e_expanded <= -1
        return ifelse(a.s == 1, -1, 0) # -1 < a < 1
    end

    n_bits = min(e_expanded, length(a.fi))
    abs_value = 2^e_expanded

    if n_bits > 0
        bits = rpad(a.fi[1:n_bits], e_expanded, '0')
        abs_value += parse(Int, bits, base=2)
    end

    has_fraction = ('1' in a.fi[n_bits+1:end])
    return a.s == 0 ? abs_value : -abs_value - (has_fraction ? 1 : 0)
end

import Base.floor

function floor(a::Posit{ps, es}) where {ps, es}
    x2Posit(int_floor(a), typeof(a))
end

import Base: real
Base.real(p::Posit{ps,es}) where {ps,es} = p

import Base: Float64
Float64(p::Posit) = Posit2x(p)

import Base: float
float(p::Posit) = Posit2x(p)
