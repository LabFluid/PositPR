function PositDecoder(BP::String, ::Val{ps}, ::Val{es}) where {ps, es}
    useed = 2.0^(2^es)
    sign = BP[1]
    sn = ornot(BP[2:end])

    if sn != 0
        return Posit{ps, es}(sign == '1' ? 1 : 0, sn, 0, 0, 0, 0, "0", 0, 0.0, useed)
    end

    if sign == '1'
        BP = twos_complement(BP)
        s = 1
        fi, fs, f = "0", 0, 0.0
    else
        s = 0
    end

    # Regime extraction
    ri = BP[2]
    if ri == '1'
        rn = leading_ones(BP[2:end])
        k = rn - 1
    else
        rn = leading_zeros(BP[2:end])
        k = -rn
    end
    rs = rn + 1

    # Exponent extraction
    ers = max(0, min(es, ps - rs - 1))

    if ers == 0
        e = 0
    else
        if rs == ps
            e = 0
        else
            e = to_decimal(BP[(rs + 2):(rs + 1 + ers)])
        end
    end

    # Fraction extraction
    fi = BP[(1 + rs + ers + 1):end]
    fs = length(fi)
    if fs == 0
        fi, f = "0", 0.0
    else
        f = to_decimal(fi) / 2.0^fs
    end

    return Posit{ps, es}(s, sn, k, rs, ers, e, fi, fs, f, useed)
end

PositDecoder(BP::String, ps::Integer, es::Integer) = PositDecoder(BP, Val(ps), Val(es))

function PositEncoder(p::Posit{ps, es}) where {ps, es}
    s, sn, k, rs, ers, e, fi, fs, f, useed = p.s, p.sn, p.k, p.rs, p.ers, p.e, p.fi, p.fs, p.f, p.useed
    BP = ["" for _ in 1:4]

    # Sign
    BP[1] = s == 0 ? "0" : "1"

    # Special cases
    if sn == 1
        BP[2] = "0" ^ (ps - 1)
        return join_BP(BP, ps)
    else
        # Regime
        if k >= (ps - 2)
            BP[2] = "1" ^ (ps - 1)
        elseif k <= -(ps - 2)
            BP[2] = leading_zeros_inverse(ps - 2)
        else
            BP[2] = map_to_bitstring(k)
        end

        # Exponent
        if e >= 1
            BP[3] = to_binary(e, ers)
        else
            BP[3] = "0" ^ ers
        end

        # Fraction
        BP[4] = f != 0 ? fi : ""

        if s == 1
            return "1" * twos_complement(join_BP(BP, ps)[2:end])
        else
            return join_BP(BP, ps)
        end
    end
end