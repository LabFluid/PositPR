function leading_ones(arr)
    count = 0
    for bit in arr
        if bit == '1'
            count += 1
        else
            break
        end
    end
    return count
end

function leading_zeros(arr)
    count = 0
    for bit in arr
        if bit == '0'
            count += 1
        else
            break
        end
    end
    return count
end

function leading_ones_inverse(count)
    return "1" ^ count * "0"
end

function map_to_bitstring(num)
    if num == 0
        return "10"
    elseif num > 0
        return "1" * "1" ^ num * "0"
    else
        return "0" ^ (-num) * "1"
    end
end

function leading_zeros_inverse(count)
    return "0" ^ count * "1"
end

function and(str1::String, str2::String)
    result = ""
    for (bit1, bit2) in zip(str1, str2)
        result *= (bit1 == '1' && bit2 == '1') ? "1" : "0"
    end
    return result
end

function or(str1::String, str2::String)
    result = ""
    for (bit1, bit2) in zip(str1, str2)
        result *= (bit1 == '1' || bit2 == '1') ? "1" : "0"
    end
    return result
end

function xor(str1::String, str2::String)
    result = ""
    for (bit1, bit2) in zip(str1, str2)
        result *= (bit1 != bit2) ? "1" : "0"
    end
    return result
end

function not(str::String)
    result = ""
    for bit in str
        result *= (bit == '1') ? "0" : "1"
    end
    return result
end

function shift_left(str::String, shift::Int)
    return str * repeat("0", shift)
end

function shift_right(str::String, shift::Int)
    return "0" ^ min(length(str), shift) * str[1:end-shift]
end

function ornot(str::String)
    for c in str
        if c == '1'
            return 0
        end
    end
    return 1
end

function to_decimal(bitstr::String)
    if length(bitstr) == 0
        return Int(0)
    end
	decimal = parse(Int, bitstr, base=2)
    return convert(Int, decimal)
end

function to_binary(decimal_value, num_bits)
    binary_str = string(decimal_value, base=2)
    
    while length(binary_str) < num_bits
        binary_str = "0$binary_str"
    end
    
    return binary_str
end

function get_x(s,k, e, f, useed) # float representation of a Posit
    return s*useed^k * 2.0^e * f
end

function twos_complement(str::String)
    result = not(str) 
	carry = true
	
    for i in reverse(1:length(result))
        if carry
            if result[i] == '1'
                result = string(result[1:i-1], '0', result[i+1:end])
            else
                result = string(result[1:i-1], '1', result[i+1:end])
                carry = false
            end
        end
    end

    return result
end

function twos_complement(p::Posit{ps, es}) where {ps, es}
    return PositDecoder(twos_complement(PositEncoder(p)), ps, es)
end

function join_BP(BP, size)
    joined = join(BP)
    
    if length(joined) >= size
        return joined[1:size]
    else
        zeros_to_add = size - length(joined)
        return joined * "0"^zeros_to_add
    end
end

function join_frac(fi, frac_size)
    joined = join(fi)
    
    if length(joined) >= frac_size
        return joined[1:frac_size]
    else
        zeros_to_add = frac_size - length(joined)
        return "0"^zeros_to_add * joined
    end
end

function sum_frac(str1::String, str2::String, e)# sum two fraction bitstrings and adjust exponent
    pad_str1 = rpad(str1, length(str2), '0')
    pad_str2 = rpad(str2, length(str1), '0')
  
    result = zeros(Int8, length(pad_str1))
  
    carry = 0
    for i in length(pad_str1):-1:1
      bit1 = pad_str1[i] - '0'
      bit2 = pad_str2[i] - '0'
  
      sum = bit1 + bit2 + carry
  
      if sum > 1
        carry = 1
        result[i] = sum - 2
      else
        carry = 0
        result[i] = sum
      end
    end
	
	sum_bin = string(join(result))
		 
  	if carry == 1
      return sum_bin, e+1
    else
      return sum_bin, e
    end
end

function sum_frac_plus(str1::String, str2::String, e)# sum two fraction bitstrings and adjust exponent (for + operation)

    pad_str1 = rpad(str1, length(str2), '0')
    pad_str2 = rpad(str2, length(str1), '0')
  
    result = zeros(Int8, length(pad_str1))
  
    carry = 0
    for i in length(pad_str1):-1:1
      bit1 = pad_str1[i] - '0'
      bit2 = pad_str2[i] - '0'
  
      sum = bit1 + bit2 + carry
  
      if sum > 1
        carry = 1
        result[i] = sum - 2
      else
        carry = 0
        result[i] = sum
      end
    end
  
    sum_bin = string(join(result))
  
    if carry == 1
        sum_bin, e+1
    else
        sum_bin[2:end], e
    end
end

function sub_frac(str1::String, str2::String, e)# subtract two fraction bitstrings and adjust exponent
    pad_str1 = rpad(str1, length(str2), '0')
    pad_str2 = rpad(str2, length(str1), '0')

    result = zeros(Int8, length(pad_str1))

    carry = 0
    for i in length(pad_str1):-1:1
        bit1 = pad_str1[i] - '0'
        bit2 = pad_str2[i] - '0'

        diff = bit1 - bit2 - carry

        if diff < 0
            carry = 1
            result[i] = diff + 2
        else
            carry = 0
            result[i] = diff
        end
    end
	
	result = string(join(result))
	
	i_first_non_zero = findfirst(x -> x != '0', result)

    if i_first_non_zero isa Nothing
        return nothing, nothing
    end
		
	sub_bin = result[i_first_non_zero+1:end]
	sub_bin, e-(i_first_non_zero-1)
end

function mul_frac(str1::String, str2::String)# multiply two fraction bitstrings and adjust exponent
      
	str1 = rpad(str1, length(str2), '0') 
	str2 = rpad(str2, length(str1), '0')


	# assign the string with more zeros to str2, fewer operations in the for loop
	if count(x -> x == '0', str1) > count(x -> x == '0', str2)
		str1, str2 = str2, str1
	end
	
	m = length(str1)

	temp = ["" for _ in 1:m]

	# multiply "1010"*"0110" = "101000" + "10100" (sum of "1010" * number of j zeros, where j is the position of bits equal to 1) 
	for j in m:-1:1
		bit2 = str2[j] - '0'
		
		if bit2 == 0
			continue
		else
			temp[j] = "0"^ j *str1
		end
	end
	
	# non empty results
	result_temp = [x for x in temp if x ≠ ""]
	
	n = length(result_temp)

	if n == 0
		return "0" ^ m
	end
	
	if n == 1
		return result_temp[1]
	end
	
	sum, e = sum_frac(result_temp[1], result_temp[2], 0)
	for i in 3:n
		sum, e = sum_frac(sum, result_temp[i], 0)
	end
	
	if e ≥ 1
		sum = "1"*sum
	end
	return sum
end

function exponent_update(e_expanded, es) 
    if e_expanded < 0
        k = div(e_expanded+1, 2^(es))-1
    else
        k = div(e_expanded, 2^(es))
    end

    if e_expanded >= 2^(es)
        e_expanded -= 2^(es)
    end

    e = mod(e_expanded, 2^(es))

    return k, e
end

# Smallest positive number (excluding zero)
function smallest(::Type{Posit{ps, es}}) where {ps, es}
    return Posit{ps, es}(0, 0, -ps + 2, ps - 1, 0, 0, "", 0, 0.0, 2.0^(2^es))
end

# Minimum representable number (most negative finite number)
import Base: typemin
function typemin(::Type{Posit{ps, es}}) where {ps, es}
    return Posit{ps, es}(1, 0, ps - 2, ps - 1, 0, 0, "", 0, 0.0, 2.0^(2^es))
end

# Maximum representable number (most positive finite number)
import Base: typemax
function typemax(::Type{Posit{ps, es}}) where {ps, es}
    return Posit{ps, es}(0, 0, ps - 2, ps - 1, 0, 0, "", 0, 0.0, 2.0^(2^es))
end

function compare_bin(bin1::Vector{Int8}, bin2::Vector{Int8})
	lengthdiff = length(bin1) - length(bin2)
	if lengthdiff > 0
		for i in 1:lengthdiff
			if bin1[i] == 1
				return 1
			end
		end
		for i in 1:length(bin2)
	        if bin1[i+lengthdiff] > bin2[i]
	            return 1			
	        elseif bin1[i+lengthdiff] < bin2[i]
	            return -1
	        end
		end
	elseif lengthdiff < 0 
		for i in 1:-lengthdiff
			if bin2[i] == 1
				return -1
			end
		end
		for i in 1:length(bin1)
	        if bin1[i] > bin2[i+lengthdiff]
	            return 1			
	        elseif bin1[i] < bin2[i+lengthdiff]
	            return -1
	        end
		end
	else
		for i in 1:length(bin1)
	        if bin1[i] > bin2[i]
	            return 1			
	        elseif bin1[i] < bin2[i]
	            return -1
	        end
		end
	end
    return 0
end


function subtract_bin(bin1::Vector{Int8}, bin2::Vector{Int8})# subtract two fraction bitstrings::Vector{Int8} and adjust exponent
	bin2 = Base.copy(bin2)::Vector{Int8}
	while length(bin1) > length(bin2)
		insert!(bin2, 1, Int8(0))
	end
    result = zeros(Int8, length(bin1))
    carry = 0
    for i in length(bin1):-1:1
        diff = bin1[i] - bin2[i] - carry
        if diff < 0
            carry = 1
            result[i] = diff + 2
        else
            carry = 0
            result[i] = diff
        end
    end
	for i in 1:length(result)
		if result[i] == 1
			return result[i:end]
		end
	end
    return Vector{Int8}()
end

function div_frac(str1::String, str2::String, qs::Integer)# divide two fraction bitstrings and adjust exponent
    pad_str1 = rpad(str1, length(str2), '0')

    quotient = zeros(Int8, qs)
    dividend = Int8.([c - '0' for c in pad_str1])
    divisor = Int8.([c - '0' for c in str2])

	function nextbit(dividend, i)::Int8
		return i <= length(dividend) ? dividend[i] : Int8(0)
	end

	remainder = Vector{Int8}(dividend[1:length(str2)])
	for i in 1:qs #qs = quotient real size after exponent_update 			
        if compare_bin(remainder, divisor) >= 0
			remainder = subtract_bin(remainder, divisor)
            quotient[i] = Int8(1)
        else
            quotient[i] = Int8(0)
		end
		if (isempty(remainder) && length(str2) + i > length(dividend))
			break
		end
		push!(remainder, nextbit(dividend,length(str2)+i))
    end
	
    quotient_bin = join(quotient)

    return quotient_bin
end

function expand_exponent(a::Posit{ps, es}) where {ps, es}
    return 2^es * a.k + a.e
end

function sqrt_frac(bin_str::String, fs)
	#println(bin_str, " ", fs)
	#https://www.reddit.com/r/math/comments/tc7lur/computing_square_roots_in_binary_by_hand_is/?rdt=55060
    bin_str = lstrip(bin_str, '0')
	
    if isempty(bin_str)
        return "0"
    end

    if length(bin_str) % 2 != 0
        bin_str = "0" * bin_str
    end

    n = length(bin_str)
    result = ""
    remainder = ""

	function get_next_bit_pair()
		if length(bin_str) > 1
			bits = bin_str[1:2]
			bin_str = bin_str[3:end]
			return bits
		end
		return "00"
	end
	
    for i in 1:fs
        remainder *= get_next_bit_pair()
        subtract_value = parse(Int, result * "01", base=2)

        if subtract_value <= parse(Int, remainder, base=2)
            result *= "1"
            remainder = string(parse(Int, remainder, base=2) - subtract_value, base=2)
        else
            result *= "0"
        end
    end

    return result
end

function Base.copy(a::Posit{ps, es}) where {ps, es}
    return Posit{ps, es}(
        a.s, a.sn, a.k, a.rs, a.ers,
        a.e, a.fi, a.fs, a.f, a.useed
    )
end

function factorize(n::Integer)
    factors = Dict{Integer, Integer}()
    d = 2
    while n > 1
        count = 0
        while n % d == 0
            n ÷= d
            count += 1
        end
        if count > 0
            factors[d] = count
        end
        d += 1
        if d*d > n
            if n > 1
                factors[n] = 1
            end
            break
        end
    end
    return factors
end

import Base: decompose

function decompose(x::Posit{ps, es}):: Tuple{Int64, Int64, Int64} where {ps, es} 
    if iszero(x)
        return 0, 0, 0
    end

    if x.sn == 1
        return ifelse(x.s == 0, 1, -1), 0, 0
    end

    num = parse(UInt64, "1" * x.fi, base=2) % Int64
    pow = expand_exponent(x) - x.fs
    den = ifelse(x.s == 0, 1, -1)
    
    return num, pow, den
end