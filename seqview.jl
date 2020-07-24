# Mention in docs no boundscheck is performed on instantiation
struct SeqView{A<:Alphabet} <: BioSequence{A}
    data::Vector{UInt64}
    part::UnitRange{Int}
end

Base.length(v::SeqView) = last(v.part) - first(v.part) + 1
encoded_data(v::SeqView) = v.data

function SeqView(seq::LongSequence{A}, part::UnitRange{Int}) where A
    @boundscheck checkbounds(seq, part)
    return SeqView{A}(seq.data, part)
end

@inline function bitindex(seq::SeqView, i::Integer)
    return bitindex(BitsPerSymbol(seq), encoded_data_eltype(seq), i + first(seq.part) - 1)
end

function Base.setindex!(seq::SeqView, x, i::Integer)
    @boundscheck checkbounds(seq, i)
    return unsafe_setindex!(seq, x, i)
end

@inline function unsafe_setindex!(seq::SeqView{A}, x, i::Integer) where {A}
    bin = enc64(seq, x)
    return encoded_setindex!(seq, bin, i)
end

function enc64(::SeqView{A}, x) where {A}
    #TODO: Resolve these two use cases of A().
    return UInt64(encode(A(), convert(eltype(A()), x)))
end

@inline function encoded_setindex!(seq::SeqView{A}, bin::UInt64, i::Integer) where {A}
    j, r = bitindex(seq, i)
    data = encoded_data(seq)
    @inbounds data[j] = (bin << r) | (data[j] & ~(bindata_mask(seq) << r))
    return seq
end

# This is also already implemented for LongSeq, I think
function reverse!(s::SeqView)
	i, j = 1, lastindex(s)
	@inbounds while i < j
		s[i], s[j] = s[j], s[i]
		i, j = i + 1, j - 1
	end
	return s
end

function complement!(s::SeqView{A}) where {A <: NucleicAcidAlphabet}
	bps = bits_per_symbol(A())
	bi = firstbitindex(s)
	i = 1
	stop = lastbitindex(s) + bps
	@inbounds while !iszero(offset(bi))
		s[i] = complement(s[i])
		bi += bps
		i += 1
	end
	@inbounds for j in index(bi):index(stop)-1
		s.data[j] = complement_bitpar(s.data[j], Alphabet(s))
		bi += 64
		i += symbols_per_data_element(s)
	end
	@inbounds while bi < stop
		s[i] = complement(s[i])
		bi += bps
		i += 1
	end
	return s
end
