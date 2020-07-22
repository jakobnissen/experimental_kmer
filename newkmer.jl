module t

using BitIntegers
using BioSequences

import Random

import BioSequences: encoded_data, extract_encoded_element, complement_bitpar,
reversebits, BitsPerSymbol, reverse, complement, reverse_complement,
MerIterResult, encoded_data_eltype, repeatpattern, canonical, inbounds_getindex,
encoded_data_type, bits_per_symbol, bitmask, decode, encode, bitindex, BitIndex

struct Kmer{K, A <: Alphabet, T <: Unsigned} <: BioSequence{A}
    data::T

    function Kmer{K, A, T}(data::Integer) where {K, A, T}
        if capacity(Kmer{K, A, T}) < K
            throw(ArgumentError("Kmer of type $T cannot support a K of $K"))
        end
        return new(data)
    end
end

const RNAKmer{K} = Kmer{K, RNAAlphabet{2}, UInt64}
const DNAKmer{K} = Kmer{K, DNAAlphabet{2}, UInt64}

function Base.convert(::Type{<:Kmer{K, A, T1}}, m::Kmer{K, A, T2}) where {K, A, T1, T2}
    return Kmer{K, A, T1}(m.data % T1)
end

# TODO: Add convertsion methods from NucleicAcidAlphabet to each other

encoded_data(m::Kmer) = m.data
Base.length(m::Kmer{K}) where K = K
mask(::Type{<:Kmer{K, A, T}}) where {K, A, T} = one(T) << (bits_per_symbol(A()) * K) - 1
capacity(::Type{<:Kmer{K, A, T}}) where {K, A, T} = div(8, bits_per_symbol(A())) * sizeof(T)
ksize(::Type{<:Kmer{K}}) where K = K
Base.summary(x::Kmer{K, A}) where {K, A} = string(K, "-mer of ", A)
encoded_data_type(::Type{<:Kmer{K, A, T}}) where {K, A, T} = T
Base.zero(::Type{T}) where {T <: Kmer} = T(zero(encoded_data_type(T)))

# TODO: Add this to Twiddle?
↞(data::Unsigned, shift::Integer) = data << (shift & (8*sizeof(data) - 1))
↠(data::Unsigned, shift::Integer) = data >>> (shift & (8*sizeof(data) - 1))

@inline function bitindex(m::Kmer, i::Integer)
    B = bits_per_symbol(m)
    T = encoded_data_type(m)
    val = (B * (ksize(typeof(m))-i))
    return BitIndex{B, encoded_data_type(m)}(val)
end

function extract_encoded_element(bitind::BitIndex, data::Unsigned)
    return (data ↠ bitind.val) & bitmask(bits_per_symbol(bitind))
end

# TODO: Fix the one in BioSequences instead or lift it here.
@inline function reversebits(x::Unsigned, ::BitsPerSymbol{2})
     mask = repeatpattern(typeof(x), 0x33)
     x = ((x >> 2) & mask) | ((x & mask) << 2)
     return reversebits(x, BitsPerSymbol{4}())
end

@inline function reversebits(x::Unsigned, ::BitsPerSymbol{4})
     mask = repeatpattern(typeof(x), 0x0F)
     x = ((x >> 4) & mask) | ((x & mask) << 4)
     return bswap(x)
end

reversebits(x::Unsigned, ::BitsPerSymbol{8}) = bswap(x)

@inline function complement(m::Kmer{K, <:NucleicAcidAlphabet}) where K
    typeof(m)(complement_bitpar(m.data, Alphabet(m)) & mask(typeof(m)))
end

@inline function reverse(m::Kmer{K, A}) where {K, A}
    T = typeof(m)
    data = reversebits(m.data, BitsPerSymbol(m))
    offset = bits_per_symbol(A()) * (capacity(T) - ksize(T))
    return typeof(m)(data ↠ offset)
end

reverse_complement(m::Kmer) = reverse(complement(m))

function canonical(m::Kmer)
    rv = reverse_complement(m)
    return ifelse(m < rv, m, rv)
end

# Add symbol to end, pushing position in kmer
@inline function push(m::Kmer, symbol_)
    symbol = convert(eltype(typeof(m)), symbol_)
    encoding = encode(Alphabet(m), symbol)
    return push(m, encoding)
end

@inline function push(m::Kmer, encoded::Unsigned)
    data = m.data ↞ bits_per_symbol(m)
    data |= encoded
    data &= mask(typeof(m))
    return typeof(m)(data)
end

@inline function pushfirst(m::Kmer, x)
    shifted = typeof(m)(m.data ↠ bits_per_symbol(m))
    return or_bits(shifted, x, 1)
end

@inline function or_bits(m::Kmer, symbol_, i::Integer)
    symbol = convert(eltype(typeof(m)), symbol_)
    encoding = encode(Alphabet(m), symbol)
    return or_bits(m, encoding, i)
end

@inline function or_bits(m::Kmer, encoding::Unsigned, i::Integer)
    add = (encoding % encoded_data_type(m)) ↞ bitindex(m, i).val
    return typeof(m)(m.data | add)
end

Base.cmp(x::T, y::T) where {T<:Kmer} = cmp(encoded_data(x), encoded_data(y))
Base.:(==)(x::T, y::T) where {T<:Kmer} = encoded_data(x) == encoded_data(y)
Base.isless(x::T, y::T) where {T<:Kmer} = isless(encoded_data(x), encoded_data(y))

# TODO: Figure out how to do this effectively for all alphabets
function Base.hash(x::Mer{<:NucleicAcidAlphabet{2},K}, h::UInt) where {K}
    return Base.hash_uint64(encoded_data(x) ⊻ K ⊻ h)
end

function Kmer{K, A, T}(s) where {K, A, T}
    kT = Kmer{K, A, T}
    m = zero(kT)
    k = 0
    for i in s
        k += 1
        k > capacity(kT) && throw(InexactError(:Kmer, kT, s))
        m = or_bits(m, i, k)
    end
    k == K || throw(ArgumentError("Input sequence not length K"))
    return m
end

# TODO: Add this to Twiddle?
# This code is weird, but it's the only way I can get it to constant fold
# for UInt1024.
function repeatpattern(::Type{T}, pattern::P) where {T <: Unsigned, P <: Unsigned}
    y = Ref(zero(T))
    GC.@preserve y begin
        ptr = Ptr{P}(pointer_from_objref(y))
        for i in 1:div(sizeof(T), sizeof(P))
            unsafe_store!(ptr, one(P), i)
        end
    end
    return y[] * pattern
end

# TODO: This only works for up to 8 bits per symbol
# It allocates, but in-register shuffles are very slow for larger kmers
function Random.shuffle(m::Kmer)
    # Convert to vector
    v = Vector{UInt8}(undef, length(m))
    data = m.data
    @inbounds for i in eachindex(m)
        v[i] = (data & bitmask(bits_per_symbol(m))) % eltype(v)
        data >>>= bits_per_symbol(m)
    end

    # Fisher-Yates shuffle
    @inbounds for i in 2:length(m)
        j = rand(1:i)
        v[i], v[j] = v[j], v[i]
    end

    # Convert back to kmer
    result = zero(typeof(m))
    @inbounds for i in eachindex(m)
        result = or_bits(result, v[i], i)
    end

    return result
end

###########################

# Generic fallback for biosequence type
@inline extract_data(s::BioSequence, i::Integer) = @inbounds s[i]

# This is more efficient since we can load data directly from one encoding to
# the other with no re-encoding/decoding
@inline function extract_data(s::Union{Kmer, LongSequence}, i::Integer)
    extract_encoded_element(bitindex(s, i), encoded_data(s))
end

struct SimpleKmerIterator{K <: Kmer, S <: BioSequence}
    seq::S

    function SimpleKmerIterator{K, S}(s::S) where {K<:Kmer{<:Any, A} where A, S}
        if Alphabet(K) !== Alphabet(S)
            throw(ArgumentError("Parameters K and S must have same alphabets"))
        end
        if capacity(K) < ksize(K)
            T = encoded_data_type(K)
            k = ksize(K)
            throw(ArgumentError("Kmer of type $T cannot support a K of $k"))
        end
        new(s)
    end
end

SimpleKmerIterator{K}(s::BioSequence) where K = SimpleKmerIterator{K, typeof(s)}(s)

Base.length(s::SimpleKmerIterator{K}) where K = length(s.seq) - ksize(K) + 1
Base.eltype(s::Type{<:SimpleKmerIterator{K}}) where K = K

function Base.iterate(it::SimpleKmerIterator{K}) where K
    length(it) < 1 && return nothing
    m = zero(K)
    for i in eachindex(m)
        data = extract_data(it.seq, i)
        m = or_bits(m, data, i)
    end
    return m, (m, ksize(K) + 1)
end

function Base.iterate(it::SimpleKmerIterator{K}, state) where K
    m, i = state
    i > length(it.seq) && return nothing
    data = extract_data(it.seq, i)
    m = push(m, data)
    return m, (m, i + 1)
end

############################

struct MerIter{M <: Kmer}
    pos::Int
    fw::M
    rv::M
end

canonical(m::MerIter) = ifelse(m.fw < m.rv, m.fw, m.rv)

struct StandardKmerIterator{K <: Kmer{<:Any, <:NucleicAcidAlphabet}, S <: BioSequence{<:NucleicAcidAlphabet}}
    seq::S

    function StandardKmerIterator{K, S}(s::S) where {K<:Kmer{<:Any, A} where A, S}
        if capacity(K) < ksize(K)
            T = encoded_data_type(K)
            k = ksize(K)
            throw(ArgumentError("Kmer of type $T cannot support a K of $k"))
        end
        new(s)
    end
end

StandardKmerIterator{K}(s::BioSequence) where K = StandardKmerIterator{K, typeof(s)}(s)
Base.IteratorSize(::Type{<:StandardKmerIterator{<:Any, <:BioSequence{<:NucleicAcidAlphabet{4}}}}) = Base.SizeUnknown()
Base.length(::StandardKmerIterator{K, BioSequence{<:NucleicAcidAlphabet{2}}}) where K = length(s.seq) - ksize(K) + 1
Base.eltype(s::Type{<:StandardKmerIterator{K}}) where K = MerIter{K}

complement_encoded(A::NucleicAcidAlphabet{2}, bits::Unsigned) = bits ⊻ typeof(bits)(3)
complement_encoded(A::NucleicAcidAlphabet{4}, bits::Unsigned) = @inbounds FOURBIT_C_LUT[bits + 1]

const TWOBIT_LUT = (0xf, 0x0, 0x1, 0x0,
                    0x2, 0xf, 0xf, 0xf,
                    0x3, 0xf, 0xf, 0xf,
                    0xf, 0xf, 0xf, 0xf)

const FOURBIT_LUT = (0x1, 0x2, 0x4, 0x8)

const FOURBIT_C_LUT = (0x00, 0x08, 0x04, 0x0c,
                       0x02, 0x0a, 0x06, 0x0e,
                       0x01, 0x09, 0x05, 0x0d,
                       0x03, 0x0b, 0x07, 0x0f)

# Fallback
convert_encoding(to::NucleicAcidAlphabet{N}, from::NucleicAcidAlphabet{N}, x) where N = x
convert_encoding(to::NucleicAcidAlphabet{2}, from::NucleicAcidAlphabet{4}, x) = @inbounds TWOBIT_LUT[x + 1]
convert_encoding(to::NucleicAcidAlphabet{4}, from::NucleicAcidAlphabet{2}, x) = @inbounds FOURBIT_LUT[x + 1]

function fill_k!(::Type{K}, seq, p) where {K <: Kmer}
    fw = zero(K)
    rv = zero(K)
    good = true
    for i in 1:ksize(K)
        data = extract_data(seq, p+i-1)
        encoding = convert_encoding(Alphabet(K), Alphabet(seq), data)
        complemented = complement_encoded(Alphabet(K), encoding)
        fw = or_bits(fw, encoding, i)
        rv = or_bits(rv, complemented, ksize(K)-i+1)
        good &= encoding < 0xf
    end
    return fw, rv, good
end

@inline function Base.iterate(it::StandardKmerIterator{K}) where K
    good = false
    len = length(it.seq)
    i = 1
    fw = rv = zero(K)
    while !good
        i + ksize(K) > len && return nothing
        fw, rv, good = fill_k!(K, it.seq, i)
        i += ksize(K)
    end
    res = MerIter(i, fw, rv)
    return res, res
end

@inline function Base.iterate(it::StandardKmerIterator{K}, state) where K
    i = state.pos + 1
    fw = state.fw
    rv = state.rv
    len = length(it.seq)
    i > len && return nothing
    data = extract_data(it.seq, i)
    encoding = convert_encoding(Alphabet(K), Alphabet(it.seq), data)
    complemented = complement_encoded(Alphabet(K), encoding)
    fw = push(fw, encoding)
    rv = pushfirst(rv, complemented)
    good = encoding < 0xf
    while !good
        i + ksize(K) > len && return nothing
        fw, rv, good = fill_k!(K, it.seq, i)
        i += ksize(K)
    end
    res = MerIter(i, fw, rv)
    return res, res
end

end # t
