# Experimental kmers

This is a proof-of-concept package. It is intended show a working example of a more general, very "hackable", yet still fast, kmer type for use in BioSequences.

Extra features compared to existing kmer types:
* Supports arbitrary alphabets
* Variable-length encoding from 8- to 1024 bits.

## Demonstration

The type is parameterized by three parameters: `Kmer{K, A, T}`:
* `K` is the value of k
* `A` is the alphabet type, a subtype of `Alphabet`
* `T` is the storage type, a subtype of `Unsigned`

__Arbitrary alphabets__

```julia
julia> m = Kmer{4, CharAlphabet, UInt128}("读写汉字")
4-mer of CharAlphabet:
读写汉字

julia> reverse(m)
4-mer of CharAlphabet:
字汉写读
```

__Supports large kmers through BitIntegers.jl__

```julia
julia> Kmer(RNAAlphabet{2}, randrnaseq(500))
500-mer of RNAAlphabet{2}:
CAUAUGAUGGAUGGGUUUGGUGCGCAGACUUUAGGACUA…GGCAGUAGAUAAAUAUUCAACGGAGUGUCUAUAGCUGUG

julia> kmer"GWYFPPNML"aa
9-mer of AminoAcidAlphabet:
GWYFPPNML
```

__Simple, efficient kmer iterator of any alphabet__
```
julia>  it = SimpleKmerIterator{AAKmer{5}}(randaaseq(10^6))
Main.t.SimpleKmerIterator{Main.t.Kmer{5,AminoAcidAlphabet,UInt64},LongSequence{AminoAcidAlphabet}}(DDHMGFAHPCAYFVNMFPTCVFSAVSSTHNVKLYFTQTY…MGRIPSPPWCWHIIKDQGDIFWNHWLDSCCCRKDYTIRL)

julia> @btime sum(i.data for i in it)
  875.581 μs (2 allocations: 32 bytes)
0x0091802f3c41a9b9
```

__Uses the high-level BioSequences v2 API internally, with sane fallbacks__
This allows very generic code, with following features being example of emergent properties:

* Creating kmers from kmers

```julia
julia> T1 = t.Kmer{3, RNAAlphabet{4}, UInt32};

julia> mer = t.kmer"TAGTCGCGAGAA"
12-mer of DNAAlphabet{2}:
TAGTCGCGAGAA

julia> it = t.StandardKmerIterator{T1, typeof(mer)}(mer);

julia> [canonical(i) for i in it]
10-element Array{Main.t.Kmer{3,RNAAlphabet{4},UInt32},1}:
 CUA
 ACU
 GAC
 CGA
 CGC
 CGC
 CGA
 CUC
 AGA
 GAA
 ```

 * Many basic properties like `getindex` just uses the `BioSequence` fallback
