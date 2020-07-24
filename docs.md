## Constructing kmers

Kmer can be constructed the following ways:

* Constructing a specified `Kmer` type from an iterable:
```
julia> Kmer{5, DNAAlphabet{4}, UInt32}("TACGC")
5-mer of DNAAlphabet{4}:
TACGC
```
* Constructing a `Kmer` from an alphabet and a sequence. Note that this is type unstable since `K` and the data parameter is determined at runtime.
```
julia> Kmer(DNAAlphabet{4}(), "TACGC")
5-mer of DNAAlphabet{4}:
TACGC
```
* Constructing using the `kmer" ... "` string literal - also type unstable:
```
julia> kmer"TAGC"
4-mer of DNAAlphabet{2}:
TAGC
```
The alphabet type can be specified using the string literal by appending one of the following flags to the literal: "aa" or "a" for `AminoAcidAlphabet`, "rna" or "r" for `RNAAlphabet{2}`, and "dna" or "d" for `DNAAlphabet{2}`. No other alphabets are supported at this time. Note that in the absence of a flag, it default to DNA.
```
julia> kmer"AUCCGC"r
6-mer of RNAAlphabet{2}:
AUCCGC
```
* Finally, you can construct them directly from the underlying data. This is efficient, but not recommended for user-facing code, because no input validation will be done:
```
julia> Kmer{1, AminoAcidAlphabet, UInt8}(0xff)
1-mer of AminoAcidAlphabet:
Error showing value of type Kmer{1,AminoAcidAlphabet,UInt8}:
ERROR: cannot decode 255 in AminoAcidAlphabet
```

## Operations on kmers
Kmers are subtypes of `BioSequence`, and so many biosequence-related methods are defined for them. However, they are immutable. The defined operations include:

* getindex
* length
* reverse
* complement (for those with an alphabet subtyped from `NucleicAcidAlphabet`)
* reverse_complement
* iterate
* canonical
* rand
* Random.shuffle

## Kmer iterators
Kmers are ususally extracted from some other sequence. Of course, kmers may be manually constructed from other sequences:

```
julia> seq = randdnaseq(10)
10nt DNA Sequence:
TATACCCTCT

julia> for i in 1:3
       println(t.Kmer(DNAAlphabet{2}(), seq[i:i+4]))
       end
TATAC
ATACC
TACCC
```

But this is not very efficient. It is faster to use the predefined kmer iterators. The simplest one is `SimpleKmerIterator`. This takes a kmer type as a parameter, and requires that the underlying sequence it iterates over is of the same alphabet as the kmers:
```
julia>  seq = randseq(DNAAlphabet{2}(), 10)
10nt DNA Sequence:
GCGGCAC

julia> it = SimpleKmerIterator{DNAKmer{5}}(seq);

julia> collect(it)
6-element Array{Kmer{5,DNAAlphabet{2},UInt64},1}:
 GCGGC
 CGGCA
 GGCAC
 ...
```

Another iterator is the `StandardKmerIterator`. Iterating over this does not yield kmers, but rather a struct containing the kmer, the reverse-complemented kmer, and the starting position of the kmer in the sequence. For this reason, it can only work with nucleotide sequences and kmers, but the sequence and kmer alphabet types can differ. By convention, all noncanonical kmers are skipped during iteration.
```
julia> it = StandardKmerIterator{DNAKmer{5}}(seq);

julia> collect(it)
6-element Array{MerIter{Kmer{5,DNAAlphabet{2},UInt64}},1}:
 MerIter{Kmer{5,DNAAlphabet{2},UInt64}}(1, GCGGC, GCCGC)
 MerIter{Kmer{5,DNAAlphabet{2},UInt64}}(2, CGGCA, TGCCG)
 MerIter{Kmer{5,DNAAlphabet{2},UInt64}}(3, GGCAC, GTGCC)
 ...
```
The canonical kmer can be obtained from a `KmerIterResult` with `canonical(i)`.
