
## The ORFI type

For convenience, the `ORFI` type is more stringent in preventing the creation of incompatible instances. As a result, attempting to create an instance with incompatible parameters will result in an error. For instance, the following code snippet will trigger an error:

```julia
ORFI{4,NaiveFinder}(1:10, '+', 4) # Or any F <: GeneFinderMethod

ERROR: MethodError: no method matching OpenReadingFrameInterval{4, NaiveFinder}(::UnitRange{Int64}, ::Char, ::Int64)

Closest candidates are:
  (::Type{OpenReadingFrameInterval{N, F}} where {N, F})(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:49
  OpenReadingFrameInterval{N, F}(::Type{F}, ::String, ::Int64, ::Int64, ::Strand, ::Int8, ::LongSubSeq{DNAAlphabet{N}}, ::NamedTuple) where {N, F<:GeneFinder.GeneFinderMethod}
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:58

Stacktrace:
 [1] top-level scope
   @ REPL[21]:1
```
 
Similar behavior will be encountered when the strand is neither `+` nor `-`. This precautionary measure helps prevent the creation of invalid ORFs, ensuring greater stability and enabling the extension of its interface. For example, after creating a specific `ORFI`, users can seamlessly iterate over a sequence of interest and verify whether the ORFI is contained within the sequence.

```julia
ORFI{4,NaiveFinder}("seq", 1, 33, STRAND_POS, 1, convert(LongSubSeq, dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), NamedTuple())
```

!!! warning
    It is still possible to create an `ORFI` and pass it to a sequence that does not necessarily contain an actual open reading frame. This will be addressed in future versions of the package. But the benefit of having it is that it will retrieve the corresponding subsequence of the sequence in a convinient way (5' to 3') regardless of the strand.


## Finding complete and overlapped ORFIs

The first implemented function is `findorfs` a very non-restrictive ORFI finder function that will catch all ORFIs in a dedicated structure. Note that this will catch random ORFIs not necesarily genes since it has no ORFIs size or overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFIs.

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```

Now lest us find the ORFs

```julia
orfs = findorfs(seq, finder=NaiveFinder)

12-element Vector{ORFI}:
 ORFI{NaiveFinder}(29:40, '+', 2)
 ORFI{NaiveFinder}(137:145, '+', 2)
 ORFI{NaiveFinder}(164:184, '+', 2)
 ORFI{NaiveFinder}(173:184, '+', 2)
 ORFI{NaiveFinder}(236:241, '+', 2)
 ORFI{NaiveFinder}(248:268, '+', 2)
 ORFI{NaiveFinder}(362:373, '+', 2)
 ORFI{NaiveFinder}(470:496, '+', 2)
 ORFI{NaiveFinder}(551:574, '+', 2)
 ORFI{NaiveFinder}(569:574, '+', 2)
 ORFI{NaiveFinder}(581:601, '+', 2)
 ORFI{NaiveFinder}(695:706, '+', 2)
```

Two other methods where implemented into `sequence` to get the ORFIs in DNA or aminoacid sequences, respectively. They use the `findorfs` function to first get the ORFIs and then get the correspondance array of `BioSequence` objects.

```julia
sequence.(orfs)

12-element Vector{NucSeq{4, DNAAlphabet{4}}}
 ATGCAACCCTGA
 ATGCGCTGA
 ATGCGTCGAATGGCACGGTGA
 ATGGCACGGTGA
 ATGTGA
 ATGTGTCCAACGGCAGTCTGA
 ATGCAACCCTGA
 ATGCACTGGCTGGTCCTGTCAATCTGA
 ATGTCACCGCACAAGGCAATGTGA
 ATGTGA
 ATGTGTCCAACGGCAGCCTGA
 ATGCAACCCTGA
```

```julia
transalate.(orfs)

12-element Vector{LongSubSeq{AminoAcidAlphabet}}:
 MQP*
 MR*
 MRRMAR*
 MAR*
 M*
 MCPTAV*
 MQP*
 MHWLVLSI*
 MSPHKAM*
 M*
 MCPTAA*
 MQP*
```