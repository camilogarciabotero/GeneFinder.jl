
## The ORF type

For convenience, the `ORF` type is more stringent in preventing the creation of incompatible instances. As a result, attempting to create an instance with incompatible parameters will result in an error. For instance, the following code snippet will trigger an error:

```julia
ORF(1:10, '+', 4)

ERROR: AssertionError: Invalid frame value. Frame must be 1, 2, or 3.
Stacktrace:
 [1] ORF(location::UnitRange{Int64}, strand::Char, frame::Int64)
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:47
 [2] top-level scope
   @ REPL[25]:1
```
 
Similar behavior will be encountered when the strand is neither `+` nor `-`. This precautionary measure helps prevent the creation of invalid ORFs, ensuring greater stability and enabling the extension of its interface. For example, after creating a specific `ORF`, users can seamlessly iterate over a sequence of interest and verify whether the ORF is contained within the sequence.

```julia
orf = ORF(137:145, '+', 2)
seq[orf]

9nt DNA Sequence:
ATGCGCTGA
```

!!! warning
    It is still possible to create an `ORF` and pass it to a sequence that does not necessarily contain an actual open reading frame. This will be addressed in future versions of the package. But the benefit of having it is that it will retrieve the corresponding subsequence of the sequence in a convinient way (5' to 3') regardless of the strand.


## Finding complete and overlapped ORFs

The first implemented function is `findorfs` a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```
Now lest us find the ORFs

```julia
orfs = findorfs(seq, finder=NaiveFinder)

12-element Vector{ORF}:
 ORF{NaiveFinder}(29:40, '+', 2)
 ORF{NaiveFinder}(137:145, '+', 2)
 ORF{NaiveFinder}(164:184, '+', 2)
 ORF{NaiveFinder}(173:184, '+', 2)
 ORF{NaiveFinder}(236:241, '+', 2)
 ORF{NaiveFinder}(248:268, '+', 2)
 ORF{NaiveFinder}(362:373, '+', 2)
 ORF{NaiveFinder}(470:496, '+', 2)
 ORF{NaiveFinder}(551:574, '+', 2)
 ORF{NaiveFinder}(569:574, '+', 2)
 ORF{NaiveFinder}(581:601, '+', 2)
 ORF{NaiveFinder}(695:706, '+', 2)
```

Two other methods where implemented into `sequence` to get the ORFs in DNA or aminoacid sequences, respectively. They use the `findorfs` function to first get the ORFs and then get the correspondance array of `BioSequence` objects.

```julia
sequece.(orfs)

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