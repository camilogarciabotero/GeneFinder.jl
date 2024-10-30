## Finding complete and overlapped ORFIs

The main package function is `findorfs`. Under the hood, the `findorfs` function is an interface for different gene finding algorithms that can be plugged using the `finder` keyword argument. By default it uses the `NaiveFinder` algorithm, which is a simple algorithm that finds all (non-outbounded) ORFIs in a DNA sequence (see the [NaiveFinder](https://camilogarciabotero.github.io/GeneFinder.jl/dev/api/#GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{N}}%20where%20N) documentation for more details).

!!! note
    The `minlen` kwarg in the `NaiveFinder` mehtod has been set to 6nt, so it will catch random ORFIs not necesarily genes thus it might consider `dna"ATGTGA"` -> `aa"M*"` as a plausible ORFI.

Here is an example of how to use the `findorfs` function with the `NaiveFinder` algorithm:

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"

orfs = findorfs(seq, finder=NaiveFinder) # use finder=NaiveCollector as an alternative

12-element Vector{ORFI{4, NaiveFinder}}:
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

The `ORFI` structure displays the location, frame, and strand, but currently does not include the sequence *per se*. To extract the sequence of an `ORFI` instance, you can use the `sequence` method directly on it, or you can also broadcast it over the `orfs` collection using the dot syntax `.`:

```julia
sequence.(orfs)

12-element Vector{LongSubSeq{DNAAlphabet{4}}}:
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

Similarly, you can extract the amino acid sequences of the ORFIs using the `translate` function.

```julia
translate.(orfs)

12-element Vector{LongAA}:
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