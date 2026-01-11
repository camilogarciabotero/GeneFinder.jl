## Finding Complete and Overlapped ORFs

The main function in the GeneFinder package is `findorfs`, which serves as an interface to various gene-finding algorithms. By default, `findorfs` uses a `NaiveFinder` algorithm, a simple approach that detects all non-outbounded Open Reading Frames (ORFs) in a DNA sequence. You can also specify a different algorithm by setting the `finder` keyword argument. For more details on the NaiveFinder algorithm, see the [NaiveFinder](https://camilogarciabotero.github.io/GeneFinder.jl/dev/api#GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{N}}%20where%20N) documentation.

!!! note
    The `minlen` keyword argument in `NaiveFinder` is set to a minimum length of 6 nucleotides (nt). As a result, it may identify short ORFs that aren't necessarily genes, such as `dna"ATGTGA"` producing the amino acid sequence `aa"M*"`.


## Usage Example

Here's an example of using `findorfs` with the `NaiveFinder` algorithm to identify ORFs in a DNA sequence:

```julia
julia> using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"

orfs = findorfs(seq, finder=NaiveFinder) # use finder=NaiveCollector as an alternative

12-element Vector{ORF{4, NaiveFinder}}:
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

## Extracting Sequences from ORFs

The `ORF` structure displays the location, frame, and strand, along with other fields (see more about the [ORF struct](https://camilogarciabotero.github.io/GeneFinder.jl/dev/orftype#The-ORF-type)). To extract the sequence of an `ORF` instance, you can use the `sequence` method directly on it, or you can also broadcast it over the `orfs` collection using the dot syntax `.`:

```julia
julia> sequence.(orfs)

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

## Translating ORFs to Amino Acid Sequences

Similarly, you can extract the amino acid sequences of the ORFs using the `translate` function:

```julia
julia> translate.(orfs)

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

This returns a vector of translated amino acid sequences, allowing for easy interpretation of each ORF's potential protein product.