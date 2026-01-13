
## Finding Complete and Overlapped ORFs {#Finding-Complete-and-Overlapped-ORFs}

The main function in the GeneFinder package is `findorfs`, which serves as an interface to various gene-finding algorithms. All methods return an `ORFCollection`—a structure that bundles the found ORFs with their source sequence, providing a clean API for sequence extraction.

::: tip Note

The `minlen` keyword argument is set to a minimum length of 6 nucleotides (nt). As a result, it may identify short ORFs that aren&#39;t necessarily genes, such as `dna"ATGTGA"` producing the amino acid sequence `aa"M*"`.

:::

## Usage Example {#Usage-Example}

Here&#39;s an example of using `findorfs` with the `NaiveFinder` algorithm:

```julia
julia> using BioSequences, GeneFinder

julia> seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC";

julia> collection = findorfs(seq, finder=NaiveFinder)
ORFCollection{NaiveFinder} with 12 ORFs in 726bp sequence:
 ORF{NaiveFinder}(29:40, '+', 2)
 ORF{NaiveFinder}(137:145, '+', 2)
 ORF{NaiveFinder}(164:184, '+', 2)
 ⋮
```


## Working with ORFCollections {#Working-with-ORFCollections}

The `ORFCollection` provides a convenient interface for working with ORFs:

```julia
julia> collection[1]  # First ORF
ORF{NaiveFinder}(29:40, '+', 2)

julia> collection[1:5]  # Range of ORFs (returns a Vector)
5-element Vector{OpenReadingFrame{NaiveFinder}}:
 ORF{NaiveFinder}(29:40, '+', 2)
 ⋮

julia> length(collection)  # Number of ORFs
12

julia> source(collection)  # Access source sequence
726nt DNA Sequence:
AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC
```


## Extracting Sequences from ORFs {#Extracting-Sequences-from-ORFs}

Extract sequences using the `sequence` function with an index or ORF:

```julia
julia> sequence(collection, 1)  # By index
12nt DNA Sequence:
ATGCAACCCTGA

julia> orf = collection[1]; sequence(collection, orf)  # By ORF
12nt DNA Sequence:
ATGCAACCCTGA

julia> sequence.(Ref(collection), collection.orfs)  # All sequences using broadcasting
12-element Vector{LongSubSeq{DNAAlphabet{4}}}:
 ATGCAACCCTGA
 ATGCGCTGA
 ⋮
```


## Translating ORF Sequences {#Translating-ORF-Sequences}

To translate ORF sequences to amino acids, use `BioSequences.translate` on extracted sequences:

```julia
julia> using BioSequences

julia> translate(sequence(collection, 1))
4aa Amino Acid Sequence:
MQP*

julia> translate.(sequence.(Ref(collection), collection.orfs))  # Translate all ORFs
12-element Vector{LongAA}:
 MQP*
 MR*
 MRRMAR*
 ⋮
```

