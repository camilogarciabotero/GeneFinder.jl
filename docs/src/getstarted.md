## Finding Complete and Overlapped ORFs

The main function in the GeneFinder package is `findorfs`, which serves as an interface to various gene-finding algorithms. All methods return an `ORFCollection`—a structure that bundles the found ORFs with their source sequence, providing a clean API for sequence extraction.

!!! note
    The `minlen` keyword argument is set to a minimum length of 6 nucleotides (nt). As a result, it may identify short ORFs that aren't necessarily genes, such as `dna"ATGTGA"` producing the amino acid sequence `aa"M*"`.


## Usage Example

Here's an example of using `findorfs` with the `NaiveFinder` algorithm:

```julia
using BioSequences, GeneFinder

seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"

collection = findorfs(seq, finder=NaiveFinder)
```

Output:
```
ORFCollection{NaiveFinder} with 12 ORFs in 726bp sequence:
 ORF{NaiveFinder}(29:40, 'PSTRAND', 2)
 ORF{NaiveFinder}(137:145, 'PSTRAND', 2)
 ORF{NaiveFinder}(164:184, 'PSTRAND', 2)
 ⋮
```

## Working with ORFCollections

The `ORFCollection` provides a convenient interface for working with ORFs:

```julia
# Iteration
for orf in collection
    println(leftposition(orf), " - ", rightposition(orf))
end

# Indexing
first_orf = collection[1]
subset = collection[1:5]

# Length
n_orfs = length(collection)

# Access source sequence
src = source(collection)
```

## Extracting Sequences from ORFs

Extract sequences using the `sequence` function with an index or ORF:

```julia
# By index
dna_seq = sequence(collection, 1)

# By ORF object
orf = collection[1]
dna_seq = sequence(collection, orf)
```

## Translating ORF Sequences

To translate ORF sequences to amino acids, use `BioSequences.translate` on extracted sequences:

```julia
using BioSequences

# Translate a single ORF
protein = translate(sequence(collection, 1))

# Translate all ORFs
proteins = [translate(sequence(collection, i)) for i in eachindex(collection)]
```