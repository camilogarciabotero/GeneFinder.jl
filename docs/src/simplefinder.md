
# A simple algorithm

The first implemented function is `simplefinder` a very non-restrictive
ORF finder function that will catch all ORFs in a dedicated structure.
Note that this will catch random ORFs not necesarily genes since it has
no ORFs size or overlapping condition contraints. Thus it might consider
`aa"M*"` a posible encoding protein from the resulting ORFs.

``` julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```

    726nt DNA Sequence:
    AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC

``` julia
simplefinder(seq)
```

    12-element Vector{ORF}:
     ORF(29:40, '+')
     ORF(137:145, '+')
     ORF(164:184, '+')
     ORF(173:184, '+')
     ORF(236:241, '+')
     ORF(248:268, '+')
     ORF(362:373, '+')
     ORF(470:496, '+')
     ORF(551:574, '+')
     ORF(569:574, '+')
     ORF(581:601, '+')
     ORF(695:706, '+')

Two other functions (`findcds` and `findproteins`) pass the sequence to
`simplefinder` take the ORFs to index search the CDS and traslate into
Protein:

``` julia
findcds(seq)
```

    12-element Vector{CDS}:
     CDS(ORF(29:40, '+'), ATGCAACCCTGA)
     CDS(ORF(137:145, '+'), ATGCGCTGA)
     CDS(ORF(164:184, '+'), ATGCGTCGAATGGCACGGTGA)
     CDS(ORF(173:184, '+'), ATGGCACGGTGA)
     CDS(ORF(236:241, '+'), ATGTGA)
     CDS(ORF(248:268, '+'), ATGTGTCCAACGGCAGTCTGA)
     CDS(ORF(362:373, '+'), ATGCAACCCTGA)
     CDS(ORF(470:496, '+'), ATGCACTGGCTGGTCCTGTCAATCTGA)
     CDS(ORF(551:574, '+'), ATGTCACCGCACAAGGCAATGTGA)
     CDS(ORF(569:574, '+'), ATGTGA)
     CDS(ORF(581:601, '+'), ATGTGTCCAACGGCAGCCTGA)
     CDS(ORF(695:706, '+'), ATGCAACCCTGA)

``` julia
findproteins(seq)
```

    12-element Vector{Protein}:
     Protein(ORF(29:40, '+'), MQP*)
     Protein(ORF(137:145, '+'), MR*)
     Protein(ORF(164:184, '+'), MRRMAR*)
     Protein(ORF(173:184, '+'), MAR*)
     Protein(ORF(236:241, '+'), M*)
     Protein(ORF(248:268, '+'), MCPTAV*)
     Protein(ORF(362:373, '+'), MQP*)
     Protein(ORF(470:496, '+'), MHWLVLSI*)
     Protein(ORF(551:574, '+'), MSPHKAM*)
     Protein(ORF(569:574, '+'), M*)
     Protein(ORF(581:601, '+'), MCPTAA*)
     Protein(ORF(695:706, '+'), MQP*)
