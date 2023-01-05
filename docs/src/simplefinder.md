
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

Two other functions (`cdsgenerator` and `proteingenerator`) pass the
sequence to `simplefinder` take the ORFs and act as generators of the
sequence, so this way the can be `collect`ed in the REPL as an standard
output or `write`en into a file more conviniently using the `FASTX` IO
system:

``` julia
cds = cdsgenerator(seq)
```

Line 1  
We actually need to collect the `Generator` to see the CDSs in the
stdout. Other way to to is simply by array comprenhension
`[for cds in cdsgenerator(seq)]`

<!-- -->

    Base.Generator{Vector{ORF}, GeneFinder.var"#1#2"{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}}(GeneFinder.var"#1#2"{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}(AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC, GCTCAGGTTACCGCCGGCAGTCAGGGTTGCATTGCCCGC…CAGGGTTGCATTGCCCGCGGTACTGATATTGCCCTGGTT), ORF[ORF(29:40, '+'), ORF(137:145, '+'), ORF(164:184, '+'), ORF(173:184, '+'), ORF(236:241, '+'), ORF(248:268, '+'), ORF(362:373, '+'), ORF(470:496, '+'), ORF(551:574, '+'), ORF(569:574, '+'), ORF(581:601, '+'), ORF(695:706, '+')])

``` julia
collect(cds)
```

    12-element Vector{LongSequence{DNAAlphabet{4}}}:
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

``` julia
protein = proteingenerator(seq)
```

    Base.Generator{Vector{ORF}, GeneFinder.var"#3#4"{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}}(GeneFinder.var"#3#4"{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}(AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC, GCTCAGGTTACCGCCGGCAGTCAGGGTTGCATTGCCCGC…CAGGGTTGCATTGCCCGCGGTACTGATATTGCCCTGGTT), ORF[ORF(29:40, '+'), ORF(137:145, '+'), ORF(164:184, '+'), ORF(173:184, '+'), ORF(236:241, '+'), ORF(248:268, '+'), ORF(362:373, '+'), ORF(470:496, '+'), ORF(551:574, '+'), ORF(569:574, '+'), ORF(581:601, '+'), ORF(695:706, '+')])

``` julia
collect(protein)
```

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
