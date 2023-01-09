
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

## Finding all ORFs

``` julia
simplefind(seq)
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

## Generting cds and proteins with its ORF

``` julia
cds = simplecds_generator(seq)
```

    Base.Generator{Base.Iterators.Filter{GeneFinder.var"#4#6"{Int64}, Vector{ORF}}, GeneFinder.var"#3#5"{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}}(GeneFinder.var"#3#5"{LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}(AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC, GCTCAGGTTACCGCCGGCAGTCAGGGTTGCATTGCCCGC…CAGGGTTGCATTGCCCGCGGTACTGATATTGCCCTGGTT), Base.Iterators.Filter{GeneFinder.var"#4#6"{Int64}, Vector{ORF}}(GeneFinder.var"#4#6"{Int64}(6), ORF[ORF(29:40, '+'), ORF(137:145, '+'), ORF(164:184, '+'), ORF(173:184, '+'), ORF(236:241, '+'), ORF(248:268, '+'), ORF(362:373, '+'), ORF(470:496, '+'), ORF(551:574, '+'), ORF(569:574, '+'), ORF(581:601, '+'), ORF(695:706, '+')]))

``` julia
[i.sequence for i in cds]
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
protein = simpleprot_generator(seq)
```

    Base.Generator{Base.Iterators.Filter{GeneFinder.var"#14#16"{Int64}, Vector{ORF}}, GeneFinder.var"#13#15"{Bool, BioSequences.GeneticCode, LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}}(GeneFinder.var"#13#15"{Bool, BioSequences.GeneticCode, LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}(false, The Standard Code, AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC, GCTCAGGTTACCGCCGGCAGTCAGGGTTGCATTGCCCGC…CAGGGTTGCATTGCCCGCGGTACTGATATTGCCCTGGTT), Base.Iterators.Filter{GeneFinder.var"#14#16"{Int64}, Vector{ORF}}(GeneFinder.var"#14#16"{Int64}(6), ORF[ORF(29:40, '+'), ORF(137:145, '+'), ORF(164:184, '+'), ORF(173:184, '+'), ORF(236:241, '+'), ORF(248:268, '+'), ORF(362:373, '+'), ORF(470:496, '+'), ORF(551:574, '+'), ORF(569:574, '+'), ORF(581:601, '+'), ORF(695:706, '+')]))

``` julia
[i.sequence for i in protein]
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

## Combining `FASTX` to read a fasta record

``` julia
using FASTX

filename = "../../test/data/KK037166.fna"
rdr = FASTA.Reader(open(filename))
record = first(rdr)
seq = sequence(record)
dnaseq = LongDNA{4}(seq)
[i.sequence for i in simpleprot_generator(dnaseq)]
```

## Writting cds and proteins fastas

``` julia
write_cds("cds.fasta", seq)
```

``` bash
cat cds.fasta
>locus=29:40 strand=+
ATGCAACCCTGA
>locus=137:145 strand=+
ATGCGCTGA
>locus=164:184 strand=+
ATGCGTCGAATGGCACGGTGA
>locus=173:184 strand=+
ATGGCACGGTGA
>locus=236:241 strand=+
ATGTGA
>locus=248:268 strand=+
ATGTGTCCAACGGCAGTCTGA
>locus=362:373 strand=+
ATGCAACCCTGA
>locus=470:496 strand=+
ATGCACTGGCTGGTCCTGTCAATCTGA
>locus=551:574 strand=+
ATGTCACCGCACAAGGCAATGTGA
>locus=569:574 strand=+
ATGTGA
>locus=581:601 strand=+
ATGTGTCCAACGGCAGCCTGA
>locus=695:706 strand=+
ATGCAACCCTGA
```

``` julia
write_proteins("proteins.fasta", seq)
```

``` bash
cat proteins.fasta
>locus=29:40 strand=+
MQP*
>locus=137:145 strand=+
MR*
>locus=164:184 strand=+
MRRMAR*
>locus=173:184 strand=+
MAR*
>locus=236:241 strand=+
M*
>locus=248:268 strand=+
MCPTAV*
>locus=362:373 strand=+
MQP*
>locus=470:496 strand=+
MHWLVLSI*
>locus=551:574 strand=+
MSPHKAM*
>locus=569:574 strand=+
M*
>locus=581:601 strand=+
MCPTAA*
>locus=695:706 strand=+
MQP*
```
