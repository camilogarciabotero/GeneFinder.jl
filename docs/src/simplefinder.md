
# A simple algorithm

The first implemented function is `orf_finder` a very non-restrictive
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
orf_finder(seq)
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

Two other functions (`get_cds` and `get_proteins`) pass the sequence to
`orf_finder` take the ORFs and act as generators of the sequence, so
this way the can be `collect`ed in the REPL as an standard output or
written into a file more conviniently using the `FASTX` IO system:

## Generting cds and proteins with its ORF

``` julia
get_cds(seq)
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
get_proteins(seq)
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

## Combining `FASTX` for reading and writing a fasta record

``` julia
using FASTX

write_proteins("../../test/data/NC_001884.fasta", "proteins.fasta")
```

``` bash
head proteins.fasta
>locus=75:113 strand=+
MKLNLRIGVISN*
>locus=144:215 strand=+
MLTITSFKTILNSSFFFSELDSM*
>locus=210:215 strand=+
M*
>locus=237:374 strand=+
MLFLTVLLSISDCVSCNPLSSFFAFWSSLNSSSNAAFLFKKSSSL*
>locus=337:402 strand=+
MQLFSSKKVHHCKCHFHIYRR*
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
