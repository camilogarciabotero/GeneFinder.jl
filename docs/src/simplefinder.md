
## Finding complete and internal (overlapped) ORFs

The first implemented function is `findorfs` a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```
Now lest us find the ORFs

```julia
findorfs(seq)

12-element Vector{ORF}:
 ORF(29:40, '+', 2)
 ORF(137:145, '+', 2)
 ORF(164:184, '+', 2)
 ORF(173:184, '+', 2)
 ORF(236:241, '+', 2)
 ORF(248:268, '+', 2)
 ORF(362:373, '+', 2)
 ORF(470:496, '+', 2)
 ORF(551:574, '+', 2)
 ORF(569:574, '+', 2)
 ORF(581:601, '+', 2)
 ORF(695:706, '+', 2)
```

Two other functions (`get_orfs_dna` and `get_orfs_aa`) pass the sequence to `findorfs` take the ORFs and act as generators of the sequence, so this way the can be `collect`ed in the REPL as an standard output or writteen into a file more conviniently using the `FASTX` IO system:

```julia
get_orfs_dna(seq)

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

```julia
get_orfs_aa(seq)

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

### Writting cds, proteins fastas, bed and gffs whether from a `LongSeq` or from a external fasta file.

```julia
write_cds("cds.fasta", seq)
```

```bash
cat cds.fasta

>location=29:40 strand=+ frame=2
ATGCAACCCTGA
>location=137:145 strand=+ frame=2
ATGCGCTGA
>location=164:184 strand=+ frame=2
ATGCGTCGAATGGCACGGTGA
>location=173:184 strand=+ frame=2
ATGGCACGGTGA
>location=236:241 strand=+ frame=2
ATGTGA
>location=248:268 strand=+ frame=2
ATGTGTCCAACGGCAGTCTGA
>location=362:373 strand=+ frame=2
ATGCAACCCTGA
>location=470:496 strand=+ frame=2
ATGCACTGGCTGGTCCTGTCAATCTGA
>location=551:574 strand=+ frame=2
ATGTCACCGCACAAGGCAATGTGA
>location=569:574 strand=+ frame=2
ATGTGA
>location=581:601 strand=+ frame=2
ATGTGTCCAACGGCAGCCTGA
>location=695:706 strand=+ frame=2
ATGCAACCCTGA
```

### Combining `FASTX` for reading and writing fastas

```julia
using FASTX

write_orfs_aa("test/data/NC_001884.fasta", "proteins.fasta")
```

```bash
head proteins.fasta

>location=41:145 strand=- frame=2
MTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
>location=41:172 strand=- frame=2
MVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
>location=41:454 strand=- frame=2
MSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
>location=41:472 strand=- frame=2
MKTKKQMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
>location=41:505 strand=- frame=2
MLSKYEDDNSNMKTKKQMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
```

## The *log-odds ratio* decision rule

The sequence probability given a transition probability model (eq. 2)
could be used as the source of a sequence classification based on a
decision rule to classify whether a sequence correspond to a model or
another. Now, imagine we got two DNA sequence transition models, a CDS
model and a No-CDS model. The *log-odds ratio* decision rule could be
establish as:

``` math
\begin{align}
S(X) = \log \frac{{P_C(X_1=i_1, \ldots, X_T=i_T)}}{{P_N(X_1=i_1, \ldots, X_T=i_T)}}  \begin{cases} > \eta & \Rightarrow \text{coding} \\ < \eta & \Rightarrow \text{noncoding} \end{cases}
\end{align}
```

Where the ``P_{C}`` is the probability of the sequence given a
CDS model, ``P_{N}`` is the probability of the sequence given a
No-CDS model, the decision rule is finally based on whether the ratio is
greater or lesser than a given threshold *η* of significance level.

In the GeneFinder we have implemented this rule and a couple of basic
transition probability models of CDS and No-CDS of *E. coli* from
Axelson-Fisk (2015) work. To check whether a random sequence could be
coding based on these decision we use the predicate `iscoding` with the
`ECOLICDS` and `ECOLINOCDS` models:

``` julia
randseq = get_orfs_dna(randdnaseq(99))[1] # this will retrieved a random coding ORF

iscoding(randseq, ECOLICDS, ECOLINOCDS)
```

    true