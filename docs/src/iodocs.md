### Writting cds, proteins fastas, bed and gffs whether from a `LongSeq` or from a external fasta file.

This package provides a set of functions to write `FASTA`, `BED` and `GFF` files with the ORF information out of a `BioSequence` (more specifically `NucleicSeqOrView{A} where A`). The `FASTA` files can be written with the `write_orfs_fna` or `write_orfs_faa` function, the `BED` files with the `write_orfs_bed` function and the `GFF` files with the `write_orfs_gff` function. All of these functions can be used with a `NucleicSeqOrView` or with a external `FASTA` file (using the `String` input as a path to a `FASTA` file containing a genome for instance). The `write_orfs_faa` function can also be used to write a `FASTA` file with the protein sequences of the ORFs.


```julia

```

```bash
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

write_orfs_faa("test/data/NC_001884.fasta", "proteins.fasta")
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