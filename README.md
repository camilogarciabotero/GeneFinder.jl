<p align="center">
  <img src="docs/src/assets/logo.svg" height="150"><br/>
  <i>A Gene Finder framework for Julia.</i>
</p>

<div align="center">

[![Documentation](https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white)](https://camilogarciabotero.github.io/GeneFinder.jl/dev/)
[![Release](https://img.shields.io/github/release/camilogarciabotero/GeneFinder.jl.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/releases/latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7519184.svg)](https://doi.org/10.5281/zenodo.7519184)<br>
[![GitHub Actions](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE)
[![Repo Status](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/GeneFinder&label=downloads)](https://pkgs.genieframework.com?packages=GeneFinder)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

</div>

***

<!-- [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) -->

## Overview

>This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) *gene finder* library framework for the Julia Language.

The main goal of `GeneFinder` is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).

## Installation

You can install GeneFinder from the julia REPL. Press `]` to enter pkg mode, and enter the following:

```julia
add GeneFinder
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

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