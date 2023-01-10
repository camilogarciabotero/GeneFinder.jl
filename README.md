<p align="center">
  <img src="docs/assets/logo.svg" height="150"><br/>
  <i>A Gene Finder framework for Julia.</i><br/><br/>
  <a href="https://camilogarciabotero.github.io/GeneFinder.jl/dev/">
    <img src="https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white">
  </a>
  <a href="https://github.com/camilogarciabotero/GeneFinder.jl/releases/latest"> 
  <img src="https://img.shields.io/github/release/camilogarciabotero/GeneFinder.jl.svg">
    <a href="https://doi.org/10.5281/zenodo.7519185"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7519185.svg" alt="DOI"></a>
  </a>
  <a href="https://app.travis-ci.com/camilogarciabotero/GeneFinder.jl">
    <img src="https://app.travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main">
   <a href="https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml">
    <img src="https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg">
  <a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green.svg">
  </a>
  <a href="https://www.repostatus.org/#wip">
    <img src="https://www.repostatus.org/badges/latest/wip.svg">
  </a>
</p>

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

## Example

The first implemented function is `simplefinder` a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```
### Finding all ORFs

```julia
orf_finder(seq)

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
```
Two other functions (`simplecds_generator` and `simpleprot_generator`) pass the sequence to `simplefinder` take the ORFs and act as generators of the sequence, so this way the can be `collect`ed in the REPL as an standard output or `write`en into a file more conviniently using the `FASTX` IO system:

```julia
get_cds(seq)

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
```

```julia
get_proteins(seq)

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
### Combining `FASTX` to read a fasta record

```julia
using FASTX

filename = "test/data/KK037166.fna"
rdr = FASTA.Reader(open(filename))
record = first(rdr)
seq = sequence(record)
dnaseq = LongDNA{4}(seq)
```

```julia
get_proteins(dnaseq)

1573-element Vector{LongAA}:
 MMP*
 MP*
 MPATGADL*
 MPCRVASLAR*
 MRSLPSSMARIAPPRTSWDRLPIMPPVRWNRYSARPTRVPASWLCRRRVDSRVATSPAQLSR*
 MARIAPPRTSWDRLPIMPPVRWNRYSARPTRVPASWLCRRRVDSRVATSPAQLSR*
 MPPVRWNRYSARPTRVPASWLCRRRVDSRVATSPAQLSR*
 MPA*
 MTGNRFCPARPGCRRRTASSWVRPSWVRAVVLPVPDAPETIRPRRAETWRWLSWTSSRPVVMTCRITGVCASGSCAW*
 MPASDREFVGSAELGQGGGLAGARRAGDDQAPARGDLAV…LPDHRRVRERELRVVVPSRLVQAEPLGIRQGQQIGVGR*
 MTCRITGVCASGSCAW*
 MTGSTSSTGRAGTCWGRGRLVIATFRLRRRRRRRRRGSR…SCPARARPGSAARAGSARPQPRSRRCRTRPPGGGRGAR*
 MARCVQGTVIQPSGVRRKTWARS*
 MRMP*
 MP*
 MPSRSQARVKVSETNTLPWSMTTWSGTMTGRAAASRIRA…XXXXXXXXXXXYSTVRSATPSSAALARCRCARLVVAGS*
 MTTWSGTMTGRAAASRIRASMASSRWCGSRDADIRRASA…XXXXXXXXXXXYSTVRSATPSSAALARCRCARLVVAGS*
 ⋮
```
### Writting cds, proteins fastas and bed file

```julia
write_cds("cds.fasta", seq)
```

```bash
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

```julia
write_proteins("proteins.fasta", seq)
```

```bash
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
```julia
write_bed("cds.bed", seq)
```
```
cat cds.bed

29	40	+
137	145	+
164	184	+
173	184	+
236	241	+
248	268	+
362	373	+
470	496	+
551	574	+
569	574	+
581	601	+
695	706	+
```

## Algorithms

### Coding genes (CDS - ORFs)

- [x] [Simple finder](https://camilogarciabotero.github.io/GeneFinder.jl/dev/simplefinder/)
- [ ] EasyGene
- [ ] GLIMER3
- [ ] Prodigal - Pyrodigal
- [ ] PHANOTATE
- [ ] k-mer based gene finders (?)
- [ ] Augustus (?)

### Non-coding genes (RNA)

- [ ] Infernal
- [ ] tRNAscan

## Contributing

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).