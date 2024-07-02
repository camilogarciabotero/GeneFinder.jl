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
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FGeneFinder&query=total_requests&suffix=%2Fmonth&label=Downloads)](http://juliapkgstats.com/pkg/GeneFinder)
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

## Finding complete and overlapped ORFs

The first implemented function is `findorfs` a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```
Now lest us find the ORFs

```julia
orfs = findorfs(seq, finder=NaiveFinder) # use finder=NaiveCollector as an alternative

12-element Vector{ORF{4, NaiveFinder}}:
 ORF{NaiveFinder}(29:40, '+', 2)
 ORF{NaiveFinder}(137:145, '+', 2)
 ORF{NaiveFinder}(164:184, '+', 2)
 ORF{NaiveFinder}(173:184, '+', 2)
 ORF{NaiveFinder}(236:241, '+', 2)
 ORF{NaiveFinder}(248:268, '+', 2)
 ORF{NaiveFinder}(362:373, '+', 2)
 ORF{NaiveFinder}(470:496, '+', 2)
 ORF{NaiveFinder}(551:574, '+', 2)
 ORF{NaiveFinder}(569:574, '+', 2)
 ORF{NaiveFinder}(581:601, '+', 2)
 ORF{NaiveFinder}(695:706, '+', 2)
```

The `ORF` structure contains also the sequence of the ORF, the frame, and several features. If you want to get the sequence of the ORFs you can broadcast the `sequence` function to safely extract the sequences of all ORFs.

```julia
sequence.(orfs)

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

Similarly, you can extract the amino acid sequences of the ORFs using the `translate` function.

```julia
translate.(orfs)

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

## Writting ORF information into bioinformatic formats

This package facilitates now the creation of `FASTA`, `BED`, and `GFF` files, specifically extracting Open Reading Frame (ORF) information from `BioSequence` instances, particularly those of type `NucleicSeqOrView{A} where A`, and then writing the information into the desired format.

Functionality:

The package provides four distinct functions for writing files in different formats:

| Function          | Description                                            |
|-------------------|--------------------------------------------------------|
| `write_orfs_fna`    | Writes nucleotide sequences in FASTA format.     |
| `write_orfs_faa`    | Writes amino acid sequences in FASTA format.  |
| `write_orfs_bed`    | Outputs information in BED format.                           |
| `write_orfs_gff`    | Generates files in GFF format.                              |


All these functions support processing both `BioSequence` instances and external `FASTA` files. In the case of a `BioSequence` instace into external files, simply provide the path to the `FASTA` file using a `String` to the path. To demonstrate the use of the `write_*` methods with a `BioSequence`, consider the following example:

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```

Once a `BioSequence` object has been instantiated, the `write_orfs_fna` function proves useful for generating a `FASTA` file containing the nucleotide sequences of the ORFs. Notably, the `write_orfs*` methods support either an `IOStream` or an `IOBuffer` as an output argument, allowing flexibility in directing the output either to a file or a buffer. In the following example, we demonstrate writing the output directly to a file.

```julia
outfile = "LFLS01000089.fna"

open(outfile, "w") do io
    write_orfs_fna(seq, io, finder=NaiveFinder) # use finder=NaiveCollector as an alternative
end
```

```bash
cat LFLS01000089.fna

>seq id=01 start=29 stop=40 strand=+ frame=2 score=0.0
ATGCAACCCTGA
>seq id=02 start=137 stop=145 strand=+ frame=2 score=0.0
ATGCGCTGA
>seq id=03 start=164 stop=184 strand=+ frame=2 score=0.0
ATGCGTCGAATGGCACGGTGA
>seq id=04 start=173 stop=184 strand=+ frame=2 score=0.0
ATGGCACGGTGA
>seq id=05 start=236 stop=241 strand=+ frame=2 score=0.0
ATGTGA
>seq id=06 start=248 stop=268 strand=+ frame=2 score=0.0
ATGTGTCCAACGGCAGTCTGA
>seq id=07 start=362 stop=373 strand=+ frame=2 score=0.0
ATGCAACCCTGA
>seq id=08 start=470 stop=496 strand=+ frame=2 score=0.0
ATGCACTGGCTGGTCCTGTCAATCTGA
>seq id=09 start=551 stop=574 strand=+ frame=2 score=0.0
ATGTCACCGCACAAGGCAATGTGA
>seq id=10 start=569 stop=574 strand=+ frame=2 score=0.0
ATGTGA
>seq id=11 start=581 stop=601 strand=+ frame=2 score=0.0
ATGTGTCCAACGGCAGCCTGA
>seq id=12 start=695 stop=706 strand=+ frame=2 score=0.0
ATGCAACCCTGA
```

This could also be done to writting a `FASTA` file with the nucleotide sequences of the ORFs using the `write_orfs_fna` function. Similarly for the `BED` and `GFF` files using the `write_orfs_bed` and `write_orfs_gff` functions respectively.