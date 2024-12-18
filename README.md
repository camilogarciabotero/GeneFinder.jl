<p align="center">
  <img src="docs/src/assets/logo.svg" height="150"><br/>
  <i>A Gene Finder Framework for the Julia Programming Language.</i>
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

The `GeneFinder` package aims to be a versatile module that enables the application of different gene finding algorithms to the `BioSequence` type, by providing a common interface and a flexible data structure to store the predicted ORFI or genes. The package is designed to be easily extensible, allowing users to implement their own algorithms and integrate them into the framework.

> [!WARNING] 
  This package is currently under development and is not yet ready for production use. The API is subject to change.

## Installation

You can install `GeneFinder` from the julia REPL. Press `]` to enter pkg mode, and enter the following command:

```julia
add GeneFinder
```

## Finding complete and overlapped ORFIs

The main package function is `findorfs`. Under the hood, the `findorfs` function is an interface for different gene finding algorithms that can be plugged using the `finder` keyword argument. By default it uses the `NaiveFinder` algorithm, which is a simple algorithm that finds all (non-outbounded) ORFIs in a DNA sequence (see the [NaiveFinder](https://camilogarciabotero.github.io/GeneFinder.jl/dev/api/#GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{N}}%20where%20N) documentation for more details).

> [!NOTE] 
  The `minlen` kwarg in the `NaiveFinder` mehtod has been set to 6nt, so it will catch random ORFIs not necesarily genes thus it might consider `dna"ATGTGA"` -> `aa"M*"` as a plausible ORFI.

Here is an example of how to use the `findorfs` function with the `NaiveFinder` algorithm:

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"

orfs = findorfs(seq, finder=NaiveFinder) # use finder=NaiveCollector as an alternative

12-element Vector{ORFI{4, NaiveFinder}}:
 ORFI{NaiveFinder}(29:40, '+', 2)
 ORFI{NaiveFinder}(137:145, '+', 2)
 ORFI{NaiveFinder}(164:184, '+', 2)
 ORFI{NaiveFinder}(173:184, '+', 2)
 ORFI{NaiveFinder}(236:241, '+', 2)
 ORFI{NaiveFinder}(248:268, '+', 2)
 ORFI{NaiveFinder}(362:373, '+', 2)
 ORFI{NaiveFinder}(470:496, '+', 2)
 ORFI{NaiveFinder}(551:574, '+', 2)
 ORFI{NaiveFinder}(569:574, '+', 2)
 ORFI{NaiveFinder}(581:601, '+', 2)
 ORFI{NaiveFinder}(695:706, '+', 2)
```

The `ORFI` structure displays the location, frame, and strand, but currently does not include the sequence *per se*. To extract the sequence of an `ORFI` instance, you can use the `sequence` method directly on it, or you can also broadcast it over the `orfs` collection using the dot syntax `.`:

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

Similarly, you can extract the amino acid sequences of the ORFIs using the `translate` function.

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

## Let's score the ORFIs

ORFIs sequences can be scored using different schemes that evaluate them under a biological context. There are two ways to make this possible: by adding a scoring method to the finder algorithm or by using a scoring method after predicting the ORFIs. The first approach is likely more efficient, but the second approach is more flexible. We will showcase the second approach in this example.

A commonly used scoring scheme for ORFIs is the *log-odds ratio* score. This score is based on the likelihood of a sequence belonging to a specific stochastic model, such as coding or non-coding. The [BioMarkovChains](https://github.com/camilogarciabotero/BioMarkovChains.jl) package provides a `log_odds_ratio_score` method (currently imported), also known as `lors`, which can be used to score ORFIs using the log-odds ratio approach.

```julia
orfs = findorfs(seq, finder=NaiveFinder)
```

The `lors` method has been overloaded to take an ORFI object and can be used later to calculate the score of the ORFIs.

```julia
lors.(orfs)

12-element Vector{Float64}:
 0.469404606944017
 1.0174520899042823
 1.5914902556997463
 0.9772187907841964
 0.6106494455192994
 0.7089167973379216
 0.469404606944017
 1.5523291911446804
 0.5282685400427601
 0.6106494455192994
 0.7405746713921604
 0.469404606944017
```

We can extend basically any method that scores a `BioSequence` to score an `ORFI` object. To see more about scoring ORFIs, check out the [Scoring ORFIs](https://camilogarciabotero.github.io/GeneFinder.jl/dev/features/) section in the documentation.

## Writting ORFIs   into bioinformatic formats

`GeneFinder` also now facilitates the generation of `FASTA`, `BED`, and `GFF` files directly from the found ORFIs. This feature is particularly useful for downstream analysis and visualization of the ORFIs. To accomplish this, the package provides the following functions: `write_orfs_fna`, `write_orfs_faa`, `write_orfs_bed`, and `write_orfs_gff`.

Functionality:

The package provides four distinct functions for writing files in different formats:

| Function          | Description                                            |
|-------------------|--------------------------------------------------------|
| `write_orfs_fna`    | Writes nucleotide sequences in FASTA format.     |
| `write_orfs_faa`    | Writes amino acid sequences in FASTA format.  |
| `write_orfs_bed`    | Outputs information in BED format.                           |
| `write_orfs_gff`    | Generates files in GFF format.                              |

All these function support processing `BioSequences` instances. To demonstrate the use of the `write_*` methods with a `BioSequence`, consider the following example:

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```

Once a `BioSequence` object has been created, the `write_orfs_fna` function proves useful for generating a `FASTA` file containing the nucleotide sequences of the ORFIs. Notably, the `write_orfs*` methods support either an `IOStream` or an `IOBuffer` as an output argument, allowing flexibility in directing the output either to a file or a buffer. In the following example, we demonstrate writing the output directly to a file.

```julia
outfile = "LFLS01000089.fna"

open(outfile, "w") do io
    write_orfs_fna(seq, io, finder=NaiveFinder) # use finder=NaiveCollector as an alternative
end
```

```bash
cat LFLS01000089.fna

>seq id=01 start=29 stop=40 strand=+ frame=2 features=[]
ATGCAACCCTGA
>seq id=02 start=137 stop=145 strand=+ frame=2 features=[]
ATGCGCTGA
>seq id=03 start=164 stop=184 strand=+ frame=2 features=[]
ATGCGTCGAATGGCACGGTGA
>seq id=04 start=173 stop=184 strand=+ frame=2 features=[]
ATGGCACGGTGA
>seq id=05 start=236 stop=241 strand=+ frame=2 features=[]
ATGTGA
>seq id=06 start=248 stop=268 strand=+ frame=2 features=[]
ATGTGTCCAACGGCAGTCTGA
>seq id=07 start=362 stop=373 strand=+ frame=2 features=[]
ATGCAACCCTGA
>seq id=08 start=470 stop=496 strand=+ frame=2 features=[]
ATGCACTGGCTGGTCCTGTCAATCTGA
>seq id=09 start=551 stop=574 strand=+ frame=2 features=[]
ATGTCACCGCACAAGGCAATGTGA
>seq id=10 start=569 stop=574 strand=+ frame=2 features=[]
ATGTGA
>seq id=11 start=581 stop=601 strand=+ frame=2 features=[]
ATGTGTCCAACGGCAGCCTGA
>seq id=12 start=695 stop=706 strand=+ frame=2 features=[]
ATGCAACCCTGA
```

This could also be done to writting a `FASTA` file with the nucleotide sequences of the ORFIs using the `write_orfs_fna` function. Similarly for the `BED` and `GFF` files using the `write_orfs_bed` and `write_orfs_gff` functions respectively.
