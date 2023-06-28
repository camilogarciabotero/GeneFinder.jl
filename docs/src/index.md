
![]("../assets/logo.svg")

[![](https://img.shields.io/github/release/camilogarciabotero/GeneFinder.jl.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/releases/latest)
[![](https://app.travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main)](https://app.travis-ci.com/camilogarciabotero/GeneFinder.jl)
[![](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE)
[![](https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white)](https://camilogarciabotero.github.io/GeneFinder.jl/dev/)
[![](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![](https://zenodo.org/badge/DOI/10.5281/zenodo.7519184.svg)](https://doi.org/10.5281/zenodo.7519184)


***
## Overview

> This is a species-agnostic, algorithm extensible, sequence-anonymous
> (genome, metagenomes) *gene finder* library for the Julia Language.

The main goal is to create a versatile module that enables apply
different implemented algorithm to DNA sequences. See, for instance,
BioAlignment implementations of different sequence alignment algorithms
(local, global, edit-distance).

## Installation

You can install GeneFinder from the julia REPL. Press `]` to enter pkg
mode, and enter the following:

    add GeneFinder

If you are interested in the cutting edge of the development, please
check out the master branch to try new features before release.

## Algorithms

### Coding genes (CDS - ORFs)

-   ☒ [Simple finder](https://camilogarciabotero.github.io/GeneFinder.jl/dev/simplefinder/)
-   ☐ EasyGene
-   ☐ GLIMMER
-   ☐ Prodigal - Pyrodigal
-   ☐ PHANOTATE
-   ☐ k-mer based gene finders (?)
-   ☐ Augustus (?)

### Non-coding genes (RNA)

-   ☐ Infernal
-   ☐ tRNAscan

## Other features

-   ☐ parallelism SIMD ?
-   ☐ memory management (?)
-   ☐ incorporate Ribosime Binding Sites (RBS)
-   ☐ incorporate Programmed Reading Frame Shifting (PRFS)
-   ☐ specialized types
    -   ☒ Gene
    -   ☒ ORF
    -   ☒ Codon
    -   ☒ CDS
    -   ☐ EukaryoticGene (?)
    -   ☐ ProkaryoticGene (?)
    -   ☐ Intron
    -   ☐ Exon
    -   ☐ GFF –\> See other packages
    -   ☐ FASTX –\> See I/O in other packages

## Compatibilities

Must interact with or extend:

- GenomicAnnotations.jl
- BioSequences.jl
- SequenceVariation.jl
- GenomicFeatures.jl
- FASTX.jl
- Kmers.jl
- Graphs.jl

## Contributing

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
