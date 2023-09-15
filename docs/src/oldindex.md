---
engine: knitr
cache: true
---

<p align="center">
  <img src="../assets/logo.svg" height="150"><br/>
  <i>A Gene Finder framework for Julia.</i><br/><br/>
  <a href="https://www.repostatus.org/#wip">
    <img src="https://www.repostatus.org/badges/latest/wip.svg">
  </a>
  <a href="https://codecov.io/gh/camilogarciabotero/GeneFinder.jl">
    <img src="https://img.shields.io/codecov/c/github/camilogarciabotero/GeneFinder.jl?logo=codecov&logoColor=white">
  </a>
  <a href="https://camilogarciabotero.github.io/GeneFinder.jl/dev/">
    <img src="https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white">
  </a>
  <a href="https://travis-ci.com/camilogarciabotero/GeneFinder.jl">
    <img src="https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main">
  <a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green.svg">
  </a>
</p>

***
## Overview

>This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) *gene finder* library for the Julia Language.

The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).

## Installation

You can install GeneFinder from the julia
REPL. Press `]` to enter pkg mode, and enter the following:

```
add GeneFinder
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

## Algorithms

### Coding genes (CDS - ORFs)

- [x] [Simple finder](https://camilogarciabotero.github.io/GeneFinder.jl/dev/simplefinder/)
- [ ] EasyGene
- [ ] GLIMMER
- [ ] Prodigal - Pyrodigal
- [ ] PHANOTATE
- [ ] k-mer based gene finders (?)
- [ ] Augustus (?)

### Non-coding genes (RNA)

- [ ] Infernal
- [ ] tRNAscan

## Other features

- [ ] parallelism SIMD ?
- [ ] memory management (?)
- [ ] specialized types
  - [x] Gene
  - [x] ORF
  - [x] CDS
  - [ ] EukaryoticGene (?)
  - [ ] ProkaryoticGene (?)
  - [ ] Codon
  - [ ] Intron
  - [ ] Exon
  - [ ] GFF --> See other packages
  - [ ] FASTX --> See I/O in other packages

## Compatibilities  

Must interact with or extend:

- GenomicAnnotations.jl
- BioSequences.jl
- SequenceVariation.jl
- GenomicFeatures.jl
- FASTX.jl
- Kmers.jl

## Contributing

## Citing



*Logo: gene analysis by Vector Points from the Noun Project*