# <img src="./sticker.svg" width="30%" align="right" /> GeneFinder

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/stable/)

[![Latest Release](https://img.shields.io/github/release/camilogarciabotero/GeneFinder.jl.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/releases/latest)

[![Build Status](https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main)](https://travis-ci.com/camilogarciabotero/GeneFinder.jl)

[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE)

<!-- [![Build Status](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml?query=branch%3Amain) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/dev/) -->

***

[![Coverage](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl)

[![CI](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml)

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

<!-- [![Unit tests](https://github.com/camilogarciabotero/GeneFinder.jl/workflows/Unit%20tests/badge.svg?branch=main)](https://github.com/camilogarciabotero/GeneFinder.jl/actions?query=workflow%3A%22Unit+tests%22+branch%3Amain) -->

***

> Warning this is a work in progress.

## Overview

>This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) *gene finder* library for the Julia Language.

The main idea is to create versatile module that enables apply different implemented algorithm to DNA sequences. See the BioAlignment implementation of different sequence alignment algorithms (local, global, edit-distance).

## Installation

You can install BioSequences from the julia
REPL. Press `]` to enter pkg mode, and enter the following:

```
add GeneFinder
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

## Algorithms

### Coding genes (CDS - ORFs)

- [x] Simple finder (a simple ORF finder)
- [ ] Prodigal - Pyrodigal
- [ ] k-mer based (?)
- [ ] Augustus (?)
- [ ] PHANOTATE

### Non-coding genes (RNA)

- [ ] Infernal
- [ ] tRNAscan

## Other features

- [ ] parallelism SIMD ?
- [ ] memory management (?)
- [ ] specialized types
  - [ ] Gene
  - [x] ORF
  - [ ] CDS
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

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
