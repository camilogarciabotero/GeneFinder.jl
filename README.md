# GeneFinder

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/dev/)
[![Build Status](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main)](https://travis-ci.com/camilogarciabotero/GeneFinder.jl)
[![Coverage](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl)


> Warning this is a work in progress.

## Overview

>This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) *gene finder* library for the Julia Language

The main idea is to create versatile module that enables apply different implemented algorithm to DNA sequences. See the BioAlignment implementation of different Seq. Alignment algorithms...(local, global, edit-distance)


## Algorithms

### Coding genes (CDS - ORFs)

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
  - [ ] ORF
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

## Some annotations



## Tests and Benches



## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
