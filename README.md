# GeneFinder

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/dev/)
[![Build Status](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main)](https://travis-ci.com/camilogarciabotero/GeneFinder.jl)
[![Coverage](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl)

WIP
>This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) *gene finder* library for the Julia Lang

The main idea is to create versatile module that enables apply different implemented algorithm to DNA sequences.

See the BioAlignment implementation of different Seq. Alignment algorithms...(local, global, edit-distance)


## Algorithms

### Coding genes (CDS - ORFs)

- [ ] Prodigal - Pyrodigal
- [ ] PHANOTATE
- [ ] k-mer based (?)
- [ ] Augustus (?)

### Non-coding genes (RNA)

- [ ] Infernal
- [ ] tRNAscan

## Other features

- [ ] parallelism SIMD ?
- [ ] memory management (?)
- [ ] ...

## Compatibilities  

Must interact with:

- GenomicAnnotations.jl
- BioSequences.jl
- SequenceVariation.jl
- GenomicFeatures.jl
- FASTX.jl
- Kmers.jl


## Some annotations

There might be several abstract types:

- ProkaryoticGene
- EukaryoticGene
- GeneModel
- Algorithm
- ORF ?
- Codons ? 

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
