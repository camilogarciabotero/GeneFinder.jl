[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://camilogarciabotero.github.io/GeneFinder.jl/dev/) 
[![Build Status](https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main)](https://travis-ci.com/camilogarciabotero/GeneFinder.jl)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE)

***

[![Coverage](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/camilogarciabotero/GeneFinder.jl)
[![CI](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml)

## Overview

>This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) *gene finder* library for the Julia Language.

The main idea is to create versatile module that enables apply different implemented algorithm to DNA sequences. See the BioAlignment implementation of different sequence alignment algorithms (local, global, edit-distance).

## Installation

You can install GeneFinder from the julia
REPL. Press `]` to enter pkg mode, and enter the following:

```
add GeneFinder
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

## Algorithms

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
