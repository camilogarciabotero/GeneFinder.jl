# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog and this project adheres to Semantic Versioning.

## [UNRELEASED](https://github.com/camilogarciabotero/GeneFinder.jl/compare/v0.5.0...main)

## [0.5.0]

- The `ORF` type has been updated to handle the new `finder` argument type.
- `findorfs` now has a new `finder` argument that can be used to select the gene finding method.
- The `getorf` method has been removed in favor of having a simple `sequence` method for `ORF` instances.
- The `write_*` API has been updated to handle the new `finder` argument.
- A `NaiveFinder` method has been implemented with a scoring scheme directly.
- `iscoding` has been updated to handle the new `finder` argument.
- Documentation has been updated and improved.


## [0.4.0]

- `findorfs` now has a new `method` argument that can be used to select the gene finding method.
- The `getorf` method has been updated to handle the new `method` argument.
- The `write_*` API has been updated to handle the new `method` argument.
- A `NaiveFinder` method has been implemented with a scoring scheme as well (`NaiveFinderScored`).
- `iscoding` has been updated to handle the new `method` argument.

## [0.3.0]

- `ORF` has a new `score` field that can be used to rank ORFs.
- `findorfs` is now felixible and can be expanded to more gene finding methods/algorithms.
- A new (and quite simple) finding method with a scoring scheme has been implemented.
- Docuemntation has been updated and improved.

## [0.2.0]

- A more stable ORF type that asserts for incorrect ORFs.
- Extended `Base.getindex` method to better handle ORF enabling: `seq[orf]` regardless of the ORF strand.
- The `write_*` API is much more stable and its been rewriten to handle IOStreams and `IOBuffers` as well.
- The new `record_orfs_fna` and `record_orfs_faa` methods intearact with the FASTX interface.
- Codebase has been reduced. Direct String manipulation is not enabled for `findorfs` and `get_*` methods.
- Docstrings have been updated following method updates.
- Better tests and now `Aqua.jl` is used for several tests as well.

## [0.1.0]

- `GeneFinder` is now more reliable and load faster.

## [0.0.23]

- Some cleanups and docs updates

## [0.0.22]

- Decoupling Markov chains to its own repo plus some general cleaning.

## [0.0.21]

- Improve the orf prediction and add frame information see [#20](https://github.com/camilogarciabotero/GeneFinder.jl/pull/20)

## [0.0.20]

- Better docs w/ MathJax

## [0.0.19]

- Improve docs on markov chains and transition models

## [0.0.18]

- Clean code and update methods for transition matrix

## [0.0.17]

- Add several methods towards HMMs in gene finding see PR [#17](https://github.com/camilogarciabotero/GeneFinder.jl/pull/17)

## [0.0.16]

- Fix a bug with importing `FASTX.sequence()` in `write` API
 
## [0.0.15]

- Add precompilation statement to `findorfs()` 
- Clean src and update docs

## [0.0.14]

- Implement a sorting method for the orfs  ([#16](https://github.com/camilogarciabotero/GeneFinder.jl/pull/16))

## [0.0.13]

- Add `write_gff` methods

## [0.0.12]

- some refactors

## [0.0.11]

- Fix a bug with `@view`
- Improve `cdsgenerator` and `getcds`

## [0.0.10]

- Fix a bug in `getcds`

## [0.0.9]

- Improvements in `findorfs` (prev orf_finder). Main ideas: [Dan](https://discourse.julialang.org/u/dan/summary) also see this [Post](https://discourse.julialang.org/t/how-to-improve-a-generator-to-be-more-memory-efficient-when-it-is-collected/92932)
- Update some docs

## [0.0.8]

- Renaming functions in simple finder implementation
- Better docs following some BluStyle guidelines

## [0.0.7]

- Clean up comments.
- New `write_bed` function.
- Document methods.
- Add Zenodo badge.

## [0.0.6]

- Work on extended start codon for generators.
- `simplefind` and its generators now can be fine tuned with alternative codons and orf length.
- Work on the IO system
- More on Docs

## [0.0.5]

- Changed in the internals of `simplefinder`.
- Added a new logo.
- Improved tests and docs.
- Added generators `orfgenerator` and `locationgenerator`.
- Included some helper functions.
- The `findcds` and `findproteins` turned into generators (`cdsgenerator` and `proteingenerator`).
- Add the Codon type and clean others.

## [0.0.4]

- Now registered: `add GeneFinder` to install from the REPL
- Better docs

## [0.0.3]

- Clean `simplefinder`
- Add references
- Clean README
- Clean initial structs
- Make it public

## [0.0.2]

- Add `findcds` and `findproteins` function for the first algorithm
- Add docs for `simplefinder` functions and structs
- Add test using `TestItems` package

## [0.0.1]

- Add some basic `structs` like `ORF`
- Add first ORF finder function `simplefinder`
- Add package overview
