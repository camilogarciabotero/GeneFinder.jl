# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog and this project adheres to Semantic Versioning.

## [UNRELEASED](https://github.com/camilogarciabotero/GeneFinder.jl/compare/v0.0.6...main)


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
