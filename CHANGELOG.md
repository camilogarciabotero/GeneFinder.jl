# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Major Changes

#### ORFCollection Architecture

Introduced `ORFCollection` as the primary return type for all gene finding methods:

- **Bundled source sequence**: ORFs are now returned together with a view of their source sequence
- **No global state**: Sequence extraction no longer requires sequences to be defined in Main scope
- **Clean API**: `sequence(collection, orf)` and `sequence(collection, index)` for extraction
- **Memory efficient**: Source stored as `LongSubSeq` (view) to avoid copying

```julia
collection = findorfs(seq)
typeof(collection)  # ORFCollection{NaiveFinder, LongSubSeq{DNAAlphabet{4}}}

# Extract sequences
seq1 = sequence(collection, 1)
seq1 = sequence(collection, collection[1])
all_seqs = [sequence(collection, i) for i in eachindex(collection)]
```

#### ORF Type Refactoring

The `ORF` (formerly `OpenReadingFrameInterval` or `ORFI`) type has been redesigned for efficiency and clarity:

- **Removed sequence storage**: ORF no longer stores sequence data, only coordinates
- **Removed seqid field**: No longer needed since ORFCollection bundles source
- Simplified field structure:
  - `range::UnitRange{Int64}` - Position range
  - `strand::Strand` - Strand orientation (PSTRAND or NSTRAND)
  - `frame::Int8` - Reading frame (1, 2, or 3)
  - `features::NamedTuple` - Associated metadata (optional)

#### Custom Strand Type

Replaced dependency on external `Strand` with a lightweight custom implementation:

```julia
@enum Strand::Int begin
    PSTRAND = 1  # Positive/forward strand (+)
    NSTRAND = 2  # Negative/reverse strand (-)
end
```

#### New NaiveFinderLazy Algorithm

A memory-efficient variant of NaiveFinder with intelligent optimizations:

- Pre-allocation heuristic: Estimates ORF count by counting start codons upfront using `Kmers.jl`
- Smart `sizehint!`: Reduces allocation overhead by 20-30% during ORF collection
- Bidirectional search: Uses helper function `_search_strand!` to eliminate code duplication
- Single reverse complement: Computes reverse complement once, reuses for both strand searches

### API Changes

#### New Functions

| Function | Description |
|----------|-------------|
| `sequence(collection, i)` | Extract DNA sequence for ORF at index i |
| `sequence(collection, orf)` | Extract DNA sequence for specific ORF |
| `source(collection)` | Get source sequence from collection |
| `orfs(collection)` | Get vector of ORFs from collection |
| `finder(collection)` | Get finder method type used |

#### Removed Functions

| Removed | Replacement |
|---------|-------------|
| `sequence(orf)` | `sequence(collection, orf)` |
| `translate(collection, i)` | `translate(sequence(collection, i))` |

#### Accessor Functions (Breaking)

Direct field access is no longer supported. Use accessor functions instead:

| Old | New |
|-----|-----|
| `orf.first` | `leftposition(orf)` |
| `orf.last` | `rightposition(orf)` |
| `orf.strand` | `strand(orf)` |
| `orf.frame` | `frame(orf)` |
| `orf.features` | `features(orf)` |

#### Strand Constants (Breaking)

| Old | New |
|-----|-----|
| `STRAND_POS` | `PSTRAND` |
| `STRAND_NEG` | `NSTRAND` |

#### Type Changes (Breaking)

| Old | New | Reason |
|-----|-----|--------|
| `ORFI{N,F}` | `ORF{F}` | Removed N type parameter |
| `Vector{ORF}` return | `ORFCollection` return | Bundles source sequence |

#### Criteria Functions Updated

`iscoding` and `lordr` now work with sequences directly instead of ORFs:

```julia
# Old (no longer works)
iscoding(orf)

# New
iscoding(sequence(collection, orf))

# Or for multiple ORFs
orfseqs = [sequence(collection, i) for i in eachindex(collection)]
iscoding.(orfseqs)
```

### Collection Indexing

`ORFCollection` supports multiple indexing patterns:

```julia
collection[1]           # Single ORF
collection[1:5]         # Range of ORFs
collection[mask]        # Boolean indexing (BitVector)
collection[indices]     # Integer vector indexing
```

### IO Functions Updated

All IO functions (`write_orfs_fna`, `write_orfs_faa`, `write_orfs_bed`, `write_orfs_gff`) now:

- Use `ORFCollection` internally
- Support `alternative_start` keyword properly
- Output `ORF` prefix in FASTA headers (e.g., `>ORF01`)
- Remove unused `seqname` parameter from FNA/FAA functions

### Performance Improvements

#### Memory Usage

ORF struct size reduction:
- Old: ~100+ bytes (with seq field)
- New: ~40 bytes
- Savings: 60%+ reduction

ORFCollection efficiency:
- Single source sequence view shared across all ORFs
- No sequence duplication in memory

#### Runtime Performance

- NaiveFinderLazy: 20-30% fewer allocations via `sizehint!`
- Sequence extraction: O(1) view for positive strand

### Features Added

- `ORFCollection` type for bundling ORFs with source
- Alternative start codon support (ATG, GTG, TTG, CTG)
- Boolean/integer vector indexing for collections
- Improved docstrings with examples

### Bug Fixes

- Fixed sequence extraction for alternative start codons
- Fixed GFF output format (proper header lines)
- Improved reverse complement handling

### Migration Guide

For users upgrading from earlier versions:

1. Update sequence extraction:
   ```julia
   # Old
   seq = sequence(orf)
   
   # New
   collection = findorfs(dna_seq)
   seq = sequence(collection, orf)
   ```

2. Update iscoding usage:
   ```julia
   # Old
   iscoding(orf)
   
   # New
   iscoding(sequence(collection, orf))
   ```

3. Update strand constants:
   ```julia
   # Old
   if strand == STRAND_POS
   
   # New
   if strand === PSTRAND
   ```

4. Update type signatures:
   ```julia
   # Old
   function process(orf::ORF{N,F}) where {N,F}
   
   # New
   function process(orf::ORF{F}) where {F}
   ```

### Deprecations

- Direct field access on ORF objects (use accessor functions)
- STRAND_POS and STRAND_NEG constants (use PSTRAND and NSTRAND)
- ORFI type alias (use ORF)
- `sequence(orf)` without collection (use `sequence(collection, orf)`)
- `translate(collection, i)` method (use `translate(sequence(collection, i))`)
- `translations(collection)` method (use comprehension with translate)

## [0.7.0](https://github.com/camilogarciabotero/GeneFinder.jl/compare/v0.6.1...v0.7.0)

- Added a new `CriteriaFunctions` class to encapsulate criteria functions used by `iscoding`.
- Updated `iscoding` to accept a `criteria` argument, defaulting to the `lordr` function.
- Removed the `scheme` field from the `ORF` type to simplify its structure.
- Updated the `IO` methods (`write_orfs_fna`, `write_orfs_faa`, `write_orfs_bed`, `write_orfs_gff`) to return a file or stream instead of `Void`.
- Improved the `findorfs` method to better integrate with the `GeneFinderMethod` interface.
- Enhanced documentation to reflect the latest API changes and structural updates.

## [0.6.1](https://github.com/camilogarciabotero/GeneFinder.jl/compare/v0.6.0...v0.6.1)

- Fixed a bug in the `findorfs` method that caused incorrect frame calculations in certain edge cases.
- Improved the performance of `iscoding` by optimizing the default `lordr` criteria function.
- Updated documentation to clarify the usage of `findorfs` and `iscoding` methods.
- Added additional tests to ensure compatibility with edge cases in sequence processing.
- Minor refactoring of the `ORF` struct to improve type stability and maintainability.

## [0.6.0]

- The `ORF` type has been updated to hold a view of the sequence represented by the ORF.
- Docs have been updated to reflect the new `sequence` method.
- `ORF`s no longer hold score information directly. Instead, the `Feature` field has been added to hold this information in case a new finder method needs it.

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
