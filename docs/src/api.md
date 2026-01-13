```@meta
CurrentModule = GeneFinder
DocTestSetup = quote
    using GeneFinder
end
```

## Core Types

The main types of the package are `ORF` (Open Reading Frame) and `ORFCollection`. 
An `ORF` stores coordinates and metadata, while `ORFCollection` bundles ORFs with 
their source sequence for clean sequence extraction.

```@autodocs
Modules = [GeneFinder]
Pages = ["types.jl"]
```

## Finding ORFs

The `findorfs` function serves as a unified interface for different gene finding methods.
All methods return an `ORFCollection`.

```@autodocs
Modules = [GeneFinder]
Pages = ["findorfs.jl"]
```

## ORF Finding Algorithms

### NaiveFinder

Uses regular expression matching to find ORFs.

```@autodocs
Modules = [GeneFinder]
Pages = ["algorithms/naivefinder.jl"]
```

## Writing ORFs to Files

Export ORFs in various formats (FASTA, BED, GFF3).

```@autodocs
Modules = [GeneFinder]
Pages = ["io.jl"]
```