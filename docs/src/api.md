```@meta
CurrentModule = GeneFinder
DocTestSetup = quote
    using GeneFinder
end
```

## The Main ORF type

The main type of the package is `ORFI` which represents an Open Reading Frame Interval. It is a subtype of the `GenomicInterval` type from the `GenomicFeatures` package.

```@autodocs
Modules = [GeneFinder]
Pages = ["types.jl"]
```

## Finding ORFIs

The function `findorfs` serves as a method interface as it is generic method that can handle different gene finding methods.

```@autodocs
Modules = [GeneFinder]
Pages = ["findorfs.jl"]
```

## Finding ORFs using BioRegex

```@autodocs
Modules = [GeneFinder]
Pages = ["algorithms/naivefinder.jl", "algorithms/naivecollector.jl"]
```

## Writing ORFs to files

```@autodocs
Modules = [GeneFinder]
Pages = ["io.jl"]
```