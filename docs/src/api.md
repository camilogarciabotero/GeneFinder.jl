```@meta
CurrentModule = GeneFinder
DocTestSetup = quote
    using GeneFinder
end
```

## The Main ORF type

The main type of the package is `ORF` which represents an Open Reading Frame.

```@autodocs
Modules = [GeneFinder]
Pages = ["types.jl"]
```

## Finding ORFs

The function `findorfs` is the main function of the package. It is generic method that can handle different gene finding methods. 

```@autodocs
Modules = [GeneFinder]
Pages = ["findorfs.jl"]
```

## Finding ORFs using BioRegex and scoring

```@autodocs
Modules = [GeneFinder]
Pages = ["algorithms/naivefinder.jl"]
```

## Writing ORFs to files

```@autodocs
Modules = [GeneFinder]
Pages = ["io.jl"]
```