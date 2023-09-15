## Roadmap

### Coding genes (CDS - ORFs)

-   ☒ [Finding ORFs](https://camilogarciabotero.github.io/GeneFinder.jl/dev/simplefinder/)
-   ☐ EasyGene
-   ☐ GLIMMER
-   ☐ Prodigal - Pyrodigal
-   ☐ PHANOTATE
-   ☐ k-mer based gene finders (?)
-   ☐ Augustus (?)

### Non-coding genes (RNA)

-   ☐ Infernal
-   ☐ tRNAscan

## Other features

-   ☐ parallelism SIMD ?
-   ☐ memory management (?)
-   ☐ incorporate Ribosime Binding Sites (RBS)
-   ☐ incorporate Programmed Reading Frame Shifting (PRFS)
-   ☐ specialized types
    -   ☒ Gene
    -   ☒ ORF
    -   ☒ Codon
    -   ☒ CDS
    -   ☐ EukaryoticGene (?)
    -   ☐ ProkaryoticGene (?)
    -   ☐ Intron
    -   ☐ Exon
    -   ☐ GFF –\> See other packages
    -   ☐ FASTX –\> See I/O in other packages

## Compatibilities

Must interact with or extend:

- GenomicAnnotations.jl
- BioSequences.jl
- SequenceVariation.jl
- GenomicFeatures.jl
- FASTX.jl
- Kmers.jl
- Graphs.jl
