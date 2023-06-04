module GeneFinder

using BioSequences
using BioSequences: GeneticCode
using FASTX: FASTA
using IterTools
using TestItems: @testitem
using PrecompileTools

include("types.jl")
export ORF, Codon, CDS, Protein

include("algorithms/findorfs.jl")
export locationiterator, findorfs

include("findgenes.jl")
export cdsgenerator, proteingenerator, getcds, getproteins

include("io.jl")
export write_cds, write_proteins, write_bed, write_gff

include("helpers.jl")
export fasta_to_dna, sort


@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    using BioSequences
    seq = randdnaseq(1000)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        findorfs(seq)
    end
end

end