module GeneFinder

using Automa: RegExp
using BioSequences:
    DNA,
    DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N,
    
    NucleicAcidAlphabet,
    DNAAlphabet,
    AminoAcidAlphabet,
    LongDNA,
    LongAA,
    LongSequence,
    LongSubSeq,
    @biore_str,
    @dna_str,
    GeneticCode,
    reverse_complement

using BioMarkovChains: BioMarkovChain
using FASTX: FASTA, sequence, FASTAReader
using IterTools: takewhile, iterated
using PrecompileTools: @setup_workload, @compile_workload
using StatsBase: countmap
using TestItems: @testitem

include("types.jl")
export ORF

include("algorithms/findorfs.jl")
export locationiterator, findorfs, getorfdna, getorfaa

include("io.jl")
export write_cds, write_proteins, write_bed, write_gff

include("utils.jl")
export fasta_to_dna, nucleotidefreqs, hasprematurestop, dnaseqprobability, iscoding

include("extended.jl")

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    using BioSequences
    seq = randdnaseq(99)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        findorfs(seq)
        getorfdna(seq)
        getorfaa(seq)
    end
end

end
