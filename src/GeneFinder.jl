module GeneFinder

using BioSequences:
    DNA,
    DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N,
    
    Alphabet,
    NucleicAcidAlphabet,
    DNAAlphabet,
    AminoAcidAlphabet,
    LongDNA,
    LongAA,
    LongSequence,
    LongSubSeq,
    NucleicSeqOrView,
    @biore_str,
    @dna_str,
    GeneticCode,
    reverse_complement,
    randdnaseq,
    ncbi_trans_table,
    translate

using BioMarkovChains: BioMarkovChain, dnaseqprobability, ECOLICDS, ECOLINOCDS, log_odds_ratio_score
using FASTX: FASTAReader, sequence
using IterTools: takewhile, iterated
using PrecompileTools: @setup_workload, @compile_workload

include("algorithms/naivefinder.jl")
include("types.jl")
include("findorfs.jl")
include("iscoding.jl")
include("getorfs.jl")
include("io.jl")
include("utils.jl")
include("extended.jl")

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    seq = randdnaseq(99)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        findorfs(NaiveFinder(), seq)
    end
end

end