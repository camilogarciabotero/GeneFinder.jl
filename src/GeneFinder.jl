module GeneFinder

using BioSequences:
    DNA,
    DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N,
    
    Alphabet,
    BioRegex,
    NucleicAcidAlphabet,
    DNAAlphabet,
    AminoAcidAlphabet,
    LongDNA,
    LongAA,
    LongSequence,
    LongSubSeq,
    SeqOrView,
    NucleicSeqOrView,
    @biore_str,
    @dna_str,
    GeneticCode,
    reverse_complement,
    randdnaseq,
    ncbi_trans_table,
    translate,
    findnext,

    ExactSearchQuery,
    # ApproximateSearchQuery,
    iscompatible

using BioMarkovChains: BioMarkovChain, ECOLICDS, ECOLINOCDS, log_odds_ratio_score
using IterTools: takewhile, iterated

using BioSymbols: encoded_data
using Kmers: @mer_str, FwDNAMers, FwRvIterator, Kmer

# Main types (must come first - defines Strand, ORF, etc.)
include("types.jl")
include("utils.jl")

# Finder Algorithms (depends on types.jl)
include("algorithms/naivefinder.jl")
include("algorithms/naivecollector.jl")

# Main functions
include("findorfs.jl")
include("io.jl")

# Extended functions
include("extended.jl")

# RBS Scoring
include("rbs.jl")

# Coding Criteria
include("criteria.jl")

# Precompiled workloads
include("workload.jl")

end