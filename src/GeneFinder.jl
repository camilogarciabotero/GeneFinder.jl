module GeneFinder

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
    NucleicSeqOrView,
    @biore_str,
    @dna_str,
    GeneticCode,
    reverse_complement,
    randdnaseq,
    ncbi_trans_table,
    translate

using FASTX: FASTAReader, FASTARecord, description, sequence
using IterTools: takewhile, iterated
using PrecompileTools: @setup_workload, @compile_workload

include("types.jl")
export ORF

include("algorithms/findorfs.jl")
export locationiterator, findorfs, get_orfs_dna, get_orfs_aa, record_orfs_fna, record_orfs_faa

include("io.jl")
export write_orfs_fna, write_orfs_faa, write_orfs_bed, write_orfs_gff

include("utils.jl")
export fasta_to_dna, hasprematurestop

include("extended.jl")

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    seq = randdnaseq(99)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        findorfs(seq)
        get_orfs_dna(seq)
        # get_orfs_aa(seq)
    end
end

end