# # this will be the main functions taking all the algorithms
using BioSequences

include("types.jl")
include("helpers.jl")

"""
    locationgenerator(sequence::LongDNA)

Generate the locations of ORFs in the given DNA `sequence`.

This function searches the sequence for start codons, and generates ranges of indices corresponding to the locations of ORFs in the `sequence`. 
    The ORFs are generated by iterating over the start codon indices and searching for the first stop codon that follows each start codon. 
        ORFs that contain premature stop codons are filtered out using the `hasprematurestop` function. 
            The `sequence` argument must be a `LongDNA` object, which is a type of DNA sequence with a longer maximum length than the `DNA` type.

Returns:
    A generator expression that yields ranges of indices corresponding to the locations of ORFs in the `sequence`.
"""
function locationgenerator(sequence::LongDNA)
    startcodon = ExactSearchQuery(Codon("ATG"), iscompatible)
    seqbound = length(sequence) - 3
    start_codon_indices = findall(startcodon, sequence)
    @inbounds begin
        return(i.start:j+2 for i in start_codon_indices for j in i.start:3:seqbound if sequence[j:j+2] ∈ stopcodons && !hasprematurestop(sequence[i.start:j+2]))
    end
end

"""
    orfgenerator(sequence::LongDNA)

Generate ORFs from the given DNA `sequence`.

This function generates ORFs from the forward and reverse complement strands of the `sequence` using the `locationgenerator` function. 
    It generates an ORF object for each range of indices returned by `locationgenerator`, and includes a `'+'` or `'-'` strand label 
        to indicate the strand from which the ORF was generated. The `sequence` argument must be a `LongDNA` object, which is a type 
        of DNA sequence with a longer maximum length than the `DNA` type.

Returns:
    A generator expression that yields `ORF` objects corresponding to the ORFs in the `sequence`.
"""
function orfgenerator(sequence::LongDNA)
    reversedseq = reverse_complement(sequence)
    @inbounds begin
        orfs = (ORF(l, s) for s in ['+', '-'] for l in locationgenerator(s == '+' ? sequence : reversedseq))
    end
    return orfs
end


# """
# FindGene struct
# """
# mutable struct FindGene{S1,S2}
#     a::GeneFinderAlgorithm{S1}
#     b::LongDNA
#     c::GeneticCode
# end

# function findgenes(::SimpleFinder, sequence::LongDNA) # type::GeneticCode
#     orfs = simplefinder(sequence)
#     seqs = Vector{CDS}()
    
# end