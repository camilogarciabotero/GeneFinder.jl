module GeneFinderTests

using Test
using BioSequences
using FASTX
using GeneFinder
using Aqua

"""
    fasta_to_dna(input::String)

Converts a FASTA formatted file (even if it is a multi-fasta) to an array of `LongSequence{DNAAlphabet{4}}` objects.
"""
function fasta_to_dna(input::AbstractString)::Vector{LongSequence{DNAAlphabet{4}}}
    FASTAReader(open(input)) do reader
        return [LongSequence{DNAAlphabet{4}}(sequence(record)) for record in reader]
    end
end

include("findorfstest.jl")
include("iotest.jl")
include("getindextest.jl")
include("aquatest.jl")

end
