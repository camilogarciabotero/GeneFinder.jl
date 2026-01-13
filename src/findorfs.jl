export findorfs

"""
    findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; finder::Type{F}, kwargs...) where {N, F<:GeneFinderMethod}

Main interface for finding Open Reading Frames (ORFs) in a DNA sequence.

Returns an `ORFCollection` containing the found ORFs bundled with the source sequence,
providing a clean API for sequence extraction.

# Arguments
- `sequence::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.

# Keywords
- `finder::Type{F}=NaiveFinder`: The algorithm to use (`NaiveFinder`, `NaiveFinderLazy`, etc.).
- `alternative_start::Bool=false`: Whether to consider alternative start codons (GTG, TTG).
- `minlen::Int=6`: Minimum ORF length in nucleotides.

# Returns
- `ORFCollection{F}`: Collection of ORFs bundled with the source sequence.

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = findorfs(seq)

# Access ORFs sequences
orfseqs = sequence.(Ref(collection), collection.orfs)  # Get sequences of all ORFs

# Use a different finder
collection = findorfs(seq, finder=NaiveFinderLazy)
```

See also: [`NaiveFinder`](@ref), [`NaiveFinderLazy`](@ref), [`ORFCollection`](@ref)
"""
function findorfs(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    finder::Type{F}=NaiveFinder,
    kwargs...
) where {N,F<:GeneFinderMethod}
    return finder(seq; kwargs...)::ORFCollection{F}
end



### Some possible ideas:

# function findorfs(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}},
#     method::Type{M};
#     kwargs...
# ) where {N, M<:GeneFinderMethod}
#     return method(sequence; kwargs...)::Vector{ORFI} #alternative_start, minlen, scoringscheme
# end

# This general implementation would allow for the addition of new finder methods without changing the findorfs method.
# The general method have some advantages:
# - It is more flexible and extensible.
# - It is more robust to changes in the method implementations.
# - It is more readable and maintainable.
# - We can set a default method to be used when no method is provided.

# The general method also have some disadvantages:
# The `kwargs` argument is not type-stable, which can lead to performance issues. (?)
# The `kwargs` are not discoverable from the method when its called.
# Will need that finder methods be created inside the struct:

# struct NaiveFinderScored <: GeneFinderMethod

#     function NaiveFinderScored(
#         sequence::NucleicSeqOrView{DNAAlphabet{N}};
#         alternative_start::Bool = false,
#         minlen::Int64 = 6,
#         scoringscheme::BioMarkovChain = ECOLICDS,
#         byf::Symbol = :location,
#         kwargs...
#     ) where {N}
#         seqlen = length(sequence)
#         framedict = Dict(0 => 3, 1 => 1, 2 => 2)
#         orfs = Vector{ORFI}()
#         for strand in ('+', '-')
#             seq = strand == '-' ? reverse_complement(sequence) : sequence

#             @inbounds for location in @views _locationiterator(seq; alternative_start)
#                 if length(location) >= minlen
#                     frame = strand == '+' ? framedict[location.start % 3] : framedict[(seqlen - location.stop + 1) % 3]
#                     start = strand == '+' ? location.start : seqlen - location.stop + 1
#                     stop = start + length(location) - 1
#                     score = -10log10(dnaseqprobability(seq[start:stop], scoringscheme)) # orfs[argmax([orf.score for orf in orfs])]
#                     push!(orfs, ORFI(start:stop, strand, frame, score))
#                 end
#             end
#         end
        
#         return sort!(orfs; by = orf -> getproperty(orf, byf), alg=QuickSort, kwargs...)
#     end

# end