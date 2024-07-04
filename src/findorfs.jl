export findorfs

"""
    findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; ::F, kwargs...) where {N, F<:GeneFinderMethod}

This is the main interface method for finding open reading frames (ORFs) in a DNA sequence.

It takes the following required arguments:

- `sequence`: The nucleic acid sequence to search for ORFs.
- `finder`: The algorithm used to find ORFs. It can be either `NaiveFinder`, `NaiveCollector` or yet other implementations.

## Keyword Arguments regardless of the finder method:
- `alternative_start::Bool`: A boolean indicating whether to consider alternative start codons. Default is `false`.
- `minlen::Int`: The minimum length of an ORF. Default is `6`.
- `scheme::Function`: The scoring scheme to use for scoring the sequence from the ORF. Default is `nothing`.

## Returns
A vector of `ORF` objects representing the found ORFs.

## Example

```julia
sequence = randdnaseq(120)

120nt DNA Sequence:
 GCCGGACAGCGAAGGCTAATAAATGCCCGTGCCAGTATC…TCTGAGTTACTGTACACCCGAAAGACGTTGTACGCATTT

findorfs(sequence, finder=NaiveFinder)

1-element Vector{ORF}:
 ORF{NaiveFinder}(77:118, '-', 2, 0.0)
```

"""
function findorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    finder::Type{F}=NaiveFinder,
    kwargs...
) where {N, F<:GeneFinderMethod}
    return finder(sequence; kwargs...)::Vector{ORF{N,F}}
end


### Some possible ideas:

# function findorfs(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}},
#     method::Type{M};
#     kwargs...
# ) where {N, M<:GeneFinderMethod}
#     return method(sequence; kwargs...)::Vector{ORF} #alternative_start, minlen, scoringscheme
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
#         orfs = Vector{ORF}()
#         for strand in ('+', '-')
#             seq = strand == '-' ? reverse_complement(sequence) : sequence

#             @inbounds for location in @views _locationiterator(seq; alternative_start)
#                 if length(location) >= minlen
#                     frame = strand == '+' ? framedict[location.start % 3] : framedict[(seqlen - location.stop + 1) % 3]
#                     start = strand == '+' ? location.start : seqlen - location.stop + 1
#                     stop = start + length(location) - 1
#                     score = -10log10(dnaseqprobability(seq[start:stop], scoringscheme)) # orfs[argmax([orf.score for orf in orfs])]
#                     push!(orfs, ORF(start:stop, strand, frame, score))
#                 end
#             end
#         end
        
#         return sort!(orfs; by = orf -> getproperty(orf, byf), alg=QuickSort, kwargs...)
#     end

# end