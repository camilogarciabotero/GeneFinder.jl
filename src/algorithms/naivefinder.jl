export naivefinder

"""
    locationiterator(sequence::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}

This is an iterator function that uses regular expressions to search the entire ORF (instead of start and stop codons) in a `LongSequence{DNAAlphabet{4}}` sequence.
    It uses an anonymous function that will find the first regularly expressed ORF. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.

!!! note
    As a note of the implementation we want to expand on how the ORFs are found:

    The expression `(?:[N]{3})*?` serves as the boundary between the start and stop codons. Within this expression, the character class `[N]{3}` captures exactly three occurrences of any character (representing nucleotides using IUPAC codes). This portion functions as the regular codon matches. Since it is enclosed within a non-capturing group `(?:)` and followed by `*?`, it allows for the matching of intermediate codons, but with a preference for the smallest number of repetitions. 
    
    In summary, the regular expression `ATG(?:[N]{3})*?T(AG|AA|GA)` identifies patterns that start with "ATG," followed by any number of three-character codons (represented by "N" in the IUPAC code), and ends with a stop codon "TAG," "TAA," or "TGA." This pattern is commonly used to identify potential protein-coding regions within genetic sequences.

    See more about the discussion [here](https://discourse.julialang.org/t/how-to-improve-a-generator-to-be-more-memory-efficient-when-it-is-collected/92932/8?u=camilogarciabotero)

"""
function _locationiterator(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false
) where {N}
    regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    # regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG([N]{3})*T(AG|AA|GA)?"dna # an attempt to make it non PCRE non-determinsitic
    finder(x) = findfirst(regorf, sequence, first(x) + 1) # + 3
    itr = takewhile(!isnothing, iterated(finder, findfirst(regorf, sequence)))
    return itr
end

"""
naivefinder(sequence::NucleicAlphabet{DNAAlphabet{N}}; alternative_start::Bool=false, min_len::Int64=6)::Vector{ORF} where {N}

A simple implementation that finds ORFs in a DNA sequence.

The `naivefinder` function takes a LongSequence{DNAAlphabet{4}} sequence and returns a Vector{ORF} containing the ORFs found in the sequence. 
    It searches entire regularly expressed CDS, adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
        Extending the starting codons with the `alternative_start = true` will search for ATG, GTG, and TTG.
    Some studies have shown that in *E. coli* (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.
!!! note
    This function has not ORFs scoring scheme. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.
    
# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function naivefinder(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    orfs = Vector{ORF}(undef, 0)
    reversedseq = reverse_complement(sequence)
    seqlen = length(sequence)
    framedict = Dict(0 => 3, 1 => 1, 2 => 2)

    for strand in ('+', '-')
        seq = strand == '-' ? reversedseq : sequence

        @inbounds for location in @views _locationiterator(seq; alternative_start)
            if length(location) >= min_len
                frame = strand == '+' ? framedict[location.start % 3] : framedict[(seqlen - location.stop + 1) % 3]
                push!(orfs, ORF(strand == '+' ? location : (seqlen - location.stop + 1):(seqlen - location.start + 1), strand, frame, 0.0))
            end
        end
    end
    return sort(orfs)
end


## Another alternative:

# function _locationiterator(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}};
#     alternative_start::Bool = false
# ) where {N}
#     regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
#     finder(x) = findfirst(regorf, sequence, first(x) + 1) # + 3
#     itr = takewhile(!isnothing, iterated(finder, findfirst(regorf, sequence)))
#     return itr
# end

# function naivefinder(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}};
#     alternative_start::Bool = false,
#     min_len::Int64 = 6
# ) where {N}
#     orfs = Vector{ORF}(undef, 0)
#     reversedseq = reverse_complement(sequence)
#     seqlen = length(sequence)
#     framedict = Dict(0 => 3, 1 => 1, 2 => 2)

#     for strand in ('+', '-')
#         seq = strand == '-' ? reversedseq : sequence

#         @inbounds for location in @views _locationiterator(seq; alternative_start)
#             if length(location) >= min_len
#                 frame = strand == '+' ? framedict[location.start % 3] : framedict[(seqlen - location.stop + 1) % 3]
#                 push!(orfs, ORF(strand == '+' ? location : (seqlen - location.stop + 1):(seqlen - location.start + 1), strand, frame, 0.0))
#             end
#         end
#     end
#     return sort(orfs)
# end

## Ideas to genefinder (this is the way window function are implemented in the DSP.jl package)

# function make_orf_finder(finderfunc::Function, sequence::NucleicSeqOrView{DNAAlphabet{N}}) where {N}
#     ...
# end

# function locationiterator(sequence::NucleicSeqOrView{DNAAlphabet{N}}) where {N}
#     make_orf_finder(sequence::NucleicSeqOrView{DNAAlphabet{N}}) do x
#         ...
#     end
# end

# function findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; finder::Function=locationiterator) where {N}
#     ...
# end