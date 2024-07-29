# Methods from main packages that expand their fuctions to this package structs
import Base: isless, iterate, sort, getindex, length, ==
import BioSequences: translate

Base.isless(a::ORF{N,F}, b::ORF{N,F}) where {N,F} = isless(a.first:a.last, b.first:b.last)

function Base.getindex(seq::NucleicSeqOrView{A}, orf::ORF{N,F}) where {A,N,F}
    # @assert @view(seq[begin:end]) == @view(source(orf)[begin:end]) "The source sequence of the ORF and the given sequence are different"
    # @view(seq[begin:end]) == @view(source(orf)[begin:end]) || throw(ArgumentError("The source sequence of the ORF and the given sequence are different"))
    # seq[begin:end] == source(orf) || throw(ArgumentError("The source sequence of the ORF and the given sequence are different"))
    
    if orf.strand == STRAND_NEG
        s = reverse_complement(seq[orf.first:orf.last])
    else
        s = seq[orf.first:orf.last]
    end

    if !occursin(biore"T(AG|AA|GA)"dna, s)
        error("There is no stop codon at the end of the sequence/orf")
    end

    return s
end

function Base.length(i::ORF{N,F}) where {N,F}
    return length(sequence(i))
end

function translate(orf::ORF{N,F}; kwargs...) where {N,F}
    return translate(sequence(orf); kwargs...)
end

# Define custom equality function for ORF structs
function ==(a::ORF{N,F}, b::ORF{N,F}) where {N,F}
    return a.first == b.first && a.last == b.last && a.strand == b.strand && a.frame == b.frame && a.groupname == b.groupname && a.seq == b.seq
end

## Methods from BioMarkovChains that expand their fuctions to this package structs
import BioMarkovChains: log_odds_ratio_score

export log_odds_ratio_score, lors

function log_odds_ratio_score(
    orf::ORF{N,F};
    kwargs...
) where {N,F}
    return log_odds_ratio_score(sequence(orf); kwargs...)
end

const  lors = log_odds_ratio_score

function lordr(orf::ORF{N,F}; kwargs...) where {N,F}
    return lordr(sequence(orf); kwargs)
end