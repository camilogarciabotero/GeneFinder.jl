# Methods from main packages that expand their fuctions to this package structs

@inline Base.isless(a::ORF{F}, b::ORF{F}) where {F} = isless(a.range, b.range)
@inline Base.length(orf::ORF{F}) where {F} = length(sequence(orf))
@inline Base.getindex(seq::SeqOrView{DNAAlphabet{N}}, orf::ORF{F}) where {N,F} = _orfseq(seq, leftposition(orf), rightposition(orf), strand(orf))
@inline Base.range(orf::ORF{F}) where {F} = orf.range

function Base.:(==)(a::ORF{F}, b::ORF{F}) where {F}
    return a.range == b.range && a.strand == b.strand && a.frame == b.frame && a.seqid == b.seqid
end

## Methods from BioSequences that expand their fuctions to this package structs
import BioSequences: translate
@inline translate(orf::ORF{F}; kwargs...) where {F} = translate(sequence(orf); kwargs...)

## Methods from BioMarkovChains that expand their fuctions to this package structs
import BioMarkovChains: log_odds_ratio_score

@inline log_odds_ratio_score(orf::ORF{F}; kwargs...) where {F} = log_odds_ratio_score(sequence(orf); kwargs...)
