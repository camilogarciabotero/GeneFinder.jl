# Methods from main packages that expand their fuctions to this package structs
# import Base: isless, getindex, length, ==

@inline Base.isless(a::ORFI{N,F}, b::ORFI{N,F}) where {N,F} = isless(a.first:a.last, b.first:b.last)
@inline Base.length(orf::ORFI{N,F}) where {N,F} = length(sequence(orf))
@inline Base.getindex(seq::SeqOrView{DNAAlphabet{N}}, orf::ORFI{N,F}) where {N,F} =  _orfseq(seq, orf.first, orf.last, orf.strand)

function Base.:(==)(a::ORFI{N,F}, b::ORFI{N,F}) where {N,F}
    return a.first == b.first && a.last == b.last && a.strand == b.strand && a.frame == b.frame && a.groupname == b.groupname && a.seq == b.seq # all(name -> getproperty(a, name) == getproperty(b, name), propertynames(a))
end

## Methods from BioSequences that expand their fuctions to this package structs
import BioSequences: translate
@inline translate(orf::ORFI{N,F}; kwargs...) where {N,F} =  translate(sequence(orf); kwargs...)

## Methods from BioMarkovChains that expand their fuctions to this package structs
import BioMarkovChains: log_odds_ratio_score
export log_odds_ratio_score, lors

@inline log_odds_ratio_score(orf::ORFI{N,F}; kwargs...) where {N,F} = log_odds_ratio_score(sequence(orf); kwargs...)
@inline log_odds_ratio_decision_rule(orf::ORFI{N,F}; kwargs...) where {N,F} = log_odds_ratio_decision_rule(sequence(orf); kwargs)

const  lors = log_odds_ratio_score
