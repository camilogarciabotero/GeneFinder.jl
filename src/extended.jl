# Methods from main packages that expand their fuctions to this package structs
import BioSequences: translate, standard_genetic_code
import Base: show, length, iterate, sort

function translate(
    ntseq::LongSubSeq{DNAAlphabet{4}};
    code::GeneticCode = standard_genetic_code,
    allow_ambiguous_codons = true,
    alternative_start = false,
)
    ntseq = copy(ntseq)
    translate(ntseq; code, allow_ambiguous_codons, alternative_start)
end

Base.sort(v::Vector{<:ORF}) = sort(v, by = _orf_sort_key)

function _orf_sort_key(orf::ORF)
    return (orf.location, orf.strand)
end