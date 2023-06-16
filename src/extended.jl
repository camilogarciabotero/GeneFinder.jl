# Methods from main packages that expand their fuctions to this package structs
import BioSequences: translate, standard_genetic_code
function translate(ntseq::LongSubSeq{DNAAlphabet{4}}; 
    code::GeneticCode=standard_genetic_code, 
    allow_ambiguous_codons=true, 
    alternative_start=false)
ntseq = copy(ntseq)
translate(ntseq; code, allow_ambiguous_codons, alternative_start)
end

Base.sort(v::Vector{<:ORF}) = sort(v, by = _orf_sort_key)

function _orf_sort_key(orf::ORF)
    return (orf.location, orf.strand)
end

function Base.show(io::IO, dtcm::DTCM)
    nucleotides = sort(collect(keys(dtcm.order)))

    # Print type
    println(io, "GeneFinder.DTCM{$(typeof(dtcm.order)), $(typeof(dtcm.counts)):")

    # Print header
    println(io, "   ", join(nucleotides, "  "))

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = dtcm.counts[i, :]
        row_str = join(row, "  ")
        println(io, "$nucleotide1  $row_str")
    end
end

function Base.show(io::IO, dtpm::DTPM)
    nucleotides = sort(collect(keys(dtpm.order)))

    # Print type
    println(io, "GeneFinder.DTPM{$(typeof(dtpm.order)), $(typeof(dtpm.probabilities)):")

    # Print header
    println(io, "   ", join(nucleotides, "      "))

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = dtpm.probabilities[i, :]
        row_str = join(row, "  ")
        println(io, "$nucleotide1  $row_str")
    end
end