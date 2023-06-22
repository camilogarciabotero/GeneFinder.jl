# Methods from main packages that expand their fuctions to this package structs
import BioSequences: translate, standard_genetic_code
import Base: show

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

function Base.show(io::IO, tcm::TCM)
    nucleotides = sort(collect(keys(tcm.order)))

    # Print type
    println(io, "GeneFinder.TCM{$(typeof(tcm.order)), $(typeof(tcm.counts)):")

    # Print header
    println(io, "   ", join(nucleotides, "  "))

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = tcm.counts[i, :]
        row_str = join(row, "  ")
        println(io, "$nucleotide1  $row_str")
    end
end

function Base.show(io::IO, tpm::TPM)
    nucleotides = sort(collect(keys(tpm.order)))

    # Print type
    println(io, "GeneFinder.tpm{$(typeof(tpm.order)), $(typeof(tpm.probabilities)):")

    # Print header
    println(io, "   ", join(nucleotides, "      "))

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = tpm.probabilities[i, :]
        row_str = join(row, "  ")
        println(io, "$nucleotide1  $row_str")
    end
end

function show(io::IO, model::TransitionModel)
    coding_rows, coding_cols = size(model.coding)
    noncoding_rows, noncoding_cols = size(model.noncoding)

    println(io, "TransitionModel with:")

    println(io, "  Coding transition probability matrix:")
    for i in 1:coding_rows
        for j in 1:coding_cols
            print(io, "    $(model.coding[i, j])")
        end
        println(io)
    end

    println(io, "  Noncoding transition ptobability matrix:")
    for i in 1:noncoding_rows
        for j in 1:noncoding_cols
            print(io, "    $(model.noncoding[i, j])")
        end
        println(io)
    end

    println(io, "  Coding initial probabilities: $(model.codinginits)")
    println(io, "  Noncoding initial probabilities: $(model.noncodinginits)")
end
