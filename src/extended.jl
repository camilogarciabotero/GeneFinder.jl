# Methods from main packages that expand their fuctions to this package structs
import Base: isless, iterate, sort, getindex, length
import BioSequences: translate

Base.isless(a::ORF{N,F}, b::ORF{N,F}) where {N,F} = isless(a.first:a.last, b.first:b.last)

function getindex(seq::NucleicSeqOrView{A}, orf::ORF{N,F}) where {A,N,F}
    @assert @view(seq[begin:end]) == source(orf) "The source sequence of the ORF and the given sequence are different"

    if orf.strand == STRAND_POS
        s = @view seq[orf.first:orf.last]
    else
        s = reverse_complement(@view seq[orf.first:orf.last])
    end

    if !occursin(biore"T(AG|AA|GA)"dna, s)
        error("There is no stop codon at the end of the sequence/orf")
    end

    return s
end

function translate(orf::ORF{N,F}; kwargs...) where {N,F}
    return translate(sequence(orf); kwargs...)
end