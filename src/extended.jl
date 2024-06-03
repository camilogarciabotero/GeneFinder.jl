# Methods from main packages that expand their fuctions to this package structs
import Base: isless, iterate, sort, getindex, length
import BioSequences: translate

Base.isless(a::ORF{N,F}, b::ORF{N,F}) where {N,F} = isless(a.first:a.last, b.first:b.last)

function getindex(sequence::NucleicSeqOrView{A}, orf::ORF{N,F}) where {A,N,F}
    if orf.strand == '+' || orf.strand == STRAND_POS
        return @view sequence[orf.first:orf.last]
    else
        return reverse_complement(@view sequence[orf.first:orf.last])
    end
end

function translate(orf::ORF{N,F}; kwargs...) where {N,F}
    return translate(sequence(orf); kwargs...)
end