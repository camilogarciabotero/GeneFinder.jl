# Methods from main packages that expand their fuctions to this package structs
import Base: length, iterate, sort, getindex

Base.sort(v::Vector{<:ORF}) = sort(v, by = _orf_sort_key)
Base.getindex(sequence::NucleicSeqOrView{A}, orf::ORF) where {A} = orf.strand == '+' ? (@view sequence[orf.location]) : reverse_complement(@view sequence[orf.location])
function _orf_sort_key(orf::ORF)
    return (orf.location, orf.strand)
end