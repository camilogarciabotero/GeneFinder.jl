# Methods from main packages that expand their fuctions to this package structs
import Base: show, length, iterate, sort

Base.sort(v::Vector{<:ORF}) = sort(v, by = _orf_sort_key)

function _orf_sort_key(orf::ORF)
    return (orf.location, orf.strand)
end