# Methods from main packages that expand their fuctions to this package structs
import Base: isless, iterate, sort, getindex


# Base.isless(a::ORF, b::ORF) = isless(a.location, b.location)
Base.isless(a::ORF, b::ORF) = isless(a.first:a.last, b.first:b.last)
# Base.isless(a::ORF, b::ORF) = isless(a.location, b.location)
# Base.isless(a::ORF, b::ORF) = isless(a.score, b.score)

# Base.sort(v::Vector{<:ORF}; kwargs...) = sort(v, by = _orf_sort_key)
# Base.sort!(v::Vector{<:ORF}; kwargs...) = sort!(v, by = _orf_sort_key)


#TODOs: how to make more robust the getindex method? confroning the frames?
# Base.getindex(sequence::NucleicSeqOrView{A}, orf::ORF) where {A} = orf.strand == '+' ? (@view sequence[orf.location]) : reverse_complement(@view sequence[orf.location])

# function getindex(sequence::NucleicSeqOrView{A}, orf::ORF) where {A}
#     if orf.strand == '+'
#         return @view sequence[orf.location]
#     else
#         return reverse_complement(@view sequence[orf.location])
#     end
# end

function getindex(sequence::NucleicSeqOrView{A}, orf::ORF) where {A}
    if orf.strand == '+' || orf.strand == STRAND_POS
        return @view sequence[orf.first:orf.last]
    else
        return reverse_complement(@view sequence[orf.first:orf.last])
    end
end

# function _orf_sort_key(orf::ORF)
#     return (orf.location, orf.strand, orf.score)
# end