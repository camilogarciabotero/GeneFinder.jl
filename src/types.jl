import GenomicFeatures: first, last, length, strand, groupname
import FASTX: sequence

export RBS, ORF
export gc, features, sequence
export groupname, finder, frame, scheme, score

#### Ribosome Binding Site (RBS) struct ####

# seq[orf.first-bin01.offset.start:orf.first-1]

# motifs = [dna"GGA", dna"GAG", dna"AGG"]
struct RBS
    motif::BioRegex{DNA}
    offset::UnitRange{Int64} # offset
    score::Float64

    function RBS(motif::BioRegex{DNA}, offset::UnitRange{Int64}, score::Float64)
        return new(motif, offset, score)
    end
end

# Structs associated with gene models

# const FEATUREDICT = Dict{Symbol,Any}(:gc => 0.0, :length => 0, :score => 0.0)
# const Feature = Union{Real, Vector{RBS}}
struct ORF{N,F} <: GenomicFeatures.AbstractGenomicInterval{F}
    groupname::String
    first::Int64
    last::Int64
    strand::Strand
    frame::Int
    seq::LongSubSeq{DNAAlphabet{N}}
    features::Dict{Symbol,Any}
    scheme::Union{Nothing,Function}
end

function ORF{N,F}(
    ::Type{F}, #finder
    groupname::String,
    first::Int64,
    last::Int64,
    strand::Union{Strand,Char},
    frame::Int,
    seq::LongSubSeq{DNAAlphabet{N}},
    features::Dict{Symbol,Any},
    scheme::Union{Nothing,Function}=nothing
) where {N,F<:GeneFinderMethod}
    # @assert frame in (1, 2, 3) && location[1] < location[end] && score isa Union{Nothing, Float64} && length(seq) == last(location) - first(location) + 1 && length(seq) % 3 == 0 "Invalid input. Please check the frame, location, score, and sequence length."
    return ORF{N,F}(groupname, first, last, strand, frame, seq, features, scheme) #finder
end

function groupname(i::ORF{N,F}) where {N,F}
    return i.groupname
end

function first(i::ORF{N,F}) where {N,F}
    return i.first
end

function last(i::ORF{N,F}) where {N,F}
    return i.last
end

function frame(i::ORF{N,F}) where {N,F}
    return i.frame
end

finder(i::ORF{N,F}) where {N,F} = F

function scheme(i::ORF{N,F}) where {N,F}
    return i.scheme
end

function sequence(i::ORF{N,F}) where {N,F}
    return i.seq #i.strand == STRAND_POS ? i.seq : reverse_complement(i.seq)
end

function score(i::ORF{N,F}) where {N,F}
    return i.features[:score]
end

function length(i::ORF{N,F}) where {N,F}
    return length(i.seq)
end

function features(i::ORF{N,F}) where {N,F}
    return i.features
end

function Base.show(io::IO, i::ORF{N,F}) where {N,F}
    if get(io, :compact, false)
        print(io, "ORF{$(finder(i))}($(leftposition(i)):$(rightposition(i)), '$(strand(i))', $(frame(i)))") #{$(typeof(finder(i)))} $(score(i))
    else
        print(io, "ORF{$(finder(i))}($(leftposition(i)):$(rightposition(i)), '$(strand(i))', $(frame(i)))") # , $(score(i))
    end
end

## Ideas for Gene struct

# struct CDS <: AbstractGene
#     orf::ORF
#     coding::Bool
# end


# frm = 1
# frst = 1
# lst = 99
# str = '+'
# seq = randdnaseq(99)

# orf = ORF{4,NaiveFinder}("seq1", frst, lst, str, frm, seq, nothing, nothing, nothing)