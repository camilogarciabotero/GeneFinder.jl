import GenomicFeatures: first, last, length, strand, groupname, metadata
import FASTX: sequence

export Features, RBS, ORF
export features, sequence, source
export groupname, finder, frame, scheme, score, strand, STRAND_BOTH, STRAND_NEG, STRAND_POS, STRAND_NA

#### Ribosome Binding Site (RBS) struct ####

# seq[orf.first-bin01.offset.start:orf.first-1]

# motifs = [dna"GGA", dna"GAG", dna"AGG"]

struct Features{K}
    fts::K
end

Features(fts::NamedTuple) = Features{typeof(fts)}(fts)
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
# features::Dict{Symbol,Any} or # features::@NamedTuple{score::Float64, rbs::Any}
# seq::LongSubSeq{DNAAlphabet{N}}
struct ORF{N,F} <: GenomicFeatures.AbstractGenomicInterval{F}
    groupname::String
    first::Int64
    last::Int64
    strand::Strand
    frame::Int
    features::Features
    scheme::Union{Nothing,Function}

    function ORF{N,F}(
        groupname::String,
        first::Int64,
        last::Int64,
        strand::Strand,
        frame::Int,
        features::Features,
        scheme::Union{Nothing,Function}
    ) where {N,F}
        @assert frame in (1, 2, 3) "Invalid frame. Please provide a frame between 1 and 3."

        return new{N,F}(groupname, first, last, strand, frame, features, scheme)
    end
end

function ORF{N,F}(
    ::Type{F}, #finder
    groupname::String,
    first::Int64,
    last::Int64,
    strand::Strand,
    frame::Int,
    features::Features, # ::Dict{Symbol,Any} or # ::@NamedTuple{score::Float64, rbs::Any} or @NamedTuple{Vararg{typeof(...)}}
    scheme::Union{Nothing,Function}=nothing
) where {N,F<:GeneFinderMethod}
    return ORF{N,F}(groupname, first, last, strand, frame, features, scheme) #finder seq
end


# function ORF{N,F}(
#     ::Type{F}, #finder
#     groupname::Union{Nothing,String},
#     first::Int64,
#     last::Int64,
#     strand::Strand,
#     frame::Int,
#     features::Features, # ::Dict{Symbol,Any} or # ::@NamedTuple{score::Float64, rbs::Any} or @NamedTuple{Vararg{typeof(...)}}
#     scheme::Union{Nothing,Function}=nothing
# ) where {N,F<:GeneFinderMethod}
#     groupname = groupname === nothing ? "seq" : groupname
#     return ORF{N,F}(groupname, first, last, strand, frame, features, scheme) #finder seq
# end

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
    seqsymb = Symbol(i.groupname)
    seq = getfield(Main, seqsymb)
    return i.strand == STRAND_POS ? @view(seq[i.first:i.last]) : reverse_complement(@view(seq[i.first:i.last])) #seq[i.first:i.last]
end

function source(i::ORF{N,F}) where {N,F}
    seqsymb = Symbol(i.groupname)
    return getfield(Main, seqsymb)
end

function score(i::ORF{N,F}) where {N,F}
    return i.features.fts[:score]
end

function length(i::ORF{N,F}) where {N,F}
    return length(sequence(i))
end

function features(i::ORF{N,F}) where {N,F}
    return i.features.fts
end

function strand(i::ORF{N,F}) where {N,F}
    return i.strand
end

function Base.show(io::IO, i::ORF{N,F}) where {N,F}
    if get(io, :compact, false)
        print(io, "ORF{$(finder(i))}($(leftposition(i)):$(rightposition(i)), '$(strand(i))', $(frame(i)))") #{$(typeof(finder(i)))} $(score(i))
    else
        print(io, "ORF{$(finder(i))}($(leftposition(i)):$(rightposition(i)), '$(strand(i))', $(frame(i)))") # , $(score(i))
    end
end

function metadata(i::ORF{N,F}) where {N,F}
    return features(i)
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