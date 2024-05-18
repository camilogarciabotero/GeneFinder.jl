export ORF
export NaiveFinder, NaiveFinderScored
export NaiveScoringScheme
# Structs associated with gene models 
abstract type AbstractGene end

"""
    struct GeneFeatures
        seqname::String
        start::Int64
        stop::Int64
        strand::Char
        frame::Int8
        score::Union{Nothing, Float64} # Add score field
        attribute::Dict
    end

This is the main Gene struct, based on the fields that could be found in a GFF3, still needs to be defined correctly,
    The idea is correct the frame and attributes that will have something like a possible list (id=Char;name=;locus_tag).
    The `write` and `get` functions should have a dedicated method for this struct.
"""
# struct GeneFeatures <: AbstractGene
#     seqname::String
#     start::Int64
#     stop::Int64
#     strand::Char
#     frame::Int8 # But maybe a Union to allow empty when reading a GFF? 
#     score::Union{Float64, Nothing} # Add score field
#     attribute::Dict # Should be a Dict perhaps or a NamedTuple
# end

# struct GeneFeatures <: AbstractGene
#     seqname::String could be a Union{String, Nothing} and named id.
#     orf::ORF
#     attribute::Dict
# end

"""
    struct ORF
        location::UnitRange{Int64}
        strand::Char
        frame::Int
        score::Union{Nothing, Float64}
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: 

- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence
- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
- `frame`: is an Int type indicating the reading frame of the ORF. The frame is the position of the first nucleotide of the codon that starts the ORF, relative to the start of the sequence. It can be 1, 2, or 3.
- `score`: is a Union{Nothing, Float64} type indicating the score of the ORF. It can be a Float64 or nothing.
"""
# struct ORF <: AbstractGene # could also be mutable so that scores can be updated later, but tests fails
#     #TODOs:  location might be a complex Union allowing UnitRange or a Join of ranges (e.g. 1..100, 200..300)?
#     location::UnitRange{Int64}  # Note that it is also called position for gene struct in GenomicAnotations
#     strand::Char
#     frame::Int # Use Int64 instead of Int
#     score::Union{Nothing, Float64} # Add score field shold be a namedtuple with the score and the method used to score
#     # rbs::Union{Nothing, Vector{UnitRange{Int64}}} # Add RBS field  see https://github.com/deprekate/PHANOTATE/blob/c77e80caef8dc3264f7dc698b087bcd486216bcb/phanotate_modules/functions.py#L48

#     function ORF(location::UnitRange{Int64}, strand::Char, frame::Int, score::Union{Nothing, Float64} = nothing)
#         @assert frame in (1, 2, 3) "Invalid frame value. Frame must be 1, 2, or 3."
#         @assert strand in ('+', '-') "Invalid strand value. Strand must be '+' or '-'."
#         @assert location[1] < location[end] "Invalid location. Start must be less than stop."
#         @assert score isa Union{Nothing, Float64} "Invalid score value. Score must be a Float64 or nothing."
    
#         return new(location, strand, frame, score)
#     end
# end

## Interfaces for gene finders

abstract type GeneFinderMethod end
struct NaiveFinder <: GeneFinderMethod end
struct NaiveFinderScored <: GeneFinderMethod end


## Interfaces for gene scoring schemes

abstract type GeneScoringScheme end
struct NaiveScoringScheme <: GeneScoringScheme end


export ORF

struct ORF{F} <: GenomicFeatures.AbstractGenomicInterval{F}
    groupname::String
    first::Int64
    last::Int64
    strand::Strand
    frame::Int
    finder::F
    scheme::Union{Nothing, Function}
    score::Union{Nothing, Float64}
    rbs::Union{Nothing, Vector{UnitRange{Int64}}} # Add RBS field 
end

function ORF(
    groupname::String,
    first::Int64,
    last::Int64,
    strand::Union{Strand,Char},
    frame::Int,
    finder::F,
    scheme::Union{Nothing, Function}=nothing,
    score::Union{Nothing, Float64}=nothing,
    rbs::Union{Nothing, Vector{UnitRange{Int64}}}=nothing
) where {F <: GeneFinderMethod}
    return ORF{F}(groupname, first, last, strand, frame, finder, scheme, score, rbs)
end

function ORF(
    id::String,
    location::UnitRange{Int64},
    strand::Union{Strand,Char},
    frame::Int, finder::F,
    scheme::Union{Nothing, Function}=nothing,
    score::Union{Nothing, Float64}=nothing,
    rbs::Union{Nothing, Vector{UnitRange{Int64}}}=nothing
) where {F <: GeneFinderMethod}
    return ORF{F}(id, first(location), last(location), strand, frame, finder, scheme, score, rbs)
end

function frame(i::ORF)
    return i.frame
end

function scheme(i::ORF)
    return i.scheme
end

function score(i::ORF)
    return i.score
end

function finder(i::ORF)
    return i.finder
end

function Base.show(io::IO, i::ORF)
    if get(io, :compact, false)
        # print(io, groupname(i), ":", leftposition(i), "-", rightposition(i))
        print(io, "ORF{$(typeof(finder(i)))}($(groupname(i)), $(leftposition(i)):$(rightposition(i)), $(strand(i)), $(frame(i)), $(score(i)))")
    else
        # print(io, "ORF{$(typeof(finder(i))),$(scheme(i))}($(groupname(i)), $(leftposition(i)):$(rightposition(i)), $(strand(i)), $(frame(i)), $(round(score(i), digits=5)))")
        print(io, "ORF{$(typeof(finder(i)))}($(groupname(i)), $(leftposition(i)):$(rightposition(i)), $(strand(i)), $(frame(i)), $(score(i)))")
        # println(io, " ", groupname(i), " ", leftposition(i), ":", rightposition(i), ", ", strand(i), ", ", frame(i), ", score: $(round(score(i), digits=2))")
    end
end