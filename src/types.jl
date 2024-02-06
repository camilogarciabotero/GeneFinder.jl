# Structs associated with gene models 
abstract type AbstractGene end

"""
    struct GeneFeatures
        seqname::String
        start::Int64
        stop::Int64
        strand::Char
        frame::Int8
        score::Union{Float64, Nothing} # Add score field
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
        score::Union{Float64, Nothing} # Add score field
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: 

- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence
- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
- `frame`: is an Int type indicating the reading frame of the ORF. The frame is the position of the first nucleotide of the codon that starts the ORF, relative to the start of the sequence. It can be 1, 2, or 3.
- `score`: is a Union{Float64, Nothing} type indicating the score of the ORF. It can be a Float64 or nothing.
"""
struct ORF <: AbstractGene
    location::UnitRange{Int64}  # Note that it is also called position for gene struct in GenomicAnotations
    strand::Char
    frame::Int # Use Int64 instead of Int
    score::Union{Float64, Nothing} # Add score field

    function ORF(location::UnitRange{Int64}, strand::Char, frame::Int, score::Union{Float64, Nothing} = nothing)
        @assert frame in (1, 2, 3) "Invalid frame value. Frame must be 1, 2, or 3."
        @assert strand in ('+', '-') "Invalid strand value. Strand must be '+' or '-'."
        @assert location[1] < location[end] "Invalid location. Start must be less than stop."
        @assert score isa Union{Float64, Nothing} "Invalid score value. Score must be a Float64 or nothing."
    
        return new(location, strand, frame, score)
    end
end

