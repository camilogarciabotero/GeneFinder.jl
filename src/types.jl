# Structs associated with gene models 
abstract type AbstractGene end

"""
    struct GeneFeatures
        seqname::String
        start::Int64
        stop::Int64
        score::Float64
        strand::Char
        frame::'.' 
        attribute::
    end

This is the main Gene struct, based on the fields that could be found in a GFF3, still needs to be defined correctly,
    The idea is correct the frame and attributes that will have something like a possible list (id=Char;name=;locus_tag).
    The `write` and `get` functions should have a dedicated method for this struct.
"""
# struct GeneFeatures <: AbstractGene
#     seqname::String
#     start::Int64
#     stop::Int64
#     score::Float64
#     strand::Char
#     frame::Int8 # But maybe a Union to allow empty when reading a GFF? 
#     attribute::Dict # Should be a Dict perhaps
# end

"""
    struct ORF
        location::UnitRange{Int64}
        strand::Char
        frame::Int8
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: 

- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence
- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
"""
struct ORF <: AbstractGene
    location::UnitRange{Int64}  # Note that it is also called position for gene struct in GenomicAnotations
    strand::Char
    frame::Int # Use Int64 instead of Int

    function ORF(location::UnitRange{Int64}, strand::Char, frame::Int)
        @assert frame in (1, 2, 3) "Invalid frame value. Frame must be 1, 2, or 3."
        @assert strand in ('+', '-') "Invalid strand value. Strand must be '+' or '-'."
        new(location, strand, frame)
    end
end