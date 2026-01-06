import GenomicFeatures: first, last, strand, groupname, metadata

export ORFI, OpenReadingFrameInterval
export features, sequence, source, finder, frame
export groupname, strand, STRAND_BOTH, STRAND_NEG, STRAND_POS, STRAND_NA

# abstract type GeneFinderMethod end # This should've been defined here

"""
    Strand

An enumeration type representing DNA strand orientation.
"""

@enum Strand::Int8 begin
    PSTRAND = 1  # Positive/forward strand (+)
    NSTRAND = 2  # Negative/reverse strand (-)
end

function Base.show(io::IO, s::Strand)
    if s === PSTRAND
        print(io, "+")
    elseif s === NSTRAND
        print(io, "-")
    else
        print(io, "unknown")
    end
end

"""
    struct ORF{F}

The `ORF` struct represents an Open Reading Frame (ORF) in genomics.

# Fields
- `seqid::Symbol`: The identifier of the sequence to which the ORF belongs.
- `range::UnitRange{Int32}`: The position range of the ORF on the sequence.
- `strand::Strand`: The strand on which the ORF is located.
- `frame::Int8`: The reading frame of the ORF (1, 2, or 3).
- `features::NamedTuple`: The features associated with the ORF.

# Example

```julia
ORF{NaiveFinder}(:seq01, 1:33, STRAND_POS, 1, (score = 0.8,))
```
"""
struct OpenReadingFrame{F<:GeneFinderMethod}
    seqid::Symbol
    range::UnitRange{Int32} # Could be more complex for Introns? But for ORFs it's fine, maybe for the more general GeneInterval
    strand::Strand
    frame::Int8
    features::NamedTuple
end

function OpenReadingFrame(
    ::Type{F},
    seqid::Union{Symbol, String},
    range::UnitRange{<:Integer},
    strand::Strand,
    frame::Int8,
    features::NamedTuple = (;)
) where {F<:GeneFinderMethod}
    seqid_sym = isa(seqid, Symbol) ? seqid : Symbol(seqid)
    return ORF{F}(seqid_sym, Int32(Base.first(range)):Int32(Base.last(range)), strand, frame, features)
end

const ORF = OpenReadingFrame

"""
    sequence(i::ORF{F})

Extracts the DNA sequence corresponding to the given open reading frame (ORF).
Uses the source sequence referenced by the ORF's seqid.

For positive strand ORFs, returns a LongSubSeq view (avoiding unnecessary copying).
For negative strand ORFs, returns the reverse complement (requires allocation).

# Arguments
- `i::ORF{F}`: The open reading frame (ORF) for which the DNA sequence needs to be extracted.

# Returns
- A LongSubSeq for positive strand ORFs, or a reverse complement sequence for negative strand ORFs.
"""
@inline function sequence(i::ORF{F}) where {F}
    src = source(i)
    sub = @view src[i.range]

    if i.strand === PSTRAND
        return sub
    elseif i.strand === NSTRAND
        return convert(LongSubSeq, reverse_complement(sub))
    else
        throw(ArgumentError("Cannot extract sequence for strand $(i.strand); expected PSTRAND (+) or NSTRAND (-)"))
    end
end

"""
    source(i::ORF{F})

Get the source sequence associated with the given `ORF` object.

# Arguments
- `i::ORF{F}`: The `ORF` object for which to retrieve the source sequence.

# Returns
The source sequence associated with the `ORF` object.

!!! warning
    The `source` method works if the sequence is defined in the global scope. Otherwise it will throw an error.
"""
function source(i::ORF{F}) where {F}
    try
        return getfield(Main, i.seqid)
    catch e
        error("The source sequence of the ORF is defined as $(i.seqid). Make sure to either define it or supply the correct source sequence.")
    end
end

"""
    features(i::ORF{F})

Extracts the features from an `ORF` object.

# Arguments
- `i::ORF{F}`: An `ORF` object.

# Returns
The features of the `ORF` object.
"""
function features(i::ORF{F}) where {F}
    return i.features
end

function leftposition(i::ORF{F}) where {F}
    return Int(Base.first(i.range))
end

function rightposition(i::ORF{F}) where {F}
    return Int(Base.last(i.range))
end

function frame(i::ORF{F}) where {F}
    return i.frame
end

function strand(i::ORF{F}) where {F}
    return i.strand
end

finder(i::ORF{F}) where {F} = F


function Base.show(io::IO, i::ORFI{N,F}) where {N,F}
    if get(io, :compact, false)
        print(io, "ORFI{$(finder(i))}($(leftposition(i)):$(rightposition(i)), '$(strand(i))', $(frame(i)))") #{$(typeof(finder(i)))} $(score(i))
    else
        print(io, "ORFI{$(finder(i))}($(leftposition(i)):$(rightposition(i)), '$(strand(i))', $(frame(i)))") # , $(score(i))
    end
end

## Ideas for Gene struct

# struct CDS <: AbstractGene
#     orf::ORFI
#     coding::Bool
#     join::Bool
# end

# #### Ribosome Binding Site (RBS) struct ####

# struct RBS
#     motif::BioRegex{DNA}
#     offset::UnitRange{Int64} # offset
#     score::Float64

#     function RBS(motif::BioRegex{DNA}, offset::UnitRange{Int64}, score::Float64)
#         return new(motif, offset, score)
#     end
#     # rbsinst = RBS(biore"RRR"dna, 3:4, 1.0)
# end

# seq[orf.first-bin01.offset.start:orf.first-1]

# motifs = [dna"GGA", dna"GAG", dna"AGG"]

####### End of RBS struct #######