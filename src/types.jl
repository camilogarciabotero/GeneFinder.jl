export ORF, OpenReadingFrame
export Strand, GeneFinderMethod, Strand, PSTRAND, NSTRAND
export features, sequence, source, finder, frame, strand, seqid, leftposition, rightposition

"""
    abstract type GeneFinderMethod

Abstract base type for different ORF finding methods/algorithms.

Subtypes should implement the calling interface to find ORFs in a sequence.
"""
abstract type GeneFinderMethod end

# abstract type GeneFinderMethod end # This should've been defined here

"""
    Strand

An enumeration type representing DNA strand orientation.
"""
@enum Strand::Int begin
    PSTRAND = 1  # Positive/forward strand (+)
    NSTRAND = 2  # Negative/reverse strand (-)
end

# # Helper function to convert Char strand to Strand enum
# function _to_strand(s::Strand)::Strand
#     return s
# end

# function _to_strand(s::Char)::Strand
#     s == '+' && return PSTRAND
#     s == '-' && return NSTRAND
#     throw(ArgumentError("Invalid strand character '$s'; expected '+' or '-'"))
# end

# Helper function for validation
function _isvalidorf(seqid::Symbol, range::UnitRange{Int64}, strand::Strand, frame::Int8, features::NamedTuple)::Bool

    # Sanity check: strand validity
    strand in (PSTRAND, NSTRAND) ||
        throw(ArgumentError("Invalid strand: $(strand), expected PSTRAND (+) or NSTRAND (-)"))

    # Sanity check: seqid definition
    isdefinedglobal(Main, seqid) ||
        @warn "The source sequence '$(seqid)' is not defined. Make sure to define it or supply the correct source sequence identifier."

    # Sanity check: frame validity
    (1 ≤ frame ≤ 3) || 
        throw(ArgumentError("Invalid frame: $(frame), expected 1, 2, or 3"))

    orflen = length(range)
    # Sanity check: range length divisible by 3
    (orflen % 3 == 0) || 
        throw(ArgumentError("ORF length ($(orflen)) is not divisible by 3, incomplete codons detected"))

    # Sanity check: minimum length (start + stop codons)
    orflen ≥ 6 ||
        throw(ArgumentError("ORF sequence too short to contain start and stop codons"))

    # Sanity check: features must be a NamedTuple
    features isa NamedTuple ||
        throw(ArgumentError("Features must be a NamedTuple, got $(typeof(features))"))

    return true
end

"""
    struct ORF{F}

The `ORF` struct represents an Open Reading Frame (ORF) in genomics.

# Fields
- `seqid::Symbol`: The identifier of the sequence to which the ORF belongs.
- `range::UnitRange{<:Int64}`: The position range of the ORF on the sequence.
- `strand::Strand`: The strand on which the ORF is located.
- `frame::Int8`: The reading frame of the ORF (1, 2, or 3).
- `features::NamedTuple`: The features associated with the ORF.

# Example

```julia
ORF{NaiveFinder}(:seq01, 1:33, PSTRAND, Int8(1), (;score = 0.8))
```
"""
struct OpenReadingFrame{F<:GeneFinderMethod}
    seqid::Symbol
    range::UnitRange{Int64}
    strand::Strand
    frame::Int8
    features::NamedTuple

    function OpenReadingFrame{F}(
        seqid::Symbol,
        range::UnitRange{Int64},
        strand::Strand,
        frame::Int8,
        features::NamedTuple = (;)
    ) where {F<:GeneFinderMethod}
        _isvalidorf(seqid, range, strand, frame, features)
        return new{F}(seqid, range, strand, frame, features)
    end
end

const ORF = OpenReadingFrame

const START = dna"ATG"
const STOPS = (dna"TAA", dna"TAG", dna"TGA")

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

    @views sub = src[i.range]  # view; no allocation

    # s ∉ (PSTRAND, NSTRAND) && throw(ArgumentError("Cannot extract sequence for strand $(s); expected PSTRAND (+) or NSTRAND (-)"))
    seq = i.strand === PSTRAND ? sub : reverse_complement(sub)  # reverse_complement allocates (can't be a view)

    @boundscheck begin
        # Frame consistency check
        fr = frame(i)
        if i.strand === PSTRAND
            exp = mod1(leftposition(i), 3)
            fr == exp || 
                @warn "Frame $(fr) may be inconsistent with start position $(leftposition(i)), expected frame $exp"
        end

        @views begin
        # Codon checks

            # Start codon check
            seq[1:3] == START || 
            throw(ArgumentError("Invalid start codon $(seq[1:3]), expected ATG"))

            # Stop codon check:
            seq[end-2:end] in STOPS ||
                throw(ArgumentError("Invalid stop codon $(seq[end-2:end]), expected TAA, TAG, or TGA"))
        end
    end

    return seq
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
    return getfield(Main, i.seqid)
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

function seqid(i::ORF{F}) where {F}
    if i.seqid === nothing #|| i.seqid == "unnamedseq"
        return :unnamedseq
    else
        return i.seqid
    end
end

finder(i::ORF{F}) where {F} = F

function Base.show(io::IO, i::ORF{F}) where {F}
    print(io, "ORF{", finder(i), "}(", leftposition(i), ":", rightposition(i), ", '")
    show(io, strand(i))
    print(io, "', ", frame(i))
    feats = features(i)
    if !isempty(feats)
        if length(feats) <= 3
            print(io, ", ", feats)
        else
            keys_subset = keys(feats)[1:3]
            vals_subset = Tuple(feats[k] for k in keys_subset)
            truncated = NamedTuple{keys_subset}(vals_subset)
            print(io, ", ", truncated, "…")
        end
    end
    print(io, ")")
end


## Ideas for Gene struct

# struct CDS <: AbstractGene
#     orf::ORF
#     join::Bool
#     coding::Bool
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