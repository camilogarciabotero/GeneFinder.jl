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

# Values
- `PSTRAND = 1`: Positive/forward strand (+)
- `NSTRAND = 2`: Negative/reverse strand (-)

# Example
```julia
strand = PSTRAND  # Positive strand
strand = NSTRAND  # Negative strand
```
"""
@enum Strand::Int begin
    PSTRAND = 1  # Positive/forward strand (+)
    NSTRAND = 2  # Negative/reverse strand (-)
end

# Helper function for validation
function _isvalidorf(seqid::Symbol, range::UnitRange{Int64}, strand::Strand, frame::Int8, features::NamedTuple)::Bool

    # Sanity check: strand validity
    strand in (PSTRAND, NSTRAND) ||
        throw(ArgumentError("Invalid strand: $(strand), expected PSTRAND (+) or NSTRAND (-)"))

    # Sanity check: seqid definition
    isdefined(Main, seqid) ||
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
    struct OpenReadingFrame{F<:GeneFinderMethod}

The `OpenReadingFrame` (aliased as `ORF`) struct represents an Open Reading Frame in genomics.

# Type Parameter
- `F<:GeneFinderMethod`: The gene finding algorithm used to identify this ORF.

# Fields
- `seqid::Symbol`: The identifier of the source sequence to which the ORF belongs.
- `range::UnitRange{Int64}`: The position range (start:stop) of the ORF on the sequence.
- `strand::Strand`: The strand orientation (`PSTRAND` or `NSTRAND`).
- `frame::Int8`: The reading frame (1, 2, or 3).
- `features::NamedTuple`: Additional features/metadata associated with the ORF.

# Validation
The constructor validates:
- Strand must be `PSTRAND` or `NSTRAND`
- Frame must be 1, 2, or 3
- Range length must be divisible by 3
- Range length must be ≥ 6 (minimum for start + stop codons)
- Features must be a `NamedTuple`

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Create an ORF with NaiveFinder method
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (;))

# Create an ORF with additional features
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (score=0.95, gc_content=0.52))
```

See also: [`sequence`](@ref), [`features`](@ref), [`strand`](@ref), [`frame`](@ref)
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
    sequence(orf::ORF{F}) where {F}

Extract the DNA sequence corresponding to the given Open Reading Frame (ORF).

# Arguments
- `orf::ORF{F}`: The ORF for which to extract the DNA sequence.

# Returns
- `LongSubSeq{DNAAlphabet{4}}`: For positive strand ORFs, returns a view (no allocation).
- `LongDNA{4}`: For negative strand ORFs, returns the reverse complement (allocates).

# Behavior
- For `PSTRAND`: Returns a subsequence view of the source sequence.
- For `NSTRAND`: Returns the reverse complement of the subsequence.

# Validation (bounds checked)
- Verifies the sequence starts with a valid start codon (ATG).
- Verifies the sequence ends with a valid stop codon (TAA, TAG, or TGA).
- Warns if frame is inconsistent with the start position.

# Example
```julia
seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (;))

dna_seq = sequence(orf)  # Returns the DNA sequence for this ORF
```

See also: [`source`](@ref), [`ORF`](@ref)
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
    source(orf::ORF{F}) where {F}

Retrieve the source sequence associated with the given ORF.

# Arguments
- `orf::ORF{F}`: The ORF for which to retrieve the source sequence.

# Returns
- The source sequence object referenced by the ORF's `seqid` field.

# Note
The source sequence must be defined in the `Main` module scope. If the sequence
identifier does not correspond to a defined variable, an error will be thrown.

# Example
```julia
seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (;))

src = source(orf)  # Returns the `seq` variable
```

See also: [`sequence`](@ref), [`seqid`](@ref)
"""
function source(i::ORF{F}) where {F}
    return getfield(Main, i.seqid)
end

"""
    features(orf::ORF{F}) where {F}

Extract the features/metadata from an ORF.

# Arguments
- `orf::ORF{F}`: The ORF from which to extract features.

# Returns
- `NamedTuple`: The features associated with the ORF (may be empty).

# Example
```julia
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (score=0.95, gc=0.52))

feats = features(orf)  # Returns (score = 0.95, gc = 0.52)
feats.score            # Access individual feature: 0.95
```

See also: [`ORF`](@ref)
"""
function features(i::ORF{F}) where {F}
    return i.features
end

"""
    leftposition(orf::ORF{F}) where {F}

Get the left (start) position of the ORF range.

# Arguments
- `orf::ORF{F}`: The ORF to query.

# Returns
- `Int`: The first position of the ORF range.

# Example
```julia
orf = ORF{NaiveFinder}(:seq, 10:42, PSTRAND, Int8(1), (;))
leftposition(orf)  # Returns 10
```

See also: [`rightposition`](@ref), [`ORF`](@ref)
"""
function leftposition(i::ORF{F}) where {F}
    return Int(Base.first(i.range))
end

"""
    rightposition(orf::ORF{F}) where {F}

Get the right (end) position of the ORF range.

# Arguments
- `orf::ORF{F}`: The ORF to query.

# Returns
- `Int`: The last position of the ORF range.

# Example
```julia
orf = ORF{NaiveFinder}(:seq, 10:42, PSTRAND, Int8(1), (;))
rightposition(orf)  # Returns 42
```

See also: [`leftposition`](@ref), [`ORF`](@ref)
"""
function rightposition(i::ORF{F}) where {F}
    return Int(Base.last(i.range))
end

"""
    frame(orf::ORF{F}) where {F}

Get the reading frame of the ORF.

# Arguments
- `orf::ORF{F}`: The ORF to query.

# Returns
- `Int8`: The reading frame (1, 2, or 3).

# Example
```julia
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(2), (;))
frame(orf)  # Returns 2
```

See also: [`strand`](@ref), [`ORF`](@ref)
"""
function frame(i::ORF{F}) where {F}
    return i.frame
end

"""
    strand(orf::ORF{F}) where {F}

Get the strand orientation of the ORF.

# Arguments
- `orf::ORF{F}`: The ORF to query.

# Returns
- `Strand`: Either `PSTRAND` (positive/forward) or `NSTRAND` (negative/reverse).

# Example
```julia
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (;))
strand(orf)  # Returns PSTRAND
```

See also: [`frame`](@ref), [`Strand`](@ref), [`ORF`](@ref)
"""
function strand(i::ORF{F}) where {F}
    return i.strand
end

"""
    seqid(orf::ORF{F}) where {F}

Get the sequence identifier of the ORF.

# Arguments
- `orf::ORF{F}`: The ORF to query.

# Returns
- `Symbol`: The identifier of the source sequence, or `:unnamedseq` if not set.

# Example
```julia
orf = ORF{NaiveFinder}(:chromosome1, 1:33, PSTRAND, Int8(1), (;))
seqid(orf)  # Returns :chromosome1
```

See also: [`source`](@ref), [`ORF`](@ref)
"""
function seqid(i::ORF{F}) where {F}
    if i.seqid === nothing #|| i.seqid == "unnamedseq"
        return :unnamedseq
    else
        return i.seqid
    end
end

"""
    finder(orf::ORF{F}) where {F}

Get the gene finding method type used to identify this ORF.

# Arguments
- `orf::ORF{F}`: The ORF to query.

# Returns
- `Type{F}`: The gene finder method type (e.g., `NaiveFinder`).

# Example
```julia
orf = ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (;))
finder(orf)  # Returns NaiveFinder
```

See also: [`GeneFinderMethod`](@ref), [`ORF`](@ref)
"""
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