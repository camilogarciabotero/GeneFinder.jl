export ORF, OpenReadingFrame, ORFCollection
export Strand, GeneFinderMethod, PSTRAND, NSTRAND
export features, sequence, finder, frame, strand, leftposition, rightposition
export source, orfvector

"""
    abstract type GeneFinderMethod

Abstract base type for different ORF finding methods/algorithms.

Subtypes should implement the calling interface to find ORFs in a sequence,
returning an `ORFCollection`.
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
    struct ORFCollection{F<:GeneFinderMethod, S<:NucleicSeqOrView}

A collection of Open Reading Frames (ORFs) bundled with a view of their source sequence.

This is the primary return type for all `GeneFinderMethod` implementations.
It provides a clean API for accessing ORFs and their corresponding sequences
without relying on global state.

The source is always stored as a view (`LongSubSeq`) to avoid unnecessary copying
while maintaining a reference to the original sequence data.

# Type Parameters
- `F<:GeneFinderMethod`: The gene finding algorithm used to identify these ORFs.
- `S<:NucleicSeqOrView`: The type of the source sequence view.

# Fields
- `source::S`: A view of the source DNA sequence containing the ORFs.
- `orfs::Vector{ORF{F}}`: The vector of ORFs found in the source sequence.

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAG"
collection = findorfs(seq, finder=NaiveFinder)

# Source is stored as a view
typeof(source(collection))  # LongSubSeq{DNAAlphabet{4}}

# Iteration
for orf in collection
    println(orf)
end

# Sequence extraction
orfseq = sequence(collection, 1)
```

See also: [`ORF`](@ref), [`sequence`](@ref), [`source`](@ref)
"""
struct ORFCollection{F<:GeneFinderMethod, S<:LongSubSeq}
    source::S
    orfs::Vector{ORF{F}}
    
    function ORFCollection(source::LongSubSeq{DNAAlphabet{N}}, orfs::Vector{ORF{F}}) where {N, F<:GeneFinderMethod}
        # Always store as a view for consistency and memory efficiency
        srcview = @view source[begin:end]
        
        # Validate all ORFs are within bounds
        seqlen = length(srcview)
        for orf in orfs
            rightposition(orf) <= seqlen || 
                throw(BoundsError(srcview, orf.range))
        end
        return new{F, typeof(srcview)}(srcview, orfs)
    end
end

# ────────────────────────────────────────────────────────────────────────────────
# ORFCollection Interface
# ────────────────────────────────────────────────────────────────────────────────

"""
    source(collection::ORFCollection)

Get the source sequence associated with an ORF collection.

# Arguments
- `orf::ORF{F}`: The ORF for which to extract the DNA sequence.

# Returns
- The source DNA sequence.

# Example
```julia
collection = findorfs(seq)
src = source(collection)  # Returns the original sequence
```
"""
source(collection::ORFCollection) = collection.source

"""
    orfvector(collection::ORFCollection)

Get the vector of ORFs from a collection.

# Arguments
- `collection::ORFCollection`: The collection to query.

# Returns
- `Vector{ORF{F}}`: The vector of ORFs.

# Example
```julia
collection = findorfs(seq)
orf_vector = orfvector(collection)
```
"""
orfvector(collection::ORFCollection) = collection.orfs

"""
    finder(collection::ORFCollection{F}) where {F}

Get the gene finding method type used for this collection.

# Returns
- `Type{F}`: The gene finder method type (e.g., `NaiveFinder`).
"""
finder(::ORFCollection{F}) where {F} = F
# orfvector(collection::ORFCollection{F}) where {F} = collection.orfs

# Iteration interface
Base.iterate(c::ORFCollection) = iterate(c.orfs)
Base.iterate(c::ORFCollection, state) = iterate(c.orfs, state)
Base.length(c::ORFCollection) = length(c.orfs)
Base.eltype(::Type{ORFCollection{F,S}}) where {F,S} = ORF{F}
Base.isempty(c::ORFCollection) = isempty(c.orfs)

# Indexing interface
Base.getindex(c::ORFCollection, i::Int) = c.orfs[i]
Base.getindex(c::ORFCollection, r::AbstractRange) = c.orfs[r]
Base.getindex(c::ORFCollection, v::AbstractVector{Bool}) = c.orfs[v]
Base.getindex(c::ORFCollection, v::AbstractVector{<:Integer}) = c.orfs[v]
Base.firstindex(c::ORFCollection) = firstindex(c.orfs)
Base.lastindex(c::ORFCollection) = lastindex(c.orfs)
Base.keys(c::ORFCollection) = keys(c.orfs)
Base.eachindex(c::ORFCollection) = eachindex(c.orfs)

# ────────────────────────────────────────────────────────────────────────────────
# Sequence Extraction
# ────────────────────────────────────────────────────────────────────────────────

"""
    sequence(collection::ORFCollection, i::Int)

Extract the DNA sequence for the ORF at index `i` in the collection.

# Arguments
- `collection::ORFCollection`: The collection containing the ORF and source sequence.
- `i::Int`: The index of the ORF.

# Returns
- `LongSubSeq{DNAAlphabet{4}}`: The DNA sequence corresponding to the ORF as a view.

# Example
```julia
collection = findorfs(seq)
orfseq = sequence(collection, 1)  # Get sequence of first ORF
```

# Example

For getting the sequence of all ORFs are several alternatives:

```julia
collection = findorfs(seq)
# Using a for loop with push!
orfseqs = Vector{LongSubSeq{DNAAlphabet{4}}}()
for orf in collection
    push!(orfseqs, sequence(collection, orf))
end

# Using broadcasting
orfseq = sequence.(Ref(collection), collection.orfs)

# Using list comprehension
orfseqs = [sequence(collection, orf) for orf in collection]

# Using map
orfseqs = map(orf -> sequence(collection, orf), collection)

# Using indices
orfseqs = [sequence(collection, i) for i in eachindex(collection)]

# Using map with indices
orfseqs = map(i -> sequence(collection, i), eachindex(collection))

# Using a generator expression
orfseqs = collect(sequence(collection, orf) for orf in collection)

# Using a generator with indices
orfseqs = collect(sequence(collection, i) for i in eachindex(collection))
```

See also: [`sequences`](@ref), [`ORFCollection`](@ref)
"""
@inline function sequence(collection::ORFCollection, i::Int; kwargs...)
    @boundscheck checkbounds(collection.orfs, i)
    orf = @inbounds collection.orfs[i]
    return _extract_sequence(collection.source, orf; kwargs...)
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