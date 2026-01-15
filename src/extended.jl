# Methods from main packages that expand their fuctions to this package structs

@inline Base.isless(a::ORF{F}, b::ORF{F}) where {F} = isless(a.range, b.range)
@inline Base.length(orf::ORF{F}) where {F} = length(orf.range)
@inline Base.getindex(seq::SeqOrView{DNAAlphabet{N}}, orf::ORF{F}) where {N,F} = _orfseq(seq, leftposition(orf), rightposition(orf), strand(orf))
@inline Base.range(orf::ORF{F}) where {F} = orf.range

function Base.:(==)(a::ORF{F}, b::ORF{F}) where {F}
    return a.range == b.range && a.strand == b.strand && a.frame == b.frame && a.seqid == b.seqid
end

# Methods extending `Strand`

function Base.show(io::IO, s::Strand)
    if s === PSTRAND
        print(io, '+')
    elseif s === NSTRAND
        print(io, '-')
    else
        print(io, "unknown")
    end
end

function Base.print(io::IO, s::Strand)
    if s === PSTRAND
        print(io, '+')
    elseif s === NSTRAND
        print(io, '-')
    else
        print(io, "unknown")
    end
end

## Methods from BioMarkovChains that expand their fuctions to this package structs
import BioMarkovChains: log_odds_ratio_score
@inline log_odds_ratio_score(orf::ORF{F}, models...; kwargs...) where {F} = log_odds_ratio_score(sequence(orf), models...; kwargs...)

@inline log_odds_ratio_score(collection::ORFCollection, i::Int, models...; kwargs...) = log_odds_ratio_score(sequence(collection, i), models...; kwargs...)
@inline log_odds_ratio_score(collection::ORFCollection, orf::ORF, models...; kwargs...) = log_odds_ratio_score(sequence(collection, orf), models...; kwargs...)
@inline log_odds_ratio_score(seq::SeqOrView{DNAAlphabet{N}}, orf::ORF, models...; kwargs...) where {N} = log_odds_ratio_score(seq[orf], models...; kwargs...)

# ════════════════════════════════════════════════════════════════════════════════
# Extended Interface for ORF and ORFCollection
# ════════════════════════════════════════════════════════════════════════════════
#
# This file contains extended methods for ORF and ORFCollection types:
# - Translation methods
# - Comparison and equality methods
# - Set operations (union, intersect, setdiff, issubset)
# - Hash methods for use in Sets and Dicts
#
# ════════════════════════════════════════════════════════════════════════════════

# export translations

# ────────────────────────────────────────────────────────────────────────────────
# Translation Methods
# ────────────────────────────────────────────────────────────────────────────────
## Methods from BioSequences that expand their fuctions to this package structs
import BioSequences: translate
"""
    translate(collection::ORFCollection, i::Int; kwargs...)

Translate the DNA sequence of the ORF at index `i` to an amino acid sequence.

# Arguments
- `collection::ORFCollection`: The collection containing the ORF and source sequence.
- `i::Int`: The index of the ORF to translate.

# Keywords
- All keyword arguments are passed to `BioSequences.translate`.

# Returns
- `LongAA`: The translated amino acid sequence.

# Example
```julia
collection = findorfs(seq)
protein = translate(collection, 1)  # Translate first ORF
```

See also: [`translations`](@ref), [`sequence`](@ref)
"""
function translate(collection::ORFCollection, i::Int; kwargs...)
    return translate(sequence(collection, i); kwargs...)
end

"""
    translate(collection::ORFCollection, orf::ORF; kwargs...)

Translate the DNA sequence of a specific ORF to an amino acid sequence.

# Arguments
- `collection::ORFCollection`: The collection containing the source sequence.
- `orf::ORF`: The ORF to translate.

# Keywords
- All keyword arguments are passed to `BioSequences.translate`.

# Returns
- `LongAA`: The translated amino acid sequence.

# Example
```julia
collection = findorfs(seq)
orf = collection[1]
protein = translate(collection, orf)
```
"""
function translate(collection::ORFCollection, orf::ORF; kwargs...)
    return translate(sequence(collection, orf); kwargs...)
end


# ────────────────────────────────────────────────────────────────────────────────
# ORF Comparison Methods
# ────────────────────────────────────────────────────────────────────────────────

"""
    Base.:(==)(a::ORF, b::ORF) -> Bool

Test equality of two ORFs based on range, strand, and frame.

Two ORFs are considered equal if they have the same range, strand, and frame,
regardless of their finder method type or features.

# Example
```julia
orf1 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
orf2 = ORF{NaiveFinderLazy}(1:33, PSTRAND, Int8(1))
orf1 == orf2  # true (same coordinates)

orf3 = ORF{NaiveFinder}(1:33, NSTRAND, Int8(1))
orf1 == orf3  # false (different strand)
```
"""
function Base.:(==)(a::ORF, b::ORF)
    return a.range == b.range && a.strand == b.strand && a.frame == b.frame
end

"""
    Base.isless(a::ORF, b::ORF) -> Bool

Compare two ORFs for sorting purposes.

ORFs are compared by:
1. Left position (start of range)
2. Right position (end of range)
3. Strand (PSTRAND < NSTRAND)
4. Frame

This ordering ensures consistent sorting of ORFs by genomic position.

# Example
```julia
orf1 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
orf2 = ORF{NaiveFinder}(10:42, PSTRAND, Int8(1))
orf1 < orf2  # true (earlier start position)
```
"""
function Base.isless(a::ORF, b::ORF)
    leftposition(a) != leftposition(b) && return leftposition(a) < leftposition(b)
    rightposition(a) != rightposition(b) && return rightposition(a) < rightposition(b)
    a.strand != b.strand && return Int(a.strand) < Int(b.strand)
    return a.frame < b.frame
end

"""
    Base.hash(orf::ORF, h::UInt) -> UInt

Compute hash for an ORF based on range, strand, and frame.

Enables use of ORFs in Sets and as Dict keys.
"""
function Base.hash(orf::ORF, h::UInt)
    h = hash(orf.range, h)
    h = hash(orf.strand, h)
    h = hash(orf.frame, h)
    return h
end

# ────────────────────────────────────────────────────────────────────────────────
# ORFCollection Comparison Methods
# ────────────────────────────────────────────────────────────────────────────────

"""
    Base.:(==)(a::ORFCollection, b::ORFCollection) -> Bool

Test equality of two ORFCollections.

Two collections are equal if they have the same ORFs (by range, strand, frame)
in the same order. Source sequences are not compared.

# Example
```julia
seq = dna"ATGATGCATGCATGCATGCTAG"
c1 = findorfs(seq, finder=NaiveFinder)
c2 = findorfs(seq, finder=NaiveFinderLazy)

# Equal if they found the same ORFs
c1 == c2
```
"""
function Base.:(==)(a::ORFCollection, b::ORFCollection)
    length(a) != length(b) && return false
    for (orf_a, orf_b) in zip(a.orfs, b.orfs)
        orf_a == orf_b || return false
    end
    return true
end

"""
    Base.isequal(a::ORFCollection, b::ORFCollection) -> Bool

Strict equality test including source sequence content.

Unlike `==`, this also verifies that the source sequences are identical.
"""
function Base.isequal(a::ORFCollection, b::ORFCollection)
    a == b || return false
    return a.source == b.source
end

"""
    Base.hash(c::ORFCollection, h::UInt) -> UInt

Compute hash for an ORFCollection based on its ORFs.
"""
function Base.hash(c::ORFCollection, h::UInt)
    for orf in c.orfs
        h = hash(orf, h)
    end
    return h
end

function Base.getindex(c::ORFCollection, orf::ORF)
    for (i, o) in enumerate(c.orfs)
        if o == orf
            return i
        end
    end
    throw(KeyError("ORF not found in collection"))
end

function Base.length(c::ORFCollection, idx::Int)
    return length(c.orfs[idx])
end

# ────────────────────────────────────────────────────────────────────────────────
# Set Operations for ORFCollection
# ────────────────────────────────────────────────────────────────────────────────

"""
    Base.issubset(a::ORFCollection, b::ORFCollection) -> Bool

Test if all ORFs in collection `a` are also in collection `b`.

# Example
```julia
c1 = findorfs(seq, finder=NaiveFinder, minlen=30)
c2 = findorfs(seq, finder=NaiveFinder, minlen=6)

issubset(c1, c2)  # true (longer ORFs are subset of all ORFs)
```
"""
function Base.issubset(a::ORFCollection, b::ORFCollection)
    b_set = Set(b.orfs)
    return all(orf -> orf in b_set, a.orfs)
end

"""
    Base.intersect(a::ORFCollection, b::ORFCollection) -> Vector{ORF}

Find ORFs that appear in both collections (by range, strand, frame).

# Example
```julia
c1 = findorfs(seq1)
c2 = findorfs(seq2)
common_orfs = intersect(c1, c2)
```
"""
function Base.intersect(a::ORFCollection, b::ORFCollection)
    b_set = Set(b.orfs)
    return filter(orf -> orf in b_set, a.orfs)
end

"""
    Base.setdiff(a::ORFCollection, b::ORFCollection) -> Vector{ORF}

Find ORFs in collection `a` that are not in collection `b`.

# Example
```julia
c1 = findorfs(seq, finder=NaiveFinder, minlen=6)
c2 = findorfs(seq, finder=NaiveFinder, minlen=30)
short_orfs = setdiff(c1, c2)  # ORFs shorter than 30
```
"""
function Base.setdiff(a::ORFCollection, b::ORFCollection)
    b_set = Set(b.orfs)
    return filter(orf -> orf ∉ b_set, a.orfs)
end

"""
    Base.union(a::ORFCollection, b::ORFCollection) -> Vector{ORF}

Combine ORFs from both collections, removing duplicates.

# Example
```julia
c1 = findorfs(seq1)
c2 = findorfs(seq2)
all_orfs = union(c1, c2)
```
"""
function Base.union(a::ORFCollection, b::ORFCollection)
    return collect(union(Set(a.orfs), Set(b.orfs)))
end
