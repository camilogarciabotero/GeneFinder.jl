import GenomicFeatures: first, last, length, strand, groupname, metadata
import FASTX: sequence

export Features, RBS, ORF
export features, sequence, source
export groupname, finder, frame, scheme, score, strand, STRAND_BOTH, STRAND_NEG, STRAND_POS, STRAND_NA

#### Ribosome Binding Site (RBS) struct ####
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

# seq[orf.first-bin01.offset.start:orf.first-1]

# motifs = [dna"GGA", dna"GAG", dna"AGG"]

####### End of RBS struct #######


"""
    struct ORF{N,F} <: GenomicFeatures.AbstractGenomicInterval{F}

The `ORF` struct represents an Open Reading Frame (ORF) in genomics.

# Fields
- `groupname::String`: The name of the group to which the ORF belongs.
- `first::Int64`: The starting position of the ORF.
- `last::Int64`: The ending position of the ORF.
- `strand::Strand`: The strand on which the ORF is located.
- `frame::Int`: The reading frame of the ORF.
- `features::Features`: The features associated with the ORF.
- `scheme::Union{Nothing,Function}`: The scheme used for the ORF.

# Main Constructor

```julia
ORF{N,F}(
    groupname::String,
    first::Int64,
    last::Int64,
    strand::Strand,
    frame::Int,
    features::Features,
    scheme::Union{Nothing,Function}
)
```
# Example

A full instance `ORF`

```julia
ORF{4,NaiveFinder}("seq01", 1, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing)
```

A partial instance `ORF`

```julia
ORF{NaiveFinder}(1:33, '+', 1)
```
"""
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
        @assert frame in (1, 2, 3) "Invalid frame value. Frame must be 1, 2, or 3."

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

function ORF{F}(
    range::UnitRange{Int64},
     strand::Char,
     frame::Int
) where {F<:GeneFinderMethod}
    groupname = "unnamedsource"
    first = range.start
    last = range.stop
    features = Features(NamedTuple())
    scheme = nothing
    strand = strand == '+' ? STRAND_POS : STRAND_NEG
    return ORF{4,F}(groupname, first, last, strand, frame, features, scheme)
end

"""
    sequence(i::ORF{N,F})

Extracts the DNA sequence corresponding to the given open reading frame (ORF).

# Arguments
- `i::ORF{N,F}`: The open reading frame (ORF) for which the DNA sequence needs to be extracted.

# Returns
- The DNA sequence corresponding to the given open reading frame (ORF).

"""
function sequence(i::ORF{N,F}) where {N,F}
    seqsymb = Symbol(i.groupname)
    seq = getfield(Main, seqsymb)
    return i.strand == STRAND_POS ? @view(seq[i.first:i.last]) : reverse_complement(@view(seq[i.first:i.last])) #seq[i.first:i.last]
end


"""
    source(i::ORF{N,F})

Get the source sequence associated with the given `ORF` object.

# Arguments
- `i::ORF{N,F}`: The `ORF` object for which to retrieve the source sequence.

# Returns
The source sequence associated with the `ORF` object.

# Examples
```julia
seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
orfs = findorfs(seq)
source(orfs[1])

44nt DNA Sequence:
ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA
```

!!! warning
    The `source` method works if the sequence is defined in the global scope. Otherwise it will throw an error. 
    For instance a common failure is to define a simple `ORF` that by defualt will have an "unnamedsource" as `groupname` 
    and then try to get the source sequence. 

    ```julia
    orf = ORF{NaiveFinder}(1:33, '+', 1)
    source(orf)

    ERROR: UndefVarError: `unnamedsource` not defined
    Stacktrace:
     [1] source(i::ORF{4, NaiveFinder})
       @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:192
     [2] top-level scope
       @ REPL[12]:1
    ```

"""
function source(i::ORF{N,F}) where {N,F}
    seqsymb = Symbol(i.groupname)
    return getfield(Main, seqsymb)
end


"""
    features(i::ORF{N,F})

Extracts the features from an `ORF` object.

# Arguments
- `i::ORF{N,F}`: An `ORF` object.

# Returns
The features of the `ORF` object. Those could be defined by each `GeneFinderMethod`.

"""
function features(i::ORF{N,F}) where {N,F}
    return i.features.fts
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

function score(i::ORF{N,F}) where {N,F}
    return i.features.fts[:score]
end

function length(i::ORF{N,F}) where {N,F}
    return length(sequence(i))
end

function strand(i::ORF{N,F}) where {N,F}
    return i.strand
end

function metadata(i::ORF{N,F}) where {N,F}
    return features(i)
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