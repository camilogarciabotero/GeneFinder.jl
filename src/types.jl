import GenomicFeatures: first, last, length, strand, groupname, metadata

export ORFI, OpenReadingFrameInterval
export features, sequence, source
export groupname, finder, frame, strand, STRAND_BOTH, STRAND_NEG, STRAND_POS, STRAND_NA

"""
    struct ORFI{N,F} <: AbstractGenomicInterval{F}

The `ORFI` struct represents an Open Reading Frame Interval (ORFI) in genomics.

# Fields
- `groupname::String`: The name of the group to which the ORFI belongs.
- `first::Int64`: The starting position of the ORFI.
- `last::Int64`: The ending position of the ORFI.
- `strand::Strand`: The strand on which the ORFI is located.
- `frame::Int`: The reading frame of the ORFI.
- `seq::LongSubSeq{DNAAlphabet{N}}`: The DNA sequence of the ORFI.
- `features::Features`: The features associated with the ORFI.

# Main Constructor

```julia
ORFI{N,F}(
    groupname::String,
    first::Int64,
    last::Int64,
    strand::Strand,
    frame::Int,
    features::Features,
    seq::LongSubSeq{DNAAlphabet{N}}
)
```
# Example

A full instance `ORFI`

```julia
ORFI{4,NaiveFinder}("seq01", 1, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing)
```

A partial instance `ORFI`

```julia
ORFI{NaiveFinder}(1:33, '+', 1)
```
"""
struct OpenReadingFrameInterval{N,F} <: AbstractGenomicInterval{F} #GenomicFeatures
    groupname::String
    first::Int64
    last::Int64
    strand::Strand
    frame::Int8
    seq::LongSubSeq{DNAAlphabet{N}}
    features::NamedTuple
end

function OpenReadingFrameInterval{N,F}(
    ::Type{F}, #finder
    groupname::String,
    first::Int64,
    last::Int64,
    strand::Strand,
    frame::Int8,
    seq::LongSubSeq{DNAAlphabet{N}},
    features::NamedTuple # ::Dict{Symbol,Any} or # ::@NamedTuple{score::Float64, rbs::Any} or @NamedTuple{Vararg{typeof(...)}} NTuple?
) where {N,F<:GeneFinderMethod}
    # @assert frame in (1, 2, 3) "Invalid frame value. Frame must be 1, 2, or 3."
    # frame in (1, 2, 3) || throw(ArgumentError("The source sequence of the ORFI and the given sequence are different"))
    return ORFI{N,F}(groupname, first, last, strand, frame, seq, features) #finder seq schemes
end

const ORFI = OpenReadingFrameInterval

"""
    sequence(i::ORFI{N,F})

Extracts the DNA sequence corresponding to the given open reading frame (ORFI).

# Arguments
- `i::ORFI{N,F}`: The open reading frame (ORFI) for which the DNA sequence needs to be extracted.

# Returns
- The DNA sequence corresponding to the given open reading frame (ORFI).

"""
function sequence(i::ORFI{N,F}) where {N,F}
    return i.seq
end

"""
    source(i::ORFI{N,F})

Get the source sequence associated with the given `ORFI` object.

# Arguments
- `i::ORFI{N,F}`: The `ORFI` object for which to retrieve the source sequence.

# Returns
The source sequence associated with the `ORFI` object.

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
    For instance a common failure is to define a simple `ORFI` that by defualt will have an "unnamedsource" as `groupname` 
    and then try to get the source sequence. 

    ```julia
    orf = ORFI{NaiveFinder}(1:33, '+', 1)
    source(orf)

    ERROR: UndefVarError: `unnamedsource` not defined
    Stacktrace:
     [1] source(i::ORFI{4, NaiveFinder})
       @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:192
     [2] top-level scope
       @ REPL[12]:1
    ```

"""
function source(i::ORFI{N,F}) where {N,F}
    seqsymb = Symbol(i.groupname)
    try
        return getfield(Main, seqsymb)
    catch e
        error("The source sequence of the ORFI is defined as $(i.groupname). Make sure to either define it or supply the correct source sequence.")
    end
end

"""
    features(i::ORFI{N,F})

Extracts the features from an `ORFI` object.

# Arguments
- `i::ORFI{N,F}`: An `ORFI` object.

# Returns
The features of the `ORFI` object. Those could be defined by each `GeneFinderMethod`.

"""
function features(i::ORFI{N,F}) where {N,F}
    return i.features
end

function groupname(i::ORFI{N,F}) where {N,F}
    return i.groupname
end

function first(i::ORFI{N,F}) where {N,F}
    return i.first
end

function last(i::ORFI{N,F}) where {N,F}
    return i.last
end

function frame(i::ORFI{N,F}) where {N,F}
    return i.frame
end

function strand(i::ORFI{N,F}) where {N,F}
    return i.strand
end

function metadata(i::ORFI{N,F}) where {N,F}
    return features(i)
end

finder(i::ORFI{N,F}) where {N,F} = F

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