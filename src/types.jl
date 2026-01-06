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