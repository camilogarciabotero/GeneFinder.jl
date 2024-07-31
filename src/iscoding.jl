export iscoding

@doc raw"""
    iscoding(sequence::LongSequence{DNAAlphabet{4}}; criteria::Function = lordr, kwargs...) -> Bool

Check if a given DNA sequence is likely to be coding based on a scoring scheme.

```
sequence = dna"ATGGCATCTAG"
iscoding(sequence)  # Returns: true or false
```
"""
function iscoding end

function iscoding(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    criteria::Function = lordr,
    kwargs...
) where {N}
    return criteria(sequence; kwargs...)
end

function iscoding(
    orf::ORFI{F};
    criteria::Function = lordr,
    kwargs...
) where {F<:GeneFinderMethod}
    return criteria(sequence(orf); kwargs...)
end

# iscoding.(seq[i for i in allorfs])