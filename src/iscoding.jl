export iscoding

@doc raw"""
    iscoding(
        sequence::LongSequence{DNAAlphabet{4}},
        scoring::NaiveScoringScheme;
        codingmodel::BioMarkovChain,
        noncodingmodel::BioMarkovChain,
        η::Float64 = 1e-5
        )

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
    sequence::NucleicSeqOrView{DNAAlphabet{N}},
    orf::ORF;
    criteria::Function = lordr,
    kwargs...
) where {N}
    return criteria(sequence[orf]; kwargs...)
end


# iscoding.(seq[i for i in allorfs])