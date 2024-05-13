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
    sequence::NucleicSeqOrView{DNAAlphabet{N}},
    ::NaiveScoringScheme;
    codingmodel::BioMarkovChain = ECOLICDS,
    noncodingmodel::BioMarkovChain = ECOLINOCDS,
    η::Float64 = 1e-5
) where {N}
    return isnaivecoding(sequence; codingmodel, noncodingmodel, η)
end


# function iscoding(
#     orf::ORF,
#     sequence::NucleicSeqOrView{DNAAlphabet{N}};
#     scoring::NaiveScoringScheme = NaiveScoringScheme(),
#     codingmodel::BioMarkovChain = ECOLICDS,
#     noncodingmodel::BioMarkovChain = ECOLINOCDS,
#     η::Float64 = 1e-5
# ) where {N}
#     return isnaivecoding(sequence[orf]; codingmodel, noncodingmodel, η)
# end


# iscoding.(seq[i for i in allorfs])