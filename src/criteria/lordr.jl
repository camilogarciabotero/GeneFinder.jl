export log_odds_ratio_decision_rule, lordr #lors

@doc raw"""
    log_odds_ratio_decision_rule(
        sequence::LongSequence{DNAAlphabet{4}};
        modela::BioMarkovChain,
        modelb::BioMarkovChain,
        η::Float64 = 1e-5
        )

Check if a given DNA sequence is likely to be coding based on a log-odds ratio.
    The log-odds ratio is a statistical measure used to assess the likelihood of a sequence being coding or non-coding. It compares the probability of the sequence generated by a coding model to the probability of the sequence generated by a non-coding model. If the log-odds ratio exceeds a given threshold (`η`), the sequence is considered likely to be coding.
    It is formally described as a decision rule:

```math
S(X) = \log \left( \frac{{P_C(X_1=i_1, \ldots, X_T=i_T)}}{{P_N(X_1=i_1, \ldots, X_T=i_T)}} \right) \begin{cases} > \eta & \Rightarrow \text{{coding}} \\ < \eta & \Rightarrow \text{{noncoding}} \end{cases}
```

# Arguments
- `sequence::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to be evaluated.

## Keyword Arguments
- `codingmodel::BioMarkovChain`: The transition model for coding regions, (default: `ECOLICDS`).
- `noncodingmodel::BioMarkovChain`: The transition model for non-coding regions, (default: `ECOLINOCDS`)
- `b::Number = 2`: The base of the logarithm used to calculate the log-odds ratio (default: 2).
- `η::Float64 = 1e-5`: The threshold value (eta) for the log-odds ratio (default: 1e-5).

# Returns
- `true` if the sequence is likely to be coding.
- `false` if the sequence is likely to be non-coding.

# Raises
- `ErrorException`: if the length of the sequence is not divisible by 3.
- `ErrorException`: if the sequence contains a premature stop codon.

# Example

```
sequence = dna"ATGGCATCTAG"
iscoding(sequence)  # Returns: true or false
```
"""
function lordr( #log_odds_ratio_decision, also lordr/cudr/kfdr/aadr
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    modela::BioMarkovChain = ECOLICDS,
    modelb::BioMarkovChain = ECOLINOCDS,
    b::Number = 2,
    η::Float64 = 5e-3
) where {N}

    scorea = log_odds_ratio_score(sequence; modela=modela, b=b)
    scoreb = log_odds_ratio_score(sequence; modela=modelb, b=b)
    
    logodds = scorea / scoreb

    if logodds > η
        return true
    else
        false
    end
end

const log_odds_ratio_decision_rule = lordr # criteria