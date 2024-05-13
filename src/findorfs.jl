export findorfs

"""
    findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; ::M, alternative_start=false, min_len=6) where {N, M<:GeneFinderMethod}

This is the main interface method for finding open reading frames (ORFs) in a DNA sequence.

It takes the following arguments:

- `sequence`: The nucleic acid sequence to search for ORFs.
- `method`: The algorithm used to find ORFs. It can be either `NaiveFinder()`, `NaiveFinderScored()` or yet other implementations.

## Keyword Arguments regardless of the finder method:
- `alternative_start`: A boolean indicating whether to consider alternative start codons. Default is `false`.
- `min_len`: The minimum length of an ORF. Default is `6`.

## Keyword Arguments for `NaiveFinderScored()`:
- `scoringscheme::BioMarkovChain`: The scoring scheme to use for the scoring algorithm. Default is `ECOLICDS`.

## Returns
A vector of `ORF` objects representing the found ORFs.

## Example

```julia
sequence = randdnaseq(120)

    120nt DNA Sequence:
    GCCGGACAGCGAAGGCTAATAAATGCCCGTGCCAGTATC…TCTGAGTTACTGTACACCCGAAAGACGTTGTACGCATTT


findorfs(sequence, NaiveFinder())

    1-element Vector{ORF}:
    ORF(77:118, '-', 2, 0.0)
```

"""
function findorfs end

function findorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}},
    ::NaiveFinder;
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}

    # println("Finding ORFs in sequence using $method...")
    return naivefinder(sequence; alternative_start, min_len)::Vector{ORF}
    
end

function findorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}},
    ::NaiveFinderScored;
    alternative_start::Bool = false,
    min_len::Int64 = 6,
    scoringscheme::BioMarkovChain = ECOLICDS
) where {N}

    return naivefinderscored(sequence; alternative_start, min_len, scoringscheme)::Vector{ORF}
    
end