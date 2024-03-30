export findorfs

"""
    findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; findermethod=naivefinder, alternative_start=false, min_len=6) where {N}

This is the main interface method for finding open reading frames (ORFs) in a DNA sequence.

    It takes the following arguments:
- `sequence`: The nucleic acid sequence to search for ORFs.
- `findermethod`: The algorithm used to find ORFs. Default is `naivefinder`.
- `alternative_start`: A boolean indicating whether to consider alternative start codons. Default is `false`.
- `min_len`: The minimum length of an ORF. Default is `6`.

Returns a vector of `ORF` objects representing the found ORFs.
"""
function findorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    findermethod::Function = naivefinder,
    alternative_start::Bool = false,
    min_len::Int64 = 6
    ) where {N}

    return findermethod(sequence; alternative_start, min_len)::Vector{ORF}
    
end