# General purposes methods supporting main functions

"""
    _hasprematurestop(seq::NucleicSeqOrView{DNAAlphabet{N}})::Bool where {N}

Determine whether the sequence contains a premature stop codon.

A premature stop codon is any stop codon (TAA, TAG, TGA) that appears before
the final codon position in the sequence.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to check.

# Returns
- `Bool`: `true` if a premature stop codon is found, `false` otherwise.

# Throws
- `ArgumentError`: If the sequence length is not divisible by 3.
- `ArgumentError`: If there is no stop codon at the end of the sequence.

See also: [`ORF`](@ref), [`sequence`](@ref)
"""
function _hasprematurestop(seq::NucleicSeqOrView{DNAAlphabet{N}})::Bool where {N}
    stopcodons = (dna"TAA", dna"TAG", dna"TGA")
    
    length(seq) % 3 == 0 || 
        throw(ArgumentError("Sequence length ($(length(seq))) is not divisible by 3"))
    
    @views seq[end-2:end] in stopcodons || 
        throw(ArgumentError("No stop codon at end of sequence"))

    @inbounds for i in 1:3:length(seq) - 3
        @views codon = seq[i:i+2]
        if codon in stopcodons
            return true
        end
    end

    return false
end
