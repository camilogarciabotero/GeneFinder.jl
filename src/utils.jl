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


function _adjust_circular_windows(windows, seqlen)
    adjusted_windows = []
    for w in windows
        w_first = first(w)
        w_last = last(w)
        
        # Normalize positions to 1-based indexing in circular space
        norm_first = mod(w_first - 1, seqlen) + 1
        norm_last = mod(w_last - 1, seqlen) + 1
        
        # Check if window wraps around sequence boundary
        if w_first < 1 || w_last > seqlen || norm_first > norm_last
            # Split into two ranges: end of sequence and beginning of sequence
            if w_first < 1
                # Window extends before sequence start
                range1 = (seqlen + w_first):seqlen
                range2 = 1:w_last
                push!(adjusted_windows, range1, range2)
            else
                # Window extends past sequence end
                range1 = w_first:seqlen
                range2 = 1:norm_last
                push!(adjusted_windows, range1, range2)
            end
        else
            # Window doesn't wrap, use normalized range
            push!(adjusted_windows, norm_first:norm_last)
        end
    end
    return Tuple(adjusted_windows)
end