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

"""
    _orfseq(seq::LongSequence{DNAAlphabet{N}}, start::Int64, stop::Int64, strand::Strand) where {N} -> LongSubSeq{DNAAlphabet{N}}

Extract the ORFI sequence from the `seq` of type `LongSequence{DNAAlphabet{N}}` given the `start` and `stop` positions and the `strand`.
"""
function _orfseq(seq::SeqOrView{DNAAlphabet{N}}, start::Int64, stop::Int64, strand::Strand) where {N}
    if strand == STRAND_POS
        return convert(LongSubSeq, seq[start:stop])
    else
        return convert(LongSubSeq, reverse_complement(seq[start:stop]))
    end
end

"""
    _varname(var::Any) -> String

Get the name of the variable `var` in the current scope.
"""
function _varname(var::Any)
    for name in names(Main)
        try
            if getfield(Main, name) === var
                return string(name)
            end
        catch e
            # Skip if getfield fails, e.g., for names that cannot be accessed directly
            continue
        end
    end
    return nothing
end

"""
    _varsymbol(var::Any) -> Symbol

Get the symbol of the variable `var` in the current scope.
"""
function _varsymbol(var::Any)
    for name in names(Main)
        try
            if getfield(Main, name) === var
                return Symbol(name)
            end
        catch e
            # Skip if getfield fails, e.g., for names that cannot be accessed directly
            continue
        end
    end
    return nothing
end

# export location
# function location(orf::ORFI{N,F})::UnitRange{Int64} where {N,F}
#     return orf.first:orf.last
# end

# export _norfs
# function _norfs(seq::NucleicSeqOrView{DNAAlphabet{N}}) where {N}
#     cpseq = copy(seq)
#     orfs = findorfs(cpseq; finder=NaiveCollector)

#     for l in location.(orfs)
#         deleteat!(cpseq, l)
#     end

#     return cpseq
# end