export hasprematurestop, _varname, _varsymbol, _orfseq #fasta2bioseq

# General purposes methods supporting main functions

"""
    hasprematurestop(seq::LongNucOrView{4})::Bool

Determine whether the `sequence` of type `LongSequence{DNAAlphabet{4}}` contains a premature stop codon.

Returns a boolean indicating whether the `sequence` has more than one stop codon.
"""
function hasprematurestop(seq::NucleicSeqOrView{DNAAlphabet{N}})::Bool where {N}
    
    stopcodons = [LongDNA{4}("TAA"), LongDNA{4}("TAG"), LongDNA{4}("TGA")]  # Create a set of stop codons
    length(seq) % 3 == 0 || error("The sequence is not divisible by 3")
    
    #TODO: this way of checking the stop codon at the end is idiosyncratic and should be improved since other stop codons are not considered
    occursin(biore"T(AG|AA|GA)"dna, @view(seq[end-2:end])) || error("There is no stop codon at the end of the sequence/orf")

    @inbounds for i in 1:3:length(seq) - 4
        codon = seq[i:i+2]
        if codon in stopcodons
            return true
        end
    end

    return false
end

"""
    _orfseq(seq::LongSequence{DNAAlphabet{N}}, start::Int64, stop::Int64, strand::Strand) where {N}

Extract the ORF sequence from the `seq` of type `LongSequence{DNAAlphabet{N}}` given the `start` and `stop` positions and the `strand`.
"""
function _orfseq(seq::LongSequence{DNAAlphabet{N}}, start::Int64, stop::Int64, strand::Strand) where {N}
    if strand == STRAND_POS
        return convert(LongSubSeq, seq[start:stop])
    else
        return convert(LongSubSeq, reverse_complement(seq[start:stop]))
    end
end

"""
    _varname(var::Any)

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
    _varsymbol(var::Any)

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