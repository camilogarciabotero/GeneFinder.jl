export hasprematurestop, fasta2bioseq, _varname, _varsymbol

# General purposes methods supporting main functions

"""
    hasprematurestop(sequence::LongNucOrView{4})::Bool

Determine whether the `sequence` of type `LongSequence{DNAAlphabet{4}}` contains a premature stop codon.

Returns a boolean indicating whether the `sequence` has more than one stop codon.
"""
function hasprematurestop(sequence::NucleicSeqOrView{DNAAlphabet{N}})::Bool where {N}
    
    stopcodons = [LongDNA{4}("TAA"), LongDNA{4}("TAG"), LongDNA{4}("TGA")]  # Create a set of stop codons
    
    length(sequence) % 3 == 0 || error("The sequence is not divisible by 3")
    
    #TODO: this way of checking the stop codon at the end is idiosyncratic and should be improved since other stop codons are not considered
    #TODO: This is currently used by the `getindex` method in `extended.jl`
    occursin(biore"T(AG|AA|GA)"dna, sequence[end-2:end]) || error("There is no stop codon at the end of the sequence/orf")

    @inbounds for i in 1:3:length(sequence) - 4
        codon = sequence[i:i+2]
        if codon in stopcodons
            return true
        end
    end

    return false
end

"""
    fasta2bioseq(input::String)

Converts a FASTA formatted file (even if it is a multi-fasta) to an array of `LongSequence{DNAAlphabet{4}}` objects.
"""
function fasta2bioseq(input::AbstractString)::Vector{LongSequence{DNAAlphabet{4}}}
    FASTAReader(open(input)) do reader
        return [LongSequence{DNAAlphabet{4}}(sequence(record)) for record in reader]
    end
end

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