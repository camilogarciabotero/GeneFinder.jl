# Structs associated with gene models 
# using GenomicFeatures
using BioSequences

abstract type Gene end

# abstract type exon end

# abstract type intron end

"""
    struct ORF
        location::UnitRange{Int64}
        strand::Char
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: 

- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence
- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
"""
struct ORF <: Gene
    location::UnitRange{Int64} # Note that it is also called position for gene struct in GenomicAnotations
    strand::Char
end

"""
    struct CDS
        orf::ORF
        sequence::LongDNA
    end

The `CDS` struct represents a coding sequence in a DNA sequence. It has three fields:

- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) displaying the location of the ORF and the associate strand: forward ('+') or reverse ('-')
- `sequence`: a `LongDNA` sequence representing the actual sequence of the CDS
"""
struct CDS <: Gene
    orf::ORF
    sequence::LongDNA
end

"""
    struct Protein
        orf::ORF
        sequence::LongDNA
    end


Similarly to the `CDS` struct, the `Protein` struct represents a encoded protein sequence in a DNA sequence. 
    It has three fields:

- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) of the sequence
- `sequence`: a `LongAA` sequence representing the actual translated sequence of the CDS
"""
struct Protein <: Gene
    orf::ORF
    sequence::LongAA
end

# const START_CODON_MATRIX = PWMSearchQuery([dna"ATG", dna"GTG", dna"TTG"], 1.0)
const STOP_CODONS = [dna"TAG", dna"TAA", dna"TGA"]


# # The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/
# struct Codon <: BioSequence{DNAAlphabet{2}}
#     x::UInt8
# end

# function Codon(iterable)
#     length(iterable) == 3 || error("Must have length 3")
#     x = zero(UInt)
#     for (i, nt) in enumerate(iterable)
#         x |= BioSequences.encode(Alphabet(Codon), convert(DNA, nt)) << (6-2i)
#     end
#     Codon(x % UInt8)
# end

# Base.length(::Codon) = 3

# function BioSequences.extract_encoded_element(x::Codon, i::Int)
#     ((x.x >>> (6-2i)) & 3) % UInt
# end

# Base.copy(seq::Codon) = Codon(seq.x)

# BioSequences.has_interface(BioSequence, Codon, [DNA_C, DNA_T, DNA_G], false)