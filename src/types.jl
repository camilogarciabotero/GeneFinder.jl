# Structs associated with gene models 
# using GenomicFeatures
using BioSequences

"""
    struct ORF
        location::UnitRange{Int64}
        strand::Char
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: location, which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence, and strand, which is a Char indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
"""
struct ORF
    location::UnitRange{Int64} # Note that it is also called position for gene struct in GenomicAnotations
    strand::Char
end

"""
    struct CDS
        location::UnitRange{Int64}
        strand::Char
        sequence::LongDNA
    end

The `CDS` struct represents a coding sequence in a DNA sequence. It has three fields:

    - `location`: a `UnitRange{Int64}` indicating the start and end location of the CDS in the sequence
    - `strand`: a `Char` indicating whether the CDS is on the forward ('+') or reverse ('-') strand of the sequence
    - `sequence`: a `LongDNA` sequence representing the actual sequence of the CDS
    
"""
struct CDS
    location::UnitRange{Int64}
    strand::Char
    sequence::LongDNA
end

"""
    struct Protein
        location::UnitRange{Int64}
        strand::Char
        sequence::LongDNA
    end


Similarly to the `CDS` struct, the `Protein` struct represents a encoded protein sequence in a DNA sequence. 
    It has three fields:

    - `location`: a `UnitRange{Int64}` indicating the start and end locations of the CDS in the sequence
    - `strand`: a `Char` indicating whether the CDS is on the forward ('+') or reverse ('-') strand of the sequence
    - `sequence`: a `LongAA` sequence representing the actual translated sequence of the CDS

"""
struct Protein
    location::UnitRange{Int64}
    strand::Char
    sequence::LongAA
end


const stopcodons = [dna"TAG", dna"TAA", dna"TGA"]

# abstract type Gene end

# abstract type exon end

# abstract type intron end


# struct Gene
#     sequence::DNA
#     orfs::Array{DNA, 1}
#     start_stops::Array{Tuple{Int64,Int64}, 1}
# end

# struct Gene
#     orfs::Vector{ORF}
#     dna_sequence::LongDNA
# end

# mutable struct ORF <: Gene
#     strand::Strand # from GenomicFeatures
#     start_position::Int
#     stop_position::Int
# end

# struct ORF <: Gene
#     start::Int
#     stop::Int
#     strand::Char # this allow simple bool elegant ORF.strand == '+' || '-',  note that the GenomicFeatures already got a field Strand!
#     frame::Int #  1 to 3 and 
#     sequence::BioSequences.LongDNA # Is this type piracy?
# end

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