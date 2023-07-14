# Structs associated with gene models 

abstract type Gene end

const LongNucOrView{N} = Union{LongSequence{<:NucleicAcidAlphabet{N}},LongSubSeq{<:NucleicAcidAlphabet{N}}}

"""
    struct GeneFeatures
        seqname::String
        start::Int64
        stop::Int64
        score::Float64
        strand::Char
        frame::'.'
        attribute::
    end

This is the main Gene struct, based on the fields that could be found in a GFF3, still needs to be defined correctly,
    The idea is correct the frame and attributes that will have something like a possible list (id=Char;name=;locus_tag).
    The `write` and `get` functions should have a dedicated method for this struct.
"""
struct GeneFeatures
    seqname::String
    start::Int64
    stop::Int64
    score::Float64
    strand::Char
    frame::String
    attribute::String # Should be a Dict perhaps
end

"""
    struct ORF
        location::UnitRange{Int64}
        strand::Char
        frame::Integer
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: 

- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence
- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
"""
struct ORF <: Gene
    location::UnitRange{Int64} # Note that it is also called position for gene struct in GenomicAnotations
    strand::Char
    frame::Integer # 1, 2, or 3
end


##### The following implementation is from https://biojulia.dev/BioSequences.jl/stable/interfaces/ #####
# """
#     Codon <: BioSequence{DNAAlphabet{2}}

# A `Struct` representing a codon, which is a subtype of `BioSequence` with
# an `Alphabet` of type `DNAAlphabet{2}`. It has a single field `x` of type
# `UInt8`. This was implemente in The following implementation is from https://biojulia.dev/BioSequences.jl/stable/interfaces/
# """
# struct Codon <: BioSequence{DNAAlphabet{2}}
#     x::UInt8

#     function Codon(iterable)
#         length(iterable) == 3 || error("Must have length 3")
#         x = zero(UInt)
#         for (i, nt) in enumerate(iterable)
#             x |= BioSequences.encode(Alphabet(Codon), convert(DNA, nt)) << (6 - 2i)
#         end
#         new(x % UInt8)
#     end 
# end

# Base.copy(seq::Codon) = Codon(seq.x)
# Base.count(codon::Codon, sequence::LongSequence{DNAAlphabet{4}}) = count(codon, sequence)
# Base.length(::Codon) = 3

# function Base.count(codons::Vector{Codon}, sequence::LongSequence{DNAAlphabet{4}})
#     a = 0
#     @inbounds for i in codons
#         a += count(i, sequence)
#     end
#     return a
# end

# # has_interface(BioSequence, Codon, [DNA_C, DNA_T, DNA_G], false)

# function extract_encoded_element(x::Codon, i::Int)
#     ((x.x >>> (6 - 2i)) & 3) % UInt
# end

# function translate(
#     ntseq::Codon;
#     code::GeneticCode = standard_genetic_code,
#     allow_ambiguous_codons = true,
#     alternative_start = false,
# )   
#     translate(ntseq; code, allow_ambiguous_codons, alternative_start)
# end

##### ---------------------------------------- #####
