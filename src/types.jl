# Structs associated with gene models 

abstract type Gene end

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
    attribute::String
end

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
        sequence::LongSubSeq{DNAAlphabet{4}}
    end

The `CDS` struct represents a coding sequence in a DNA sequence. It has three fields:

- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) displaying the location of the ORF and the associate strand: forward ('+') or reverse ('-')
- `sequence`: a `LongDNA` sequence representing the actual sequence of the CDS
"""
struct CDS
    sequence::LongSubSeq{DNAAlphabet{4}} #LongDNA
    orf::ORF
end

"""
    struct Protein
        sequence::LongSequence
        orf::ORF
    end
    
Similarly to the `CDS` struct, the `Protein` struct represents a encoded protein sequence in a DNA sequence. 
    It has three fields:

- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) of the sequence
- `sequence`: a `LongSequence` sequence representing the actual translated sequence of the CDS
"""
struct Protein
    sequence::LongSubSeq{AminoAcidAlphabet}
    orf::ORF
end


"""
    DTCM(alphabet::Vector{DNA})

A data structure for storing a DNA Transition Count Matrix (DTCM). The DTCM is a square matrix where each row and column corresponds to a nucleotide in the given `alphabet`. The value at position (i, j) in the matrix represents the number of times that nucleotide i is immediately followed by nucleotide j in a DNA sequence. 

Fields:
- `order::Dict{DNA, Int64}`: A dictionary that maps each nucleotide in the `alphabet` to its corresponding index in the matrix.
- `counts::Matrix{Int64}`: The actual matrix of counts.

Internal function:
- `DTCM(alphabet::Vector{DNA})`: Constructs a new `DTCM` object with the given `alphabet`. This function initializes the `order` field by creating a dictionary that maps each nucleotide in the `alphabet` to its corresponding index in the matrix. It also initializes the `counts` field to a matrix of zeros with dimensions `len x len`, where `len` is the length of the `alphabet`.

Example usage:
```julia
alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]
dtcm = DTCM(alphabet)
```
"""
struct DTCM
    order::Dict{DNA,Int64}
    counts::Matrix{Int64}

    function DTCM(alphabet::Vector{DNA})

        len = length(alphabet)

        order = Dict{DNA,Int}()
        for (i, nucleotide) in enumerate(sort(alphabet))
            order[nucleotide] = i
        end
        counts = zeros(Int64, len, len)
        new(order, counts)
    end
end


"""
    DTPM(alphabet::Vector{DNA})

A data structure for storing a DNA Transition Probability Matrix (DTPM). The DTPM is a square matrix where each row and column corresponds to a nucleotide in the given `alphabet`. The value at position (i, j) in the matrix represents the probability of transitioning from nucleotide i to nucleotide j in a DNA sequence. 

Fields:
- `order::Dict{DNA, Int64}`: A dictionary that maps each nucleotide in the `alphabet` to its corresponding index in the matrix.
- `probabilities::Matrix{Float64}`: The actual matrix of probabilities.

Example usage:
```julia
alphabet = [DNA_A, DNA_C, DNA_G, DNA_T]
dtpm = DTPM(alphabet)
```
"""
struct DTPM
    order::Dict{DNA,Int64}
    probabilities::Matrix{Float64}
end

##### The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/ #####
# """
#     Codon <: BioSequence{DNAAlphabet{2}}

# A `Struct` representing a codon, which is a subtype of `BioSequence` with
# an `Alphabet` of type `DNAAlphabet{2}`. It has a single field `x` of type
# `UInt8`. This was implemente in The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/
# """
# struct Codon <: BioSequence{DNAAlphabet{2}}
#     x::UInt8
# end

# function Codon(iterable)
#     length(iterable) == 3 || error("Must have length 3")
#     x = zero(UInt)
#     for (i, nt) in enumerate(iterable)
#         x |= BioSequences.encode(Alphabet(Codon), convert(DNA, nt)) << (6 - 2i)
#     end
#     Codon(x % UInt8)
# end

# Base.length(::Codon) = 3

# function BioSequences.extract_encoded_element(x::Codon, i::Int)
#     ((x.x >>> (6 - 2i)) & 3) % UInt
# end

# Base.copy(seq::Codon) = Codon(seq.x)

# BioSequences.has_interface(BioSequence, Codon, [DNA_C, DNA_T, DNA_G], false)

# Base.count(codon::Codon, sequence::LongDNA) = count(codon, sequence)

# function Base.count(codons::Vector{Codon}, sequence::LongDNA)
#     a = 0
#     @inbounds for i in codons
#         a += count(i, sequence)
#     end
#     return a
# end

##### ---------------------------------------- #####
