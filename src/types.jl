# Structs associated with gene models 

abstract type Gene end

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

##### The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/ #####
"""
    Codon <: BioSequence{DNAAlphabet{2}}

A `Struct` representing a codon, which is a subtype of `BioSequence` with
an `Alphabet` of type `DNAAlphabet{2}`. It has a single field `x` of type
`UInt8`. This was implemente in The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/
"""
struct Codon <: BioSequence{DNAAlphabet{2}}
    x::UInt8
end

function Codon(iterable)
    length(iterable) == 3 || error("Must have length 3")
    x = zero(UInt)
    for (i, nt) in enumerate(iterable)
        x |= BioSequences.encode(Alphabet(Codon), convert(DNA, nt)) << (6 - 2i)
    end
    Codon(x % UInt8)
end

Base.length(::Codon) = 3

function BioSequences.extract_encoded_element(x::Codon, i::Int)
    ((x.x >>> (6 - 2i)) & 3) % UInt
end

Base.copy(seq::Codon) = Codon(seq.x)

BioSequences.has_interface(BioSequence, Codon, [DNA_C, DNA_T, DNA_G], false)

Base.count(codon::Codon, sequence::LongDNA) = count(codon, sequence)

function Base.count(codons::Vector{Codon}, sequence::LongDNA)
    a = 0
    @inbounds for i in codons
        a += count(i, sequence)
    end
    return a
end

##### ---------------------------------------- #####

function BioSequences.translate(ntseq::LongSubSeq{DNAAlphabet{4}}; 
    code::GeneticCode=BioSequences.standard_genetic_code, 
    allow_ambiguous_codons=true, 
    alternative_start=false)
ntseq = copy(ntseq)
translate(ntseq; code, allow_ambiguous_codons, alternative_start)
end