# Structs associated with gene models 
using GenomicFeatures
using BioSequences

abstract type Gene end

abstract type exon end

abstract type intron end


struct Gene
    sequence::DNA
    orfs::Array{DNA, 1}
    start_stops::Array{Tuple{Int64,Int64}, 1}
end

struct Gene
    orfs::Vector{ORF}
    dna_sequence::LongDNA
end

mutable struct ORF <: Gene
    strand::Strand # from GenomicFeatures
    start_position::Int
    stop_position::Int
end

struct ORF <: Gene
    start::Int
    stop::Int
    strand::Char # this allow simple bool elegant ORF.strand == '+' || '-',  note that the GenomicFeatures already got a field Strand!
    frame::Int #  1 to 3 and 
    sequence::BioSequences.LongDNA # Is this type piracy?
end

# The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/
struct Codon <: BioSequence{DNAAlphabet{2}}
    x::UInt8
end

function Codon(iterable)
    length(iterable) == 3 || error("Must have length 3")
    x = zero(UInt)
    for (i, nt) in enumerate(iterable)
        x |= BioSequences.encode(Alphabet(Codon), convert(DNA, nt)) << (6-2i)
    end
    Codon(x % UInt8)
end

Base.length(::Codon) = 3

function BioSequences.extract_encoded_element(x::Codon, i::Int)
    ((x.x >>> (6-2i)) & 3) % UInt
end

Base.copy(seq::Codon) = Codon(seq.x)

BioSequences.has_interface(BioSequence, Codon, [DNA_C, DNA_T, DNA_G], false)