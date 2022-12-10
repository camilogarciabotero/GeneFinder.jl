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

