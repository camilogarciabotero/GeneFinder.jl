# this will be the main functions taking all the algorithms
using BioSequences

include("algorithms/simplefinder.jl")

"""
    struct CDS
        position::UnitRange{Int64}
        strand::Char
        sequence::LongDNA
    end

The `CDS` struct represents a coding sequence in a DNA sequence. It has three fields:

    - `position`: a `UnitRange{Int64}` indicating the start and end positions of the CDS in the sequence
    - `strand`: a `Char` indicating whether the CDS is on the forward ('+') or reverse ('-') strand of the sequence
    - `sequence`: a `LongDNA` sequence representing the actual sequence of the CDS
    
"""
struct CDS
    position::UnitRange{Int64}
    strand::Char
    sequence::LongDNA
end

"""
    `findcds(sequence::LongDNA)`

A function to generete CDSs sequence out of a DNA sequence.

The `findcds` function takes a `LongDNA` sequence and returns a `Vector{CDS}` containing the coding sequences (CDSs) found in the sequence. It uses the `simplefinder` function to find open reading frames (ORFs) in the sequence, and then it extracts the actual CDS sequence from each ORF. The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.

# Examples
```jldoctest
julia> seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA";

findcds(seq)
5-element Vector{CDS}:
 CDS(1:33, '+', ATGATGCATGCATGCATGCTAGTAACTAGCTAG)
 CDS(4:33, '+', ATGCATGCATGCATGCTAGTAACTAGCTAG)
 CDS(8:22, '+', ATGCATGCATGCTAG)
 CDS(12:29, '+', ATGCATGCTAGTAACTAG)
 CDS(16:33, '+', ATGCTAGTAACTAGCTAG)
```
"""
function findcds(sequence::LongDNA)

    orfs = simplefinder(sequence)
    seqs = Vector{CDS}()
    
    for i in orfs
        if i.strand == '-'
            reversedsequence = reverse_complement(sequence)
            seq = reversedsequence[i.position]
        else
            seq = sequence[i.position]
        end
        cds = CDS(i.position, i.strand, seq)
        push!(seqs, cds)
    end
    return seqs
end


"""
    struct Protein
        position::UnitRange{Int64}
        strand::Char
        sequence::LongDNA
    end


Similarly to the `CDS` struct, the `Protein` struct represents a encoded protein sequence in a DNA sequence. It has three fields:

    - `position`: a `UnitRange{Int64}` indicating the start and end positions of the CDS in the sequence
    - `strand`: a `Char` indicating whether the CDS is on the forward ('+') or reverse ('-') strand of the sequence
    - `sequence`: a `LongAA` sequence representing the actual translated sequence of the CDS

"""
struct Protein
    position::UnitRange{Int64}
    strand::Char
    sequence::LongAA
end


"""
    `findproteins(sequence::LongDNA)`

As its name suggest this function generate the possible proteins directly from a DNA sequence. The 
"""
function findproteins(sequence::LongDNA)
    cds = findcds(sequence)
    proteins = Vector{Protein}()
    for i in cds
        proteinseq = translate(i.sequence)
        protein = Protein(i.position, i.strand, proteinseq)
        push!(proteins, protein)
    end
    return proteins
end




#function findgenes(
#    sequence::LongDNA,
#    algorithm::PredictionAlgorithms,
#    type::GeneticCode)