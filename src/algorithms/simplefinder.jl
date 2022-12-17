using BioSequences
using TestItems

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

"""
    `simplefinder(sequence::LongDNA)`

The simplest algorithm that finds ORFs in a DNA sequence.

The simplefinder function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
    This function has not ORFs size and overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.

# Examples
```jldoctest
julia> using BioSequences

julia> seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA";

julia> simplefinder(seq)
5-element Vector{ORF}:
 ORF(1:33, '+')
 ORF(4:33, '+')
 ORF(8:22, '+')
 ORF(12:29, '+')
 ORF(16:33, '+')

```
"""
function simplefinder(sequence::LongDNA)
    orfs = Vector{ORF}()
    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence
        
        start_codon = ExactSearchQuery(dna"ATG", iscompatible)
        start_codon_indices = findall(start_codon, seq)

        for i in start_codon_indices
            orf = nothing
            j = i.start
            while j < length(seq) - 3
                if seq[j:j+2] ∈ [dna"TAG", dna"TAA", dna"TGA"]
                    push!(orfs, orf)
                    break
                end
                orf = ORF(i.start:j+5, strand)
                j += 3 
            end
        end
    end
    return orfs
end

@testitem "simplefinder test" default_imports=true begin
    using BioSequences
    
    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

    @test simplefinder(seq) == [ORF(1:33, '+'), ORF(4:33, '+'), ORF(8:22, '+'), ORF(12:29, '+'), ORF(16:33, '+')]
end

"""
    `findcds(sequence::LongDNA)`

A function to generete CDSs sequence out of a DNA sequence.

The `findcds` function takes a `LongDNA` sequence and returns a `Vector{CDS}` 
    containing the coding sequences (CDSs) found in the sequence. 
    It uses the `simplefinder` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.

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
            seq = reversedsequence[i.location]
        else
            seq = sequence[i.location]
        end
        cds = CDS(i.location, i.strand, seq)
        push!(seqs, cds)
    end
    return seqs
end


"""
    `findproteins(sequence::LongDNA)`

As its name suggest this function generate the possible proteins directly from a DNA sequence. 
    The `findcds` function takes a `LongDNA` sequence and returns a `Vector{CDS}` containing the 
    coding sequences (CDSs) found in the sequence. 
"""
function findproteins(sequence::LongDNA)
    cds = findcds(sequence)
    proteins = Vector{Protein}()
    for i in cds
        proteinseq = translate(i.sequence)
        protein = Protein(i.location, i.strand, proteinseq)
        push!(proteins, protein)
    end
    return proteins
end