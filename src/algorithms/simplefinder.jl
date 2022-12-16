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
    `simplefinder(sequence::LongDNA)`

The simplest algorithm that finds ORFs in a DNA sequence.

The simplefinder function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
    This function has not ORFs size and overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.

# Examples
```jldoctest
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