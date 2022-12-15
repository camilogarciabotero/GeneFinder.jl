using BioSequences

"""
    struct ORF
        position::UnitRange{Int64}
        strand::Char
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: position, which is a UnitRange{Int64} indicating the start and end positions of the ORF in the sequence, and strand, which is a Char indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
"""
struct ORF
    position::UnitRange{Int64}
    strand::Char
end

const stopcodons = [dna"TAG", dna"TAA", dna"TGA"]

"""
    simplefinder(sequence::LongDNA)

The simplest algorithm that finds ORFs in a DNA sequence.

The simplefinder function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.

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
                if seq[j:j+2] in stopcodons && orf != nothing
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