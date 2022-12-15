using BioSequences
struct ORF
    position::UnitRange{Int64}
    strand::Char
end

const stopcodons = [dna"TAG", dna"TAA", dna"TGA"]

function simplefinder(sequence::LongDNA; orfsize::Int64=12)
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