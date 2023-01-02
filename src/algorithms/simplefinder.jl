using BioSequences
using TestItems
include("../types.jl")

"""
    simplefinder(sequence::LongDNA)

The simplest algorithm that finds ORFs in a DNA sequence.

The simplefinder function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
    This function has not ORFs size and overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function simplefinder(sequence::LongDNA)
    orf = nothing
    orfs = Vector{ORF}()
    startcodon = ExactSearchQuery(Codon("ATG"), iscompatible)
    seqbound = length(sequence) - 3

    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence

        start_codon_indices = findall(startcodon, seq)

        for i in start_codon_indices
            for j in i.start:3:seqbound
                if seq[j:j+2] ∈ stopcodons
                    push!(orfs, orf)
                    break
                end
                orf = ORF(i.start:j+5, strand)
            end
        end
    end
    return orfs
end

@testitem "simplefinder test" default_imports = true begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs01 = simplefinder(seq01)

    @test simplefinder(seq01) == [ORF(1:33, '+'), ORF(4:33, '+'), ORF(8:22, '+'), ORF(12:29, '+'), ORF(16:33, '+')]
    @test length(orfs01) == 5

    seq02 = dna"GATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs02 = simplefinder(seq02)

    @test simplefinder(seq02) == [ORF(2:34, '+'), ORF(5:34, '+'), ORF(9:23, '+'), ORF(13:30, '+'), ORF(17:34, '+')]
    @test length(orfs02) == 5
end


"""
    findcds(sequence::LongDNA)

A function to generete CDSs sequence out of a DNA sequence.

The `findcds` is a generator function that takes a `LongDNA` sequence and returns an iterator over the given sequence,
    containing the coding sequences (CDSs) found in the sequence and the ORF. 
    It uses the `simplefinder` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.
"""
function findcds(sequence::LongDNA)
    orfs = simplefinder(sequence)
    reversedseq = reverse_complement(sequence)
    cds = (i.strand == '+' ? sequence[i.location] : reversedseq[i.location] for i in orfs)
    return collect(cds)
end

# function findcds(sequence::LongDNA)

#     orfs = simplefinder(sequence)
#     seqs = Vector{CDS}()

#     @simd for i in orfs
#         if i.strand == '-'
#             reversedsequence = reverse_complement(sequence)
#             seq = reversedsequence[i.location]
#         else
#             seq = sequence[i.location]
#         end
#         cds = CDS(i, seq)
#         push!(seqs, cds)
#     end
#     return seqs
# end

@testitem "findcds test" default_imports = true begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    cds01 = findcds(seq01)

    @test length(cds01) == 5
    @test cds01 == [dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCTAG", dna"ATGCATGCTAGTAACTAG", dna"ATGCTAGTAACTAGCTAG"]
    @test cds01[1] == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
end

"""
    findproteins(sequence::LongDNA)

As its name suggest this function generate the possible proteins directly from a DNA sequence. 
    The `findcds` function takes a `LongDNA` sequence and returns a `Vector{CDS}` containing the 
    coding sequences (CDSs) found in the sequence. 
"""
function findproteins(sequence::LongDNA)
    orfs = simplefinder(sequence)
    reversedseq = reverse_complement(sequence)
    protein = (i.strand == '+' ? translate(sequence[i.location]) : translate(reversedseq[i.location]) for i in orfs)
    return collect(protein)
end

# function findproteins(sequence::LongDNA)
#     cds = findcds(sequence)
#     proteins = Vector{Protein}()
#     @simd for i in cds
#         proteinseq = translate(i.sequence)
#         protein = Protein(i.orf, proteinseq)
#         push!(proteins, protein)
#     end
#     return proteins
# end

@testitem "findproteins test" default_imports = true begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    proteins01 = findproteins(seq01)
    
    @test length(proteins01) == 5
    @test proteins01 == [aa"MMHACMLVTS*", aa"MHACMLVTS*", aa"MHAC*", aa"MHASN*", aa"MLVTS*"]
    @test proteins01[1] == aa"MMHACMLVTS*"

end