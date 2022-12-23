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
    orfs = Vector{ORF}()
    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence
        
        start_codon = ExactSearchQuery(dna"ATG", iscompatible)
        start_codon_indices = findall(start_codon, seq)

        for i in start_codon_indices
            orf = nothing
            j = i.start
            while j < length(seq) - 3
                if seq[j:j+2] ∈ stopcodons
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
    orfs = simplefinder(seq)

    @test simplefinder(seq) == orfs
    @test length(orfs) == 5
end

"""
    findcds(sequence::LongDNA)

A function to generete CDSs sequence out of a DNA sequence.

The `findcds` function takes a `LongDNA` sequence and returns a `Vector{CDS}` 
    containing the coding sequences (CDSs) found in the sequence. 
    It uses the `simplefinder` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.
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
        cds = CDS(i, seq)
        push!(seqs, cds)
    end
    return seqs
end


@testitem "findcds test" default_imports=true begin
    using BioSequences
    
    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    cds = findcds(seq)
    # @test findcds(seq) == [CDS(1:33, '+', dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), CDS(4:33, '+', dna"ATGCATGCATGCATGCTAGTAACTAGCTAG"), CDS(8:22, '+', dna"ATGCATGCATGCTAG"), CDS(12:29, '+', dna"ATGCATGCTAGTAACTAG"), CDS(16:33, '+', dna"ATGCTAGTAACTAGCTAG")]

    @test cds[1].orf == ORF(1:33, '+')
    @test cds[1].sequence == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    
    @test length(cds) == 5
end


"""
    findproteins(sequence::LongDNA)

As its name suggest this function generate the possible proteins directly from a DNA sequence. 
    The `findcds` function takes a `LongDNA` sequence and returns a `Vector{CDS}` containing the 
    coding sequences (CDSs) found in the sequence. 
"""
function findproteins(sequence::LongDNA)
    cds = findcds(sequence)
    proteins = Vector{Protein}()
    for i in cds
        proteinseq = translate(i.sequence)
        protein = Protein(i.orf, proteinseq)
        push!(proteins, protein)
    end
    return proteins
end


@testitem "findproteins test" default_imports=true begin
    using BioSequences
    
    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    proteins = findproteins(seq)
    # @test findcds(seq) == [CDS(1:33, '+', dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), CDS(4:33, '+', dna"ATGCATGCATGCATGCTAGTAACTAGCTAG"), CDS(8:22, '+', dna"ATGCATGCATGCTAG"), CDS(12:29, '+', dna"ATGCATGCTAGTAACTAG"), CDS(16:33, '+', dna"ATGCTAGTAACTAGCTAG")]
    @test proteins[1].orf == ORF(1:33, '+')
    @test proteins[1].sequence == aa"MMHACMLVTS*"

    @test length(proteins) == 5
end