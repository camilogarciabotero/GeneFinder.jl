using BioSequences
using TestItems
include("../types.jl")
include("../helpers.jl")

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
    seqbound = length(sequence) - 2 #3

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

@testitem "simplefinder test" begin
    using BioSequences

    # A random seq to start
    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs01 = simplefinder(seq01)

    @test simplefinder(seq01) == [ORF(1:33, '+'), ORF(4:33, '+'), ORF(8:22, '+'), ORF(12:29, '+'), ORF(16:33, '+')]
    @test length(orfs01) == 5

    # > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
    seq02 = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs02 = simplefinder(seq02)

    @test length(orfs02) == 12
    @test simplefinder(seq02) == [ORF(29:40, '+'), ORF(137:145, '+'), ORF(164:184, '+'), ORF(173:184, '+'), ORF(236:241, '+'), ORF(248:268, '+'), ORF(362:373, '+'), ORF(470:496, '+'), ORF(551:574, '+'), ORF(569:574, '+'), ORF(581:601, '+'), ORF(695:706, '+')]
end

"""
    cdsgenerator(sequence::LongDNA)

A function to generete CDSs sequence out of a DNA sequence.

The `cdsgenerator` is a generator function that takes a `LongDNA` sequence and returns an iterator over the given sequence,
    containing the coding sequences (CDSs) found in the sequence and the ORF. 
    It uses the `simplefinder` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.
"""
function cdsgenerator(sequence::LongDNA)
    orfs = simplefinder(sequence)
    reversedseq = reverse_complement(sequence)
    return(i.strand == '+' ? sequence[i.location] : reversedseq[i.location] for i in orfs)
end

function cdsgenerator(sequence::String)
    sequence = LongDNA{4}(sequence)
    orfs = simplefinder(sequence)
    reversedseq = reverse_complement(sequence)
    return(i.strand == '+' ? sequence[i.location] : reversedseq[i.location] for i in orfs)
end

@testitem "cdsgenerator test" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    cds01 = collect(cdsgenerator(seq01))

    @test length(cds01) == 5
    @test cds01 == [dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCTAG", dna"ATGCATGCTAGTAACTAG", dna"ATGCTAGTAACTAGCTAG"]
end

"""
    proteingenerator(sequence::LongDNA)

As its name suggest this generator function that iterates over the sequence to find proteins directly from a DNA sequence. 
    The `cdsgenerator` function takes a `LongDNA` sequence and returns a `Vector{CDS}` containing the 
    coding sequences (CDSs) found in the sequence. 
"""
function proteingenerator(sequence::LongDNA)
    orfs = simplefinder(sequence)
    reversedseq = reverse_complement(sequence)
    return(i.strand == '+' ? translate(sequence[i.location]) : translate(reversedseq[i.location]) for i in orfs)
end

function proteingenerator(sequence::String)
    sequence = LongDNA{4}(sequence)
    orfs = simplefinder(sequence)
    reversedseq = reverse_complement(sequence)
    return(i.strand == '+' ? translate(sequence[i.location]) : translate(reversedseq[i.location]) for i in orfs)
end

@testitem "proteingenerator test" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    proteins01 = collect(proteingenerator(seq01))
    
    @test length(proteins01) == 5
    @test proteins01 == [aa"MMHACMLVTS*", aa"MHACMLVTS*", aa"MHAC*", aa"MHASN*", aa"MLVTS*"]
end