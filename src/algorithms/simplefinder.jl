"""
    orf_finder(sequence::LongDNA; kwargs...)::Vector{ORF}
    orf_finder(sequence::String; kwargs...)::Vector{ORF} 

The simplest algorithm that finds ORFs in a DNA sequence.

The `orf_finder` function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
    This function has not ORFs size and overlapping condition contraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.
        Extending the starting codons with the `alternative_start = true` will search for ATG, GTG, and TTG.
    Some studies have shown that in *E. coli* (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function orf_finder(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    orf = nothing
    orfs = Vector{ORF}()
    seqbound = length(sequence) - 2 #3
    startcodons = alternative_start ? EXTENDED_STARTCODONS : STARTCODON

    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence
        start_codon_indices = findall(startcodons, seq)

        for i in start_codon_indices
            for j in i.start:3:seqbound
                if seq[j:j+2] ∈ STOPCODONS
                    push!(orfs, orf)
                    break
                end
                orf = ORF(i.start:j+5, strand)
            end
        end
    end
    return filter(i -> length(i.location) >= min_len, orfs)
end

function orf_finder(sequence::String; alternative_start::Bool=false, min_len::Int64=6)
    sequence = LongDNA{4}(sequence)
    orf = nothing
    orfs = Vector{ORF}()
    seqbound = length(sequence) - 2 #3
    startcodons = alternative_start ? EXTENDED_STARTCODONS : STARTCODON

    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence
        start_codon_indices = findall(startcodons, seq)

        for i in start_codon_indices
            for j in i.start:3:seqbound
                if seq[j:j+2] ∈ STOPCODONS
                    push!(orfs, orf)
                    break
                end
                orf = ORF(i.start:j+5, strand)
            end
        end
    end
    return filter(i -> length(i.location) >= min_len, orfs)
end

@testitem "orf_finder test" begin
    using BioSequences

    # A random seq to start
    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs01 = orf_finder(seq01)

    @test orf_finder(seq01) == [ORF(1:33, '+'), ORF(4:33, '+'), ORF(8:22, '+'), ORF(12:29, '+'), ORF(16:33, '+')]
    @test length(orfs01) == 5

    # > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
    seq02 = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs02 = orf_finder(seq02)

    @test length(orfs02) == 12
    @test orf_finder(seq02) == [ORF(29:40, '+'), ORF(137:145, '+'), ORF(164:184, '+'), ORF(173:184, '+'), ORF(236:241, '+'), ORF(248:268, '+'), ORF(362:373, '+'), ORF(470:496, '+'), ORF(551:574, '+'), ORF(569:574, '+'), ORF(581:601, '+'), ORF(695:706, '+')]
end

"""
    cds_generator(sequence::LongDNA; kwargs...)
    cds_generator(sequence::String; kwargs...)

A function to generete CDSs sequence out of a DNA sequence.

The `cds_generator` is a generator function that takes a `LongDNA` sequence and returns an iterator over the given sequence,
    containing the coding sequences (CDSs) found in the sequence and the ORF. 
    It uses the `orf_finder` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF, returining both. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function cds_generator(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    orfs = orf_finder(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    cds = (i.strand == '+' ? CDS(sequence[i.location], i) : CDS(reversedseq[i.location], i) for i in orfs)
    return cds
end

function cds_generator(sequence::String; alternative_start::Bool=false, min_len::Int64=6)
    sequence = LongDNA{4}(sequence)
    orfs = orf_finder(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    cds = (i.strand == '+' ? CDS(sequence[i.location], i) : CDS(reversedseq[i.location], i) for i in orfs)
    return cds
end

function get_cds(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    cds = [i.sequence for i in cds_generator(sequence; alternative_start, min_len)]
    return cds
end

@testitem "cds_generator test" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    cds01 = collect(cds_generator(seq01))

    @test length(cds01) == 5
    # @test cds01 == [dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCTAG", dna"ATGCATGCTAGTAACTAG", dna"ATGCTAGTAACTAGCTAG"]
    # @test cds01[1] == [CDS(dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG", ORF(1:33, '+')), CDS(dna"ATGCATGCATGCATGCTAGTAACTAGCTAG", ORF(4:33, '+')), CDS(dna"ATGCATGCATGCTAG", ORF(8:22, '+')), CDS(dna"ATGCATGCTAGTAACTAG", ORF(12:29, '+')), CDS(dna"ATGCTAGTAACTAGCTAG", ORF(16:33, '+'))]
    @test cds01[1].sequence == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    @test cds01[1].orf == ORF(1:33, '+')
end

"""
    protein_generator(sequence::LongDNA; kwargs...)
    protein_generator(sequence::String; kwargs...)

As its name suggest this generator function iterates over the sequence to find proteins directly from a DNA sequence. 
    The `protein_generator` function takes a `LongDNA` sequence and returns a `Vector{CDS}` containing the 
    coding sequences (CDSs) found in the sequence and the associated ORF.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function protein_generator(sequence::LongDNA; alternative_start::Bool=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len::Int64=6)
    orfs = orf_finder(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (i.strand == '+' ? Protein(translate(sequence[i.location]; alternative_start, code), i) : Protein(translate(reversedseq[i.location]; alternative_start, code), i) for i in orfs)
    return proteins
end

function protein_generator(sequence::String; alternative_start::Bool=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len::Int64=6)
    sequence = LongDNA{4}(sequence)
    orfs = orf_finder(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (i.strand == '+' ? Protein(translate(sequence[i.location]; alternative_start, code), i) : Protein(translate(reversedseq[i.location]; alternative_start, code), i) for i in orfs)
    return proteins
end

function get_proteins(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    proteins = [i.sequence for i in protein_generator(sequence; alternative_start, min_len)]
    return proteins
end

@testitem "protein_generator test" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    proteins01 = collect(protein_generator(seq01))

    @test length(proteins01) == 5
    # @test proteins01 == [aa"MMHACMLVTS*", aa"MHACMLVTS*", aa"MHAC*", aa"MHASN*", aa"MLVTS*"]

    @test proteins01[1].sequence == aa"MMHACMLVTS*"
    @test proteins01[1].orf == ORF(1:33, '+')
end