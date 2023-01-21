# # this will be the main functions taking all the algorithms

"""
    cds_generator(sequence::LongDNA; kwargs...)
    cds_generator(sequence::String; kwargs...)

A function to generete CDSs sequence out of a DNA sequence.

The `cds_generator` is a generator function that takes a `LongDNA` sequence and returns an iterator over the given sequence,
    containing the coding sequences (CDSs) found in the sequence and the ORF. 
    It uses the `findorfs` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF, returining both. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function cds_generator(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    cds = (i.strand == '+' ? CDS(sequence[i.location], i) : CDS(reversedseq[i.location], i) for i in orfs)
    return cds
end

function cds_generator(sequence::String; alternative_start::Bool=false, min_len::Int64=6)
    sequence = LongDNA{4}(sequence)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    cds = (i.strand == '+' ? CDS(sequence[i.location], i) : CDS(reversedseq[i.location], i) for i in orfs)
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

function get_cds(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    orfs = orfiterator(sequence; alternative_start, min_len)
    revseq = reverse_complement(sequence)
    cds = [i.strand == '+' ? StringView(@view sequence[i.location]) : StringView(@view revseq[i.location]) for i in orfs]
    return cds
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
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (i.strand == '+' ? Protein(translate(sequence[i.location]; alternative_start, code), i) : Protein(translate(reversedseq[i.location]; alternative_start, code), i) for i in orfs)
    return proteins
end

function protein_generator(sequence::String; alternative_start::Bool=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len::Int64=6)
    sequence = LongDNA{4}(sequence)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (i.strand == '+' ? Protein(translate(sequence[i.location]; alternative_start, code), i) : Protein(translate(reversedseq[i.location]; alternative_start, code), i) for i in orfs)
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


function get_proteins(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    orfs = orfiterator(sequence; alternative_start, min_len)
    revseq = reverse_complement(sequence)
    proteins = [i.strand == '+' ? translate(sequence[i.location]) : translate(revseq[i.location]) for i in orfs]
    return proteins
end


# """
# FindGene struct
# """
# mutable struct FindGene{S1,S2}
#     a::GeneFinderAlgorithm{S1}
#     b::LongDNA
#     c::GeneticCode
# end

# function findgenes(::SimpleFinder, sequence::LongDNA) # type::GeneticCode
#     orfs = simplefinder(sequence)
#     seqs = Vector{CDS}()
# end