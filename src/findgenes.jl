# # this will be the main functions taking all the algorithms

"""
    cdsgenerator(sequence::LongSequence{DNAAlphabet{4}}; kwargs...)
    cdsgenerator(sequence::String; kwargs...)

A function to generete CDSs sequence out of a DNA sequence.

The `cdsgenerator` is a generator function that takes a `LongSequence{DNAAlphabet{4}}` sequence and returns an iterator over the given sequence,
    containing the coding sequences (CDSs) found in the sequence and the ORF. 
    It uses the `findorfs` function to find open reading frames (ORFs) in the sequence, 
    and then it extracts the actual CDS sequence from each ORF, returining both. 
    The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function cdsgenerator(
    sequence::LongSequence{DNAAlphabet{4}};
    alternative_start::Bool = false,
    min_len::Int64 = 6,
)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    cds = (
        i.strand == '+' ? CDS(@view(sequence[i.location]), i) :
        CDS(@view(reversedseq[i.location]), i) for i in orfs
    )
    return cds
end

function cdsgenerator(sequence::String; alternative_start::Bool = false, min_len::Int64 = 6)
    sequence = LongSequence{DNAAlphabet{4}}(sequence)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    cds = (
        i.strand == '+' ? CDS(@view(sequence[i.location]), i) :
        CDS(@view(reversedseq[i.location]), i) for i in orfs
    )
    return cds
end

@testitem "cdsgenerator test" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    cds01 = collect(cdsgenerator(seq01))

    @test length(cds01) == 5
    # @test cds01 == [dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCATGCTAGTAACTAGCTAG", dna"ATGCATGCATGCTAG", dna"ATGCATGCTAGTAACTAG", dna"ATGCTAGTAACTAGCTAG"]
    # @test cds01[1] == [CDS(dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG", ORF(1:33, '+')), CDS(dna"ATGCATGCATGCATGCTAGTAACTAGCTAG", ORF(4:33, '+')), CDS(dna"ATGCATGCATGCTAG", ORF(8:22, '+')), CDS(dna"ATGCATGCTAGTAACTAG", ORF(12:29, '+')), CDS(dna"ATGCTAGTAACTAGCTAG", ORF(16:33, '+'))]
    @test cds01[1].sequence == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    @test cds01[1].orf == ORF(1:33, '+')
end

"""
    getcds(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...)
    getcds(input::String, output::String; kwargs...) ## for strings per se

This function will take a `LongSequence{DNAAlphabet{4}}` or `String` sequence and by means of the `findorfs()` function will push `LongSubSeq{DNAAlphabet{4}}` into a `Vector{}`

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function getcds(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool = false, min_len::Int64 = 6)
    cds = Vector{LongSubSeq{DNAAlphabet{4}}}()
    for i in cdsgenerator(sequence; alternative_start, min_len)
        push!(cds, i.sequence)
    end
    return cds
end

function getcds(sequence::String; alternative_start::Bool = false, min_len = 6)
    dnaseq = LongSequence{DNAAlphabet{4}}(sequence)
    cds = Vector{LongSubSeq{DNAAlphabet{4}}}()
    for i in cdsgenerator(dnaseq; alternative_start, min_len)
        push!(cds, i.sequence)
    end
    return cds
end


## Need to replace by a IO input
# function getcds(sequence::String; alternative_start::Bool = false, min_len = 6)
#     dnaseq = fasta_to_dna(sequence)[1]
#     cds = Vector{LongSubSeq{DNAAlphabet{4}}}()
#     for i in cdsgenerator(dnaseq; alternative_start, min_len)
#         push!(cds, i.sequence)
#     end
#     return cds
# end

"""
    proteingenerator(sequence::LongSequence{DNAAlphabet{4}}; kwargs...)
    proteingenerator(sequence::String; kwargs...)

As its name suggest this generator function iterates over the sequence to find proteins directly from a DNA sequence. 
    The `proteingenerator` function takes a `LongSequence{DNAAlphabet{4}}` sequence and returns a `Vector{CDS}` containing the 
    coding sequences (CDSs) found in the sequence and the associated ORF.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function proteingenerator(
    sequence::LongSequence{DNAAlphabet{4}};
    alternative_start::Bool = false,
    code::GeneticCode = BioSequences.standard_genetic_code,
    min_len::Int64 = 6,
)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (
        i.strand == '+' ?
        Protein(translate(@view(sequence[i.location]); alternative_start, code), i) :
        Protein(translate(@view(reversedseq[i.location]); alternative_start, code), i) for
        i in orfs
    )
    return proteins
end

function proteingenerator(
    sequence::String;
    alternative_start::Bool = false,
    code::GeneticCode = BioSequences.standard_genetic_code,
    min_len::Int64 = 6,
)
    sequence = LongSequence{DNAAlphabet{4}}(sequence)
    orfs = findorfs(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (
        i.strand == '+' ?
        Protein(translate(@view(sequence[i.location]); alternative_start, code), i) :
        Protein(translate(@view(reversedseq[i.location]); alternative_start, code), i) for
        i in orfs
    )
    return proteins
end

@testitem "proteingenerator test" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    proteins01 = collect(proteingenerator(seq01))

    @test length(proteins01) == 5
    # @test proteins01 == [aa"MMHACMLVTS*", aa"MHACMLVTS*", aa"MHAC*", aa"MHASN*", aa"MLVTS*"]

    @test proteins01[1].sequence == aa"MMHACMLVTS*"
    @test proteins01[1].orf == ORF(1:33, '+')
end

"""
    getproteins(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...)
    getproteins(input::String, output::String; kwargs...)

Similar to `getcds()` function, it will take a `LongSequence{DNAAlphabet{4}}` or `String` sequence and by means of the `findorfs()` and the `translate()` function will push `LongAA`s into a `Vector`

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function getproteins(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool = false, min_len::Int64 = 6)
    orfs = findorfs(sequence; alternative_start, min_len)
    revseq = reverse_complement(sequence)
    proteins = [
        i.strand == '+' ? translate(@view(sequence[i.location])) :
        translate(@view(revseq[i.location])) for i in orfs
    ]
    return proteins
end

function getproteins(sequence::String; alternative_start::Bool = false, min_len = 6)
    dnaseq = LongSequence{DNAAlphabet{4}}(sequence)
    orfs = findorfs(dnaseq; alternative_start, min_len)
    revseq = reverse_complement(dnaseq)
    proteins = [
        i.strand == '+' ? translate(@view(dnaseq[i.location])) :
        translate(@view(revseq[i.location])) for i in orfs
    ]
    return proteins
end

# function getproteins(sequence::String; alternative_start::Bool = false, min_len = 6)
#     dnaseq = fasta_to_dna(sequence)[1]
#     orfs = findorfs(dnaseq; alternative_start, min_len)
#     revseq = reverse_complement(dnaseq)
#     proteins = [
#         i.strand == '+' ? translate(@view(dnaseq[i.location])) :
#         translate(@view(revseq[i.location])) for i in orfs
#     ]
#     return proteins
# end

# """
# FindGene struct
# """
# mutable struct FindGene{S1,S2}
#     a::GeneFinderAlgorithm{S1}
#     b::LongSequence{DNAAlphabet{4}}
#     c::GeneticCode
# end

# function findgenes(::SimpleFinder, sequence::LongSequence{DNAAlphabet{4}}) # type::GeneticCode
#     orfs = simplefinder(sequence)
#     seqs = Vector{CDS}()
# end
