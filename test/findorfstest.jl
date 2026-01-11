using Test
using BioSequences, GeneFinder, FASTX

# @testset "ORF instances" begin
#     seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    
#     orf_large = ORF{NaiveFinder}(:seq, 1:33, GeneFinder.PSTRAND, 1, NamedTuple())
#     orf_short = ORF{NaiveFinder}(:seq, 16:33, GeneFinder.PSTRAND, 1, NamedTuple())
    
#     @test typeof(GeneFinder.sequence(orf_large)) == typeof(convert(LongSubSeq{DNAAlphabet{4}}, dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"))
#     @test typeof(GeneFinder.sequence(orf_short)) == typeof(convert(LongSubSeq{DNAAlphabet{4}}, dna"ATGCTAGTAACTAGCTAG"))
# end

@testset "NaiveFinder overlaps" begin
    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    
    orfs = findorfs(seq)

    @test length(orfs) == 5
    
    @test leftposition(orfs[1]) == 1
    @test rightposition(orfs[1]) == 33
    @test strand(orfs[1]) === GeneFinder.PSTRAND
    @test frame(orfs[1]) == 1
    
    @test leftposition(orfs[5]) == 16
    @test rightposition(orfs[5]) == 33
    @test strand(orfs[5]) === GeneFinder.PSTRAND
    @test frame(orfs[5]) == 1
end

@testset "NaiveFinder outbounds" begin
    seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    
    orfs = findorfs(seq, minlen=30)

    @test length(orfs) == 0
end

@testset "NaiveFinder lambda" begin
    lambdafile = "data/NC_001416.1.fasta"
    lambdarecord = open(collect, FASTA.Reader, lambdafile)[1]
    lambdaseq = FASTX.sequence(LongDNA{4}, lambdarecord)
    
    lambdaorfs = findorfs(lambdaseq, finder=NaiveFinder, minlen=75)
    
    @test length(lambdaorfs) == 885
    
    # @test all(length(GeneFinder.sequence(orf)) >= 75 for orf in lambdaorfs)
    @test all(GeneFinder.strand(orf) in (GeneFinder.PSTRAND, GeneFinder.NSTRAND) for orf in lambdaorfs)
    @test all(GeneFinder.frame(orf) in (1, 2, 3) for orf in lambdaorfs)
end