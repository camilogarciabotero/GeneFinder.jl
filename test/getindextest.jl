# Import the function to be tested

# Run the tests
@testitem "getindex" default_imports=false begin
    
    using BioSequences: @dna_str
    using GeneFinder: ORF, NaiveFinder
    using Test: @test

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    
    orfs01 = [ORF("seq01", 1:33, '+', 1, NaiveFinder()), ORF("seq01", 4:33, '+', 1, NaiveFinder()), ORF("seq01", 8:22, '+', 2, NaiveFinder()), ORF("seq01", 12:29, '+', 3, NaiveFinder()), ORF("seq01", 16:33, '+', 1, NaiveFinder())]

    @test getindex(seq01, orfs01[1]) == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    @test getindex(seq01, orfs01[2]) == dna"ATGCATGCATGCATGCTAGTAACTAGCTAG"
end