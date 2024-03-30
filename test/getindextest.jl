# Import the function to be tested

# Run the tests
@testset "getindex" begin
    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    
    orfs01 = [ORF(1:33, '+', 1, 0.0), ORF(4:33, '+', 1, 0.0), ORF(8:22, '+', 2, 0.0), ORF(12:29, '+', 3, 0.0), ORF(16:33, '+', 1, 0.0)]

    @test getindex(seq01, orfs01[1]) == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    @test getindex(seq01, orfs01[2]) == dna"ATGCATGCATGCATGCTAGTAACTAGCTAG"
end