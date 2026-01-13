@testset "Writing methods" begin
    cd(@__DIR__)
    using BioSequences, GeneFinder, FASTX

    seq = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"
    
    # Test ORFCollection API
    collection = NaiveFinder(seq)
    
    @test length(collection) > 0
    @test !isempty(collection)
    
    # Test sequence extraction
    first_seq = GeneFinder.sequence(collection, 1)
    @test first_seq isa LongSubSeq{DNAAlphabet{4}}
    
    # Test iteration
    for orf in collection
        @test orf isa ORF
        @test leftposition(orf) >= 1
        @test rightposition(orf) <= length(seq)
    end

    # Test FNA writing
    seqfna = "data/out-seq.fna"
    open(seqfna, "w") do io
        write_orfs_fna(seq, io; finder=NaiveFinder)
    end
    
    fnarecords = open(collect, FASTA.Reader, seqfna)
    @test length(fnarecords) == length(collection)
    @test !isempty(FASTX.identifier(fnarecords[1]))
    
    # Test FAA writing
    seqfaa = "data/out-seq.faa"
    open(seqfaa, "w") do io
        write_orfs_faa(seq, io; finder=NaiveFinder)
    end
    
    faarecords = open(collect, FASTA.Reader, seqfaa)
    @test length(faarecords) == length(collection)
    
    # Test BED writing
    seqbed = "data/out-seq.bed"
    open(seqbed, "w") do io
        write_orfs_bed(seq, io; finder=NaiveFinder)
    end
    
    bedlines = readlines(seqbed)
    @test length(bedlines) == length(collection)
    
    # Test GFF writing
    seqgff = "data/out-seq.gff"
    open(seqgff, "w") do io
        write_orfs_gff(seq, io; finder=NaiveFinder)
    end
    
    gfflines = readlines(seqgff)
    # GFF has 2 header lines + ORF entries
    @test length(gfflines) == length(collection) + 2
    @test startswith(gfflines[1], "##gff-version")

    # # Cleanup
    rm(seqfna, force=true)
    rm(seqfaa, force=true)
    rm(seqbed, force=true)
    rm(seqgff, force=true)
end