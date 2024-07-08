@testitem "writting methods" begin
    cd(@__DIR__) # Required to find the fasta file    
    using BioSequences, GeneFinder, FASTX
# 
    # Test case 1 fna
# 
    # From pyrodigal issue #13 link: https://github.com/althonos/pyrodigal/blob/1f939b0913b48dbaa55d574b20e124f1b8323825/pyrodigal/tests/test_orf_finder.py#L271
    # Pyrodigal predicts 2 genes from this sequence:
    # 1) An alternative start codon (GTG) sequence at 48:347
    # 2) A common start codon sequence at 426:590
    # On the other hand, the NCBI ORFfinder program predicts 9 ORFs whose length is greater than 75 nt, from which one has an "outbound" stop
    seq = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"
    # 
    seqfna = "data/out-seq.fna"
    open(seqfna, "w") do io
        write_orfs_fna(seq, io, finder=NaiveFinder)
    end
# 
    seqfnarecords = open(collect, FASTAReader, "data/out-seq.fna")
# 
    @test seqfnarecords[1] == FASTX.FASTA.Record("unnamedseq id=01 start=5 stop=22 strand=- frame=2 score=0.0", "ATGAAACAGAACGACTGA")
    @test length(seqfnarecords) == 32
    @test identifier(seqfnarecords[1]) == "unnamedseq"
    @test description(seqfnarecords[1]) == "unnamedseq id=01 start=5 stop=22 strand=- frame=2 score=0.0"
    @test FASTX.sequence(seqfnarecords[1]) == "ATGAAACAGAACGACTGA"
# 
    # Test case 2 faa

    seqfaa = "data/out-seq.faa"
    open(seqfaa, "w") do io
        write_orfs_faa(seq, io, finder=NaiveFinder)
    end
# 
    seqfaarecords = open(collect, FASTAReader, "data/out-seq.faa")
# 
    @test seqfaarecords[2] == FASTX.FASTA.Record("unnamedseq id=02 start=37 stop=156 strand=+ frame=1 score=0.0", "MYFSCISGGIVPSSIASRYPLIEVRGERKSWDTLATNFF*")
    @test length(seqfaarecords) == 32
    @test identifier(seqfaarecords[2]) == "unnamedseq"
    @test description(seqfaarecords[2]) == "unnamedseq id=02 start=37 stop=156 strand=+ frame=1 score=0.0"
    @test FASTX.sequence(seqfaarecords[2]) == "MYFSCISGGIVPSSIASRYPLIEVRGERKSWDTLATNFF*"

end