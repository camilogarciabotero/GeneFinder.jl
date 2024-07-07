
@testitem "ORF instances" begin
    # cd(@__DIR__)
    using BioSequences, GeneFinder
    # A random seq to start

    orflarge = ORF{4,NaiveFinder}("seq01", 1, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing)
    orfshtort = ORF{NaiveFinder}(1:33, '+', 1)

    # Missing tests:
    # getindex(seq, orfs.first:orfs.last)
    # source(orf01)
    # seq01[orf01]
    # sequence(orf01) == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    # groupname(orf01) == "seq01"
    # @view(seq01[begin:end])
    # GeneFinder.source(orf01) 
end

@testitem "NaiveFinder overlaps" begin
    using BioSequences, GeneFinder

    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs = findorfs(seq)


    orfstest = [
        ORF{4,NaiveFinder}("seq01", 1, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 4, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 8, 22, STRAND_POS, 2, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 12,29, STRAND_POS, 3, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 16,33, STRAND_POS, 1, Features((score = 0.0,)), nothing)
    ]
    
    @test println(orfs) == println(orfstest)
    @test length(orfs) == 5
end

@testitem "NaiveFinder outbounds" begin
    # Currently not supported by the NaiveFinder method, but it should be implemented in the future
    using BioSequences, GeneFinder

    # > 180195.SAMN03785337.LFLS01000089 (726nt) -> finds with minlen=30:
    # Prodigal (from Pyrodigal tests): 1 gene
    # NCBI ORFfinder: 1 ORF (start = 3 stop = 452, strand = -, frame = ?) = 450nt (no stop codon)
    seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs = findorfs(seq, finder=NaiveFinder, minlen=30)

    @test length(orfs) == 0

    orfstest = Vector{ORF{4,NaiveFinder}}()
    @test println(orfs) == println(orfstest)

end

@testitem "NaiveFinder alternative starts" begin

    using BioSequences, GeneFinder
    # From pyrodigal issue #13 link: https://github.com/althonos/pyrodigal/blob/1f939b0913b48dbaa55d574b20e124f1b8323825/pyrodigal/tests/test_orf_finder.py#L271
    # Pyrodigal predicts 2 genes from this sequence:
    # 1) An alternative start codon (GTG) sequence at 48:347
    # 2) A common start codon sequence at 426:590
    # On the other hand, the NCBI ORFfinder program predicts 13 ORFs whose length is greater than 75 nt and allowing alternative starts, from which one has an "outbound" stop
    seq = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"
    orfs = findorfs(seq, finder=NaiveFinder, minlen=75, alternative_start=true)
    @test length(orfs) == 20

    orfstest = [
                        # |- oubounds notation for the NCBI ORFfinder 
        #               (>3:890, '-', 3)      NCBI-ORFfinder
        ORF{NaiveFinder}(37:156, '+', 1),   # NCBI-ORFfinder
        ORF{NaiveFinder}(48:347, '+', 3),   # Pyrodigal & NCBI-ORFfinder
        ORF{NaiveFinder}(67:156, '+', 1),   
        #               (153:347, '+', 3)     NCBI-ORFfinder 
        ORF{NaiveFinder}(126:347, '+', 3),  
        ORF{NaiveFinder}(182:289, '+', 2),  # NCBI-ORFfinder
        ORF{NaiveFinder}(194:268, '-', 2),   
        ORF{NaiveFinder}(194:283, '-', 2),  # NCBI-ORFfinder
        ORF{NaiveFinder}(249:347, '+', 3),  
        ORF{NaiveFinder}(286:375, '+', 1),  
        #               (405:590, '+', 3)     NCBI-ORFfinder
        ORF{NaiveFinder}(426:590, '+', 3),  # Pyrodigal
        ORF{NaiveFinder}(446:520, '-', 2),
        ORF{NaiveFinder}(446:523, '-', 2),
        #               (538:657, '+', 1)     NCBI-ORFfinder
        ORF{NaiveFinder}(565:657, '+', 1),
        ORF{NaiveFinder}(650:727, '-', 2),
        #               (675:872, '+', 3)     NCBI-ORFfinder
        ORF{NaiveFinder}(698:820, '+', 2),
        ORF{NaiveFinder}(746:820, '+', 2),
        ORF{NaiveFinder}(786:872, '+', 3),
        ORF{NaiveFinder}(793:876, '+', 1),
        ORF{NaiveFinder}(802:876, '+', 1),
        ORF{NaiveFinder}(887:976, '-', 2),
    ]

    @test println(orfs) == println(orfstest)
end

@testitem "NaiveFinder lambda" begin
    cd(@__DIR__) # Required to find the fasta file

    using BioSequences, GeneFinder
    # Lambda phage tests
    # Compare to https://github.com/jonas-fuchs/viral_orf_finder/blob/master/orf_finder.py 
    # Salisbury and Tsorukas (2019) paper used the Lambda phage genome with 73 CDS and 545 non-CDS ORFs (a total of 618).
    
    # For a minimal length of 75 nt the following ORFs are predicted: 
    # orf_finder.py --> 885 (222 complete)
    # findorfs (GeneFinder.jl) --> 885
    # NCBI ORFfinder --> 375 ORFs
    # orfipy --> 375 (`orfipy NC_001416.1.fasta --start ATG --include-stop --min 75`)
    
    lambda = fasta2bioseq("data/NC_001416.1.fasta")[1]
    lambdaorfs = findorfs(lambda, finder=NaiveFinder, minlen=75)
    
    @test length(lambdaorfs) == 885
end