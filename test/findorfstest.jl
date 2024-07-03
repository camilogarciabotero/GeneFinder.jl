
@testitem "findorfs tests" default_imports=false begin
    # cd(@__DIR__)
    # using GeneFinder: findorfs, NaiveFinder, ORF
    # A random seq to start
    using Test, BioSequences, GeneFinder
    
    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

    orf01 = ORF{4,NaiveFinder}("seq01", 1, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing) # ORF{4,NaiveFinder}("phi", 7, 15, STRAND_POS, 1, Features((score = 0,)))
    # source(orf01)
    # GeneFinder.source(orf01)
    orfs01 = findorfs(seq01)
    # @view(seq01[begin:end])
    getindex(seq01, orf01.first:orf01.last)

    orfs01test = [
        ORF{4,NaiveFinder}("seq01", 1, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 4, 33, STRAND_POS, 1, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 8, 22, STRAND_POS, 2, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 12,29, STRAND_POS, 3, Features((score = 0.0,)), nothing),
        ORF{4,NaiveFinder}("seq01", 16,33, STRAND_POS, 1, Features((score = 0.0,)), nothing)
    ]
    
    @test println(orfs01) == println(orfs01test)
    @test length(orfs01) == 5

    # known failures
    # seq01[orf01]
    # sequence(orf01) == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    # groupname(orf01) == "seq01"

    # > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
    seq02 = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs02 = findorfs(seq02, finder=NaiveFinder)

    @test length(orfs02) == 12

    orfs02test = [
        ORF{4,NaiveFinder}("seq02", 29, 40, STRAND_POS, 2 , Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 137,145, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 164,184, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 173,184, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 236,241, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 248,268, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 362,373, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 470,496, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 551,574, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 569,574, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 581,601, STRAND_POS, 2, Features((score = 0.0,)), nothing), 
        ORF{4,NaiveFinder}("seq02", 695,706, STRAND_POS, 2, Features((score = 0.0,)), nothing)
    ]
    
    @test println(orfs02) == println(orfs02test)

    # From pyrodigal issue #13 link: https://github.com/althonos/pyrodigal/blob/1f939b0913b48dbaa55d574b20e124f1b8323825/pyrodigal/tests/test_orf_finder.py#L271
    # Pyrodigal predicts 2 genes from this sequence:
    # 1) An alternative start codon (GTG) sequence at 48:347
    # 2) A common start codon sequence at 426:590
    # On the other hand, the NCBI ORFfinder program predicts 9 ORFs whose length is greater than 75 nt, from which one has an "outbound" stop
    # seq03 = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"
    # orfs03 = findorfs(seq03, NaiveFinder(), minlen=75)
    # @test length(orfs03) == 9
    # @test orfs03 == [ORF{NaiveFinder}(37:156, STRAND_POS, 1), ORF{NaiveFinder}(194:268, '-', 2), ORF{NaiveFinder}(194:283, '-', 2), ORF{NaiveFinder}(249:347, STRAND_POS, 3), ORF{NaiveFinder}(426:590, STRAND_POS, 3), ORF{NaiveFinder}(565:657, STRAND_POS, 1), ORF{NaiveFinder}(650:727, '-', 2), ORF{NaiveFinder}(786:872, STRAND_POS, 3), ORF{NaiveFinder}(887:976, '-', 2)]
                                                                                                           #|->  This occured in Pyrodigal
    # Lambda phage tests
    # Compare to https://github.com/jonas-fuchs/viral_orf_finder/blob/master/orf_finder.py 
    # Salisbury and Tsorukas (2019) paper used the Lambda phage genome with 73 CDS and 545 non-CDS ORFs (a total of 618) to compare predictions between several Gene Finder programs
    # For a minimal length of 75 nt the following ORFs are predicted: 
    # orf_finder.py --> 885 (222 complete)
    # findorfs (GeneFinder.jl) --> 885
    # NCBI ORFfinder --> 375 ORFs
    # orfipy --> 375 (`orfipy NC_001416.1.fasta --start ATG --include-stop --min 75`)
    # NC_001416 = fasta2bioseq("data/NC_001416.1.fasta")[1]
    # NC_001416_orfs = findorfs(NC_001416, NaiveFinder(), minlen=75)
    # @test length(NC_001416_orfs) == 885
end