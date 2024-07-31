@testitem "ORFI instances" begin
    # cd(@__DIR__)
    using BioSequences, GeneFinder
    # A random seq to start

    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

    orflarge = ORFI{4,NaiveFinder}("seq", 1, 33, STRAND_POS, 1, convert(LongSubSeq, dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), NamedTuple())
    orfshtort = ORFI{4,NaiveFinder}("seq", 16,33, STRAND_POS, 1, _orfseq(seq, 16, 33, STRAND_POS), NamedTuple()) # requires _orient function and source sequence

    @test orflarge.seq == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    @test orfshtort.seq == dna"ATGCTAGTAACTAGCTAG"
    
    @test sequence(orflarge) == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
    @test sequence(orfshtort) == dna"ATGCTAGTAACTAGCTAG"

    @test getindex(seq, orflarge.first:orflarge.last) == orflarge.seq
    @test getindex(seq, orfshtort.first:orfshtort.last) == orfshtort.seq

    @test translate(orflarge) == aa"MMHACMLVTS*"
    @test translate(orfshtort) == aa"MLVTS*"

    @test length(orflarge) == 33
    @test length(orfshtort) == 18
    
end

@testitem "NaiveFinder overlaps" begin
    using BioSequences, GeneFinder

    seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs = findorfs(seq)

    orfstest = [
        ORFI{4,NaiveFinder}("seq", 1, 33, STRAND_POS, 1, convert(LongSubSeq, dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), NamedTuple()),
        ORFI{4,NaiveFinder}("seq", 4, 33, STRAND_POS, 1, convert(LongSubSeq, dna"ATGCATGCATGCATGCTAGTAACTAGCTAG"), NamedTuple()),
        ORFI{4,NaiveFinder}("seq", 8, 22, STRAND_POS, 2, convert(LongSubSeq, dna"ATGCATGCATGCTAG"), NamedTuple()),
        ORFI{4,NaiveFinder}("seq", 12,29, STRAND_POS, 3, convert(LongSubSeq, dna"ATGCATGCTAGTAACTAG"), NamedTuple()),
        ORFI{4,NaiveFinder}("seq", 16,33, STRAND_POS, 1, convert(LongSubSeq, dna"ATGCTAGTAACTAGCTAG"), NamedTuple()) #
    ]
    
    @test println(orfs[1]) == println(orfstest[1])
    @test length(orfs) == 5
end

@testitem "NaiveFinder outbounds" begin
    # Currently not supported by the NaiveFinder method, but it should be implemented in the future
    using BioSequences, GeneFinder

    # > 180195.SAMN03785337.LFLS01000089 (726nt) -> finds with minlen=30:
    # Prodigal (from Pyrodigal tests): 1 gene
    # NCBI ORFIfinder: 1 ORFI (start = 3 stop = 452, strand = -, frame = ?) = 450nt (no stop codon)
    seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs = findorfs(seq, finder=NaiveFinder, minlen=30)

    @test length(orfs) == 0

    orfstest = Vector{ORFI{4,NaiveFinder}}()
    @test println(orfs) == println(orfstest)

end

# @testitem "NaiveFinder alternative starts" begin

#     using BioSequences, GeneFinder
#     # From pyrodigal issue #13 link: https://github.com/althonos/pyrodigal/blob/1f939b0913b48dbaa55d574b20e124f1b8323825/pyrodigal/tests/test_orf_finder.py#L271
#     # Pyrodigal predicts 2 genes from this sequence:
#     # 1) An alternative start codon (GTG) sequence at 48:347
#     # 2) A common start codon sequence at 426:590
#     # On the other hand, the NCBI ORFIfinder program predicts 13 ORFIs whose length is greater than 75 nt and allowing alternative starts, from which one has an "outbound" stop
#     seq = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"
#     orfs = findorfs(seq, finder=NaiveFinder, minlen=75, alternative_start=true)
#     @test length(orfs) == 26 #20 

#     orfstest = [
#                          "seq", # |- oubounds notation for the NCBI ORFIfinder 
#         #                "seq", (>3:890, '-', 3)     NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 37:156, '+', 1),   # NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 48:347, '+', 3),   # Pyrodigal & NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 67:156, '+', 1),   
#         #               ("seq", 153:347, '+', 3)     NCBI-ORFIfinder 
#         ORFI{NaiveFinder}("seq", 126:347, '+', 3),  
#         ORFI{NaiveFinder}("seq", 182:289, '+', 2),  # NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 194:268, '-', 2),   
#         ORFI{NaiveFinder}("seq", 194:283, '-', 2),  # NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 249:347, '+', 3),  
#         ORFI{NaiveFinder}("seq", 286:375, '+', 1),  
#         #               ("seq", 405:590, '+', 3)     NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 426:590, '+', 3),  # Pyrodigal
#         ORFI{NaiveFinder}("seq", 446:520, '-', 2),
#         ORFI{NaiveFinder}("seq", 446:523, '-', 2),
#         #               ("seq", 538:657, '+', 1)     NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 565:657, '+', 1),
#         ORFI{NaiveFinder}("seq", 650:727, '-', 2),
#         #               ("seq", 675:872, '+', 3)     NCBI-ORFIfinder
#         ORFI{NaiveFinder}("seq", 698:820, '+', 2),
#         ORFI{NaiveFinder}("seq", 746:820, '+', 2),
#         ORFI{NaiveFinder}("seq", 786:872, '+', 3),
#         ORFI{NaiveFinder}("seq", 793:876, '+', 1),
#         ORFI{NaiveFinder}("seq", 802:876, '+', 1),
#         ORFI{NaiveFinder}("seq", 887:976, '-', 2),
#     ]

#     @test println(orfs) == println(orfstest)
# end

@testitem "NaiveFinder lambda" begin
    cd(@__DIR__) # Required to find the fasta file
    using BioSequences, GeneFinder, FASTX

    # Lambda phage tests
    # Compare to https://github.com/jonas-fuchs/viral_orf_finder/blob/master/orf_finder.py 
    # Salisbury and Tsorukas (2019) paper used the Lambda phage genome with 73 CDS and 545 non-CDS ORFIs (a total of 618).
    
    # For a minimal length of 75 nt the following ORFIs are predicted: 
    # orf_finder.py --> 885 (222 complete)
    # findorfs (GeneFinder.jl) --> 885
    # NCBI ORFIfinder --> 375 ORFIs
    # orfipy --> 375 (`orfipy NC_001416.1.fasta --start ATG --include-stop --min 75`)
    
    lambdafile = "data/NC_001416.1.fasta"
    lambdarecord = open(collect, FASTA.Reader, lambdafile)[1]
    lambdaseq = FASTX.sequence(LongDNA{4}, lambdarecord)
    lambdaorfs = findorfs(lambdaseq, finder=NaiveFinder, minlen=75)
    
    @test length(lambdaorfs) == 885
end