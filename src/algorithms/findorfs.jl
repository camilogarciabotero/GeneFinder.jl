"""
    locationiterator(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool=false)

This is an iterator function that uses regular expressions to search the entire ORF (instead of start and stop codons) in a `LongSequence{DNAAlphabet{4}}` sequence.
    It uses an anonymous function that will find the first regularly expressed ORF. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.

!!! note
    As a note of the implementation we want to expand on how the ORFs are found:

    The expression `(?:[N]{3})*?` serves as the boundary between the start and stop codons. Within this expression, the character class `[N]{3}` captures exactly three occurrences of any character (representing nucleotides using IUPAC codes). This portion functions as the regular codon matches. Since it is enclosed within a non-capturing group `(?:)` and followed by `*?`, it allows for the matching of intermediate codons, but with a preference for the smallest number of repetitions. 
    
    In summary, the regular expression `ATG(?:[N]{3})*?T(AG|AA|GA)` identifies patterns that start with "ATG," followed by any number of three-character codons (represented by "N" in the IUPAC code), and ends with a stop codon "TAG," "TAA," or "TGA." This pattern is commonly used to identify potential protein-coding regions within genetic sequences.

    See more about the discussion [here](https://discourse.julialang.org/t/how-to-improve-a-generator-to-be-more-memory-efficient-when-it-is-collected/92932/8?u=camilogarciabotero)

"""
function locationiterator(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool = false)
    regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    # regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG([N]{3})*T(AG|AA|GA)?"dna
    finder(x) = findfirst(regorf, sequence, first(x) + 1) # + 3
    itr = takewhile(!isnothing, iterated(finder, findfirst(regorf, sequence)))
    return itr
end

"""
    findorfs(sequence::LongSequence{DNAAlphabet{4}}; kwargs...)::Vector{ORF}
    findorfs(sequence::String; kwargs...)::Vector{ORF} 

A simple implementation that finds ORFs in a DNA sequence.

The `findorfs` function takes a LongSequence{DNAAlphabet{4}} sequence and returns a Vector{ORF} containing the ORFs found in the sequence. 
    It searches entire regularly expressed CDS, adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
        Extending the starting codons with the `alternative_start = true` will search for ATG, GTG, and TTG.
    Some studies have shown that in *E. coli* (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.
!!! note
    This function has not ORFs scoring scheme. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.
    
# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function findorfs(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool = false, min_len::Int64 = 6)
    orfs = Vector{ORF}()
    reversedseq = reverse_complement(sequence)
    seqlen = length(sequence)

    frames = Dict(0 => 3, 1 => 1, 2 => 2)

    for strand in ('+', '-')
        seq = strand == '-' ? reversedseq : sequence

        @inbounds for location in locationiterator(seq; alternative_start)
            if length(location) >= min_len
                frame = strand == '+' ? frames[location.start % 3] : frames[(seqlen - location.stop + 1) % 3]
                push!(orfs, ORF(strand == '+' ? location : (seqlen - location.stop + 1):(seqlen - location.start + 1), strand, frame))
            end
        end
    end
    return sort(orfs)
end

function findorfs(sequence::String; alternative_start::Bool = false, min_len::Int64 = 6)
    sequence = LongSequence{DNAAlphabet{4}}(sequence)
    orfs = Vector{ORF}()
    reversedseq = reverse_complement(sequence)
    seqlen = length(sequence)

    frames = Dict(0 => 3, 1 => 1, 2 => 2)

    for strand in ('+', '-')
        seq = strand == '-' ? reversedseq : sequence

        @inbounds for location in locationiterator(seq; alternative_start)
            if length(location) >= min_len
                frame = strand == '+' ? frames[location.start % 3] : frames[(seqlen - location.stop + 1) % 3]
                push!(orfs, ORF(strand == '+' ? location : (seqlen - location.stop + 1):(seqlen - location.start + 1), strand, frame))
            end
        end
    end
    return sort(orfs)
end

@testitem "findorfs test" begin
    using BioSequences, FASTX
    cd(@__DIR__)

    # A random seq to start
    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs01 = findorfs(seq01)

    @test orfs01 == [ORF(1:33, '+', 1), ORF(4:33, '+', 1), ORF(8:22, '+', 2), ORF(12:29, '+', 3), ORF(16:33, '+', 1)]
    @test length(orfs01) == 5

    # > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
    seq02 = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs02 = findorfs(seq02)

    @test length(orfs02) == 12
    @test orfs02 == [ORF(29:40, '+', 2), ORF(137:145, '+', 2), ORF(164:184, '+', 2), ORF(173:184, '+', 2), ORF(236:241, '+', 2), ORF(248:268, '+', 2), ORF(362:373, '+', 2), ORF(470:496, '+', 2), ORF(551:574, '+', 2), ORF(569:574, '+', 2), ORF(581:601, '+', 2), ORF(695:706, '+', 2)]

    # From pyrodigal issue #13 link: https://github.com/althonos/pyrodigal/blob/1f939b0913b48dbaa55d574b20e124f1b8323825/pyrodigal/tests/test_orf_finder.py#L271
    # Pyrodigal predicts 2 genes from this sequence:
    # 1) An alternative start codon (GTG) sequence at 48:347
    # 2) A common start codon sequence at 426:590
    # On the other hand, the NCBI ORFfinder program predicts 9 ORFs whose length is greater than 75 nt, from which one has an "outbound" stop
    seq03 = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"
    orfs03 = findorfs(seq03, min_len=75)
    @test length(orfs03) == 9
    @test orfs03 == [ORF(37:156, '+', 1), ORF(194:268, '-', 2), ORF(194:283, '-', 2), ORF(249:347, '+', 3), ORF(426:590, '+', 3), ORF(565:657, '+', 1), ORF(650:727, '-', 2), ORF(786:872, '+', 3), ORF(887:976, '-', 2)]
                                                                                                           #|->  This occured in Pyrodigal
    # Lambda phage tests
    # Compare to https://github.com/jonas-fuchs/viral_orf_finder/blob/master/orf_finder.py 
    # Salisbury and Tsorukas (2019) paper used the Lambda phage genome with 73 CDS and 545 non-CDS ORFs (a total of 618) to compare predictions between several Gene Finder programs
    # For a minimal length of 75 nt the following ORFs are predicted: 
    # orf_finder.py --> 885 (222 complete)
    # findorfs (GeneFinder.jl) --> 885
    # NCBI ORFfinder --> 375 ORFs
    # orfipy --> 375 (`orfipy NC_001416.1.fasta --start ATG --include-stop --min 75`)
    NC_001416 = fasta_to_dna("../../test/data/NC_001416.1.fasta")[1]

    NC_001416_orfs = findorfs(NC_001416, min_len=75)
    @test length(NC_001416_orfs) == 885
end


"""
    getorfdna(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...)
    getorfdna(input::String; kwargs...) ## for strings per se

This function takes a `LongSequence{DNAAlphabet{4}}` or `String` sequence and identifies the open reading frames (ORFs) using the `findorfs()` function. The function then extracts the DNA sequence of each ORF and stores it in a `Vector{LongSubSeq{DNAAlphabet{4}}}`.

# Arguments

- `input`: The input sequence as a `LongSequence{DNAAlphabet{4}}` or `String`.

# Keyword Arguments

- `alternative_start::Bool=false`: If set to `true`, the function considers alternative start codons when searching for ORFs. This increases the execution time by approximately 3x.
- `min_len::Int64=6`: The minimum length of the allowed ORF. By default, it allows ORFs that can encode at least one amino acid (e.g., `aa"M*"`).

"""
function getorfdna(
    sequence::LongSequence{DNAAlphabet{4}};
    alternative_start::Bool = false,
    min_len::Int64 = 6
)
    orfs = findorfs(sequence; alternative_start, min_len)
    seqs = Vector{LongSubSeq{DNAAlphabet{4}}}()
    @inbounds for i in orfs
        if i.strand == '+'
            push!(seqs, @view sequence[i.location])
        else
            newseq = reverse_complement(@view sequence[i.location])
            push!(seqs, newseq)
        end
    end
    return seqs
end

@testitem "getorfdna tests" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfseqs = getorfdna(seq01)

    @test length(orfseqs) == 5
    @test orfseqs[1] == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
end

"""
    getorfaa(input::LongSequence{DNAAlphabet{4}}; kwargs...)
    getorfaa(input::String; kwargs...) ## for strings per se

This function takes a `LongSequence{DNAAlphabet{4}}` or `String` sequence and identifies the open reading frames (ORFs) using the `findorfs()` function. The function then translates each ORF into an amino acid sequence and stores it in a `Vector{LongSubSeq{AminoAcidAlphabet}}`.

# Arguments

- `input`: The input sequence as a `LongSequence{DNAAlphabet{4}}` or `String`.

# Keyword Arguments

- `alternative_start::Bool=false`: If set to `true`, the function considers alternative start codons when searching for ORFs. This increases the execution time by approximately 3x.
- `min_len::Int64=6`: The minimum length of the allowed ORF. By default, it allows ORFs that can encode at least one amino acid (e.g., `aa"M*"`).

"""
function getorfaa(
    sequence::LongSequence{DNAAlphabet{4}};
    alternative_start::Bool = false,
    code::GeneticCode = BioSequences.standard_genetic_code,
    min_len::Int64 = 6
)
    orfs = findorfs(sequence; alternative_start, min_len)
    aas = Vector{LongSubSeq{AminoAcidAlphabet}}()
    @inbounds for i in orfs
        if i.strand == '+'
            push!(aas, translate(@view sequence[i.location]))
        else
            newaa = translate(reverse_complement(@view sequence[i.location]); code)
            push!(aas, newaa)
        end
    end
    return aas
end

@testitem "getorfaa tests" begin
    using BioSequences

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    aas = getorfaa(seq01)

    @test length(aas) == 5
    @test aas[1] == aa"MMHACMLVTS*"
end