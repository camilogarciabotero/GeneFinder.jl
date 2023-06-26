"""
    locationiterator(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool=false)

This is an iterator function that uses regular expressions to search the entire CDS (instead of start and stop codons) in a `LongSequence{DNAAlphabet{4}}` sequence.
    It uses an anonymous function that will find the first regularly expressed CDS. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.
"""
function locationiterator(sequence::LongSequence{DNAAlphabet{4}}; alternative_start::Bool = false)
    regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
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

    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence

        @inbounds for location in locationiterator(seq; alternative_start)
            if length(location) >= min_len
                push!(orfs, ORF(location, strand))
            end
        end
    end
    return sort(orfs)
end

function findorfs(sequence::String; alternative_start::Bool = false, min_len::Int64 = 6)
    sequence = LongSequence{DNAAlphabet{4}}(sequence)
    orfs = Vector{ORF}()

    for strand in ['+', '-']
        seq = strand == '-' ? reverse_complement(sequence) : sequence

        @inbounds for location in locationiterator(seq; alternative_start)
            if length(location) >= min_len
                push!(orfs, ORF(location, strand))
            end
        end
    end
    return sort(orfs)
end

@testitem "findorfs test" begin
    using BioSequences

    # A random seq to start
    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfs01 = findorfs(seq01)

    @test findorfs(seq01) == [ORF(1:33, '+'), ORF(4:33, '+'), ORF(8:22, '+'), ORF(12:29, '+'), ORF(16:33, '+')]
    @test length(orfs01) == 5

    # > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
    seq02 = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
    orfs02 = findorfs(seq02)

    @test length(orfs02) == 12
    @test findorfs(seq02) == [
        ORF(29:40, '+'),
        ORF(137:145, '+'),
        ORF(164:184, '+'),
        ORF(173:184, '+'),
        ORF(236:241, '+'),
        ORF(248:268, '+'),
        ORF(362:373, '+'),
        ORF(470:496, '+'),
        ORF(551:574, '+'),
        ORF(569:574, '+'),
        ORF(581:601, '+'),
        ORF(695:706, '+'),
    ]
end
