"""
    eachcodon(sequence::LongDNA)

Iterate through the codons in the `sequence` of type `LongDNA`.

Returns an iterator yielding `Codon` objects for each codon in the `sequence`.
"""
function eachcodon(sequence::LongDNA)
    seqbound = length(sequence) - 2
    return(Codon(sequence[i:i+2]) for i in 1:3:seqbound)
end

@testitem "eachcodon test" begin
    using BioSequences

    seq = dna"ATGGCGTA"

    @test collect(eachcodon(seq)) == [Codon("ATG"), Codon("GCG")]
end

"""
    hasprematurestop(sequence::LongDNA)::Bool

Determine whether the `sequence` of type `LongDNA` contains a premature stop codon.

Returns a boolean indicating whether the `sequence` has more than one stop codon.
"""
function hasprematurestop(sequence::LongDNA)::Bool
    stop_codon_count = 0
    @inbounds for codon in eachcodon(sequence)
        if codon ∈ STOPCODONS
            stop_codon_count += 1
        end
    end
    stop_codon_count > 1
end

@testitem "hasprematurestop test" begin
    using BioSequences

    seq01 = dna"ATGGCGTA"
    @test hasprematurestop(seq01) == false

    seq02 = dna"ATGTCGTAATAA"
    @test hasprematurestop(seq02) == true
end