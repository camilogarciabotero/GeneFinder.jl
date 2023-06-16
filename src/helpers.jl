# General purposes methods supporting main functions
"""
    fasta_to_dna(input::String)

Converts a FASTA formatted file (even if it is a multi-fasta) to an array of `LongDNA` objects.
"""
function fasta_to_dna(input::String)::Vector{LongDNA{4}}
    FASTAReader(open(input)) do reader
        return [LongDNA{4}(sequence(record)) for record in reader]
    end
end

"""
    eachcodon(sequence::LongDNA)

Iterate through the codons in the `sequence` of type `LongDNA`.

Returns an iterator yielding `Codon` objects for each codon in the `sequence`.
"""
function eachcodon(sequence::LongDNA)
    seqbound = length(sequence) - 2
    return(Codon(sequence[i:i+2]) for i in 1:3:seqbound)
end

# @testitem "eachcodon test" begin
#     using BioSequences

#     seq = dna"ATGGCGTA"

#     @test collect(eachcodon(seq)) == [Codon("ATG"), Codon("GCG")]
# end

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

# @testitem "hasprematurestop test" begin
#     using BioSequences

#     seq01 = dna"ATGGCGTA"
#     @test hasprematurestop(seq01) == false

#     seq02 = dna"ATGTCGTAATAA"
#     @test hasprematurestop(seq02) == true
# end

"""
    nucleotidefreqs(sequence::LongDNA) -> Dict{DNA, Float64}

Calculate the frequency of each nucleotide in a DNA sequence.

# Arguments
- `sequence::LongDNA`: A `LongDNA` sequence.

# Returns
A dictionary with each nucleotide in the sequence as a key, and its frequency as a value.

# Example
```julia
    
    seq = dna"CCTCCCGGACCCTGGGCTCGGGAC"

    nucleotidefreqs(seq)

    Dict{DNA, Float64} with 4 entries:
    DNA_T => 0.125
    DNA_A => 0.0833333
    DNA_G => 0.333333
    DNA_C => 0.458333
```
"""
function nucleotidefreqs(sequence::LongDNA)::Dict{DNA, Float64}
    counts = countmap(sequence)
    T = length(sequence)
    F = Dict(i => counts[i] / T for i in keys(counts))
    return F
end

function dinucleotidetrans(sequence::LongDNA)
    alphabet = unique(sequence)
    dinucleotides = vec([LongSequence{DNAAlphabet{4}}([n1, n2]) for n1 in alphabet, n2 in alphabet])
    
    pairsdict = Dict{LongDNA{4}, Int64}()
    for pair in dinucleotides
        instances = findall(ExactSearchQuery(pair), sequence)
        pairsdict[pair] = length(instances)
    end
    return pairsdict
end

# dinucleotidetrans01(seq)

# function dinucleotidetrans(sequence::LongDNA)
#     transitions = Dict{LongDNA{4}, Int}()
#     for i in 1:length(sequence)-1
#         dinucleotides = LongSequence{DNAAlphabet{4}}([sequence[i], sequence[i+1]])
#         transitions[dinucleotides] = get(transitions, dinucleotides, 0) + 1
#     end
#     return transitions
# end

function transition_count_matrix(sequence::LongDNA)
    alphabet = unique(sequence)
    dtcm = DTCM(alphabet)

    transitions = dinucleotidetrans(sequence)

    for (dinucleotide, count) in transitions
        nucleotide1 = dinucleotide[1]
        nucleotide2 = dinucleotide[2]
        i = dtcm.order[nucleotide1]
        j = dtcm.order[nucleotide2]
        dtcm.counts[i, j] = count
    end

    return dtcm
end

function transition_probability_matrix(sequence::LongDNA)
    dtcm = transition_count_matrix(sequence)
    rowsums = sum(dtcm.counts, dims=2)
    freqs = round.(dtcm.counts ./ rowsums, digits=3)
    # freqsform = [ println("%.2f", freqs[i,j]) for i in 1:size(freqs,1), j in 1:size(freqs,2) ]

    return DTPM(dtcm.order, freqs)
end



# function count_codons(seq::LongDNA)
#     codons = Vector{Codon}()
#     for i in 1:3:length(seq)
#         if i+2 <= length(seq)
#             codon = Codon(seq[i:i+2])
#             push!(codons, codon)
#         end
#     end
#     codons
# end