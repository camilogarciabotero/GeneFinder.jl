# General purposes methods supporting main functions
"""
    fasta_to_dna(input::String)

Converts a FASTA formatted file (even if it is a multi-fasta) to an array of `LongSequence{DNAAlphabet{4}}` objects.
"""
function fasta_to_dna(input::String)::Vector{LongSequence{DNAAlphabet{4}}{4}}
    FASTAReader(open(input)) do reader
        return [LongSequence{DNAAlphabet{4}}(sequence(record)) for record in reader]
    end
end

"""
    nucleotidefreqs(sequence::LongSequence{DNAAlphabet{4}}) -> Dict{DNA, Float64}

Calculate the frequency of each nucleotide in a DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: A `LongSequence{DNAAlphabet{4}}` sequence.

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
function nucleotidefreqs(sequence::LongNucOrView{4})::Dict{DNA,Float64}
    counts = countmap(sequence)
    T = length(sequence)
    F = Dict(i => counts[i] / T for i in keys(counts))
    return F
end

"""
    dinucleotidetrans(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition counts of each dinucleotide in a given DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A dictionary with keys being `LongSequence{DNAAlphabet{4}}` objects representing
the dinucleotides, and values being the number of occurrences of each dinucleotide
in the sequence.

# Example
```julia
    seq = dna"AGCTAGCTAGCT"

    dinucleotidetrans(seq)

    Dict{LongSequence{DNAAlphabet{4}}, Int64} with 16 entries:
      GG => 0
      TC => 0
      GC => 3
      CG => 0
      CC => 0
      AG => 3
      TT => 0
      AC => 0
      TA => 2
      GT => 0
      GA => 0
      CT => 3
      CA => 0
      AT => 0
      AA => 0
      TG => 0
```
"""
function dinucleotidetrans(sequence::LongNucOrView{4}; extended_alphabet::Bool=false)
    shrtpalph = [DNA_A, DNA_C, DNA_G, DNA_T] # unique(sequence)
    extalph = collect(alphabet(DNA))
    alph = extended_alphabet ? extalph : shrtpalph
    dinucleotides = vec([LongSequence{DNAAlphabet{4}}([n1, n2]) for n1 in alph, n2 in alph])

    pairsdict = Dict{LongSequence{DNAAlphabet{4}},Int64}()
    for pair in dinucleotides
        instances = findall(ExactSearchQuery(pair), sequence)
        pairsdict[pair] = length(instances)
    end
    return pairsdict
end

"""
    transition_count_matrix(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition count matrix (TCM) of a given DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A `DTCM` object representing the transition count matrix of the sequence.

# Example
```julia
    seq = dna"AGCTAGCTAGCT"
    
    tcm = transition_count_matrix(seq)

    tcm.counts

    4×4 Matrix{Int64}:
     0  0  3  0
     0  0  0  3
     0  3  0  0
     2  0  0  0
 ```
 """
function transition_count_matrix(sequence::LongNucOrView{4}; extended_alphabet::Bool=false)

    transitions = dinucleotidetrans(sequence; extended_alphabet)
    
    alph = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]

    dtcm = DTCM(alph)

    for (dinucleotide, count) in transitions
        nucleotide1 = dinucleotide[1]
        nucleotide2 = dinucleotide[2]
        i = dtcm.order[nucleotide1]
        j = dtcm.order[nucleotide2]
        dtcm.counts[i, j] = count
    end

    return dtcm
end

"""
    transition_probability_matrix(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition probability matrix (TPM) of a given DNA sequence.

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A `DTPM` object representing the transition probability matrix of the sequence.

# Example
```julia
    seq = dna"AGCTAGCTAGCT"

    tpm = transition_probability_matrix(seq)

    tpm.probs

    4×4 Matrix{Float64}:
    0.0  0.0  1.0  0.0
    0.0  0.0  0.0  1.0
    0.0  1.0  0.0  0.0
    1.0  0.0  0.0  0.0
```
"""
function transition_probability_matrix(sequence::LongNucOrView{4}; extended_alphabet::Bool=false)
    dtcm = transition_count_matrix(sequence; extended_alphabet)
    rowsums = sum(dtcm.counts, dims = 2)
    freqs = round.(dtcm.counts ./ rowsums, digits = 3)

    freqs[isinf.(freqs)] .= 0.0
    freqs[isnan.(freqs)] .= 0.0

    return DTPM(dtcm.order, freqs)
end

@testitem "tpm" begin
    using BioSequences
    seq = dna"CCTCCCGGACCCTGGGCTCGGGAC"
    tpm = transition_probability_matrix(seq)

    @test tpm.probabilities == [0.0 1.0 0.0 0.0; 0.0 0.5 0.2 0.3; 0.25 0.125 0.625 0.0; 0.0 0.667 0.333 0.0]
end

function _int_to_dna(index; extended_alphabet::Bool=false)
    alph = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    return LongSequence{DNAAlphabet{4}}([alph[index]])
end

function generatednaseq(tpm::Matrix{Float64}, steps::Int64; extended_alphabet::Bool=false)
    newseq = LongSequence{DNAAlphabet{4}}()
    pf = transpose(tpm) # the Perron-Frobenius matrix
    trajectory = generate(pf, steps)
    for i in trajectory
        newseq = append!(newseq, _int_to_dna(i; extended_alphabet))
    end
    return newseq
end

# """
#     eachcodon(sequence::LongSequence{DNAAlphabet{4}})

# Iterate through the codons in the `sequence` of type `LongSequence{DNAAlphabet{4}}`.

# Returns an iterator yielding `Codon` objects for each codon in the `sequence`.
# """
# function eachcodon(sequence::LongSequence{DNAAlphabet{4}})
#     seqbound = length(sequence) - 2
#     return (Codon(sequence[i:i+2]) for i = 1:3:seqbound)
# end

# @testitem "eachcodon test" begin
#     using BioSequences

#     seq = dna"ATGGCGTA"

#     @test collect(eachcodon(seq)) == [Codon("ATG"), Codon("GCG")]
# end

# """
#     hasprematurestop(sequence::LongSequence{DNAAlphabet{4}})::Bool

# Determine whether the `sequence` of type `LongSequence{DNAAlphabet{4}}` contains a premature stop codon.

# Returns a boolean indicating whether the `sequence` has more than one stop codon.
# """
# function hasprematurestop(sequence::LongSequence{DNAAlphabet{4}})::Bool
#     stop_codon_count = 0
#     @inbounds for codon in eachcodon(sequence)
#         if codon ∈ STOPCODONS
#             stop_codon_count += 1
#         end
#     end
#     stop_codon_count > 1
# end

# @testitem "hasprematurestop test" begin
#     using BioSequences

#     seq01 = dna"ATGGCGTA"
#     @test hasprematurestop(seq01) == false

#     seq02 = dna"ATGTCGTAATAA"
#     @test hasprematurestop(seq02) == true
# end

# function count_codons(seq::LongSequence{DNAAlphabet{4}})
#     codons = Vector{Codon}()
#     for i in 1:3:length(seq)
#         if i+2 <= length(seq)
#             codon = Codon(seq[i:i+2])
#             push!(codons, codon)
#         end
#     end
#     codons
# end
