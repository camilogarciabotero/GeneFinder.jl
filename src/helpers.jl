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
    nucleotidetrans(sequence::LongNucOrView{4}, order::Int64=2; extended_alphabet::Bool=false)

Count the occurrences of nucleotide sequences of a given order in a LongDNA sequence.

# Arguments
- `sequence::LongNucOrView{4}`: The LongDNA sequence in which to count the nucleotide sequences.
- `order::Int64`: The order of the nucleotide sequences to count. Default is 2.
- `extended_alphabet::Bool`: Specify whether to use an extended alphabet that includes non-standard nucleotides. Default is `false`.

# Returns
A dictionary where the keys are the nucleotide sequences and the values are the counts.

"""
function nucleotidetrans(sequence::LongNucOrView{4}, order::Int64=2; extended_alphabet::Bool=false)
    A = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    nucleotides = [LongSequence{DNAAlphabet{4}}([nuc...]) for nuc in Iterators.product(fill(A, order)...)]

    nucleicdict = Dict{LongSequence{DNAAlphabet{4}}, Int64}()

    for nuc in nucleotides
        nucleicdict[nuc] = length(findall(ExactSearchQuery(nuc), sequence))
    end

    return nucleicdict
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
```
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
function transition_count_matrix(sequence::LongNucOrView{4}, order::Int64=2; extended_alphabet::Bool=false)

    transitions = nucleotidetrans(sequence, order; extended_alphabet)
    
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

@doc raw"""
    transition_probability_matrix(sequence::LongSequence{DNAAlphabet{4}})

Compute the transition probability matrix (TPM) of a given DNA sequence. Formally it construct `` \hat{A}`` where: 
```math
\[a_{ij} = P(X_t = j \mid X_{t-1} = i) = \frac{{P(X_{t-1} = i, X_t = j)}}{{P(X_{t-1} = i)}}\]

```

# Arguments
- `sequence::LongSequence{DNAAlphabet{4}}`: a `LongSequence{DNAAlphabet{4}}` object representing the DNA sequence.

# Keywords

- `extended_alphabet::Bool=false`: If true will pass the extended alphabet of DNA to search

# Returns
A `DTPM` object representing the transition probability matrix of the sequence.

# Example
```
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
function transition_probability_matrix(sequence::LongNucOrView{4}, order::Int64=2; extended_alphabet::Bool=false)
    dtcm = transition_count_matrix(sequence, order; extended_alphabet)
    rowsums = sum(dtcm.counts, dims = 2)
    freqs = round.(dtcm.counts ./ rowsums, digits = 3)

    freqs[isinf.(freqs)] .= 0.0
    freqs[isnan.(freqs)] .= 0.0

    return DTPM(dtcm.order, freqs)
end

@testitem "tpm" begin
    using BioSequences, GeneFinder
    seq = dna"CCTCCCGGACCCTGGGCTCGGGAC"
    tpm = transition_probability_matrix(seq)

    @test tpm.probabilities == [0.0 1.0 0.0 0.0; 0.0 0.5 0.2 0.3; 0.25 0.125 0.625 0.0; 0.0 0.667 0.333 0.0]
end

function initial_distribution(sequence::LongNucOrView{4}, order::Int64=2) ## π̂ estimates of the initial probabilies
    initials = Vector{Float64}()
    counts = transition_count_matrix(sequence, order).counts
    initials = sum(counts, dims = 1) ./ sum(counts)
    return initials
end

@doc raw"""
    sequenceprobability(sequence::LongNucOrView{4}, tpm::Matrix{Float64}, initials=Vector{Float64})

Compute the probability of a given sequence using a transition probability matrix and the initial probabilities distributions.

```math
P(X_1 = i_1, \ldots, X_T = i_T) = \pi_{i_1}^{T-1} \prod_{t=1}^{T-1} a_{i_t, i_{t+1}}
```

# Arguments
- `sequence::LongNucOrView{4}`: The input sequence of nucleotides.
- `tpm::Matrix{Float64}`: The transition probability matrix.
- `initials=Vector{Float64}`: Optional initial state probabilities. Default is an empty vector.

# Returns
- `probability::Float64`: The probability of the input sequence.

# Example
```

tpm = transition_probability_matrix(dna"CCTCCCGGACCCTGGGCTCGGGAC")

tpm
4×4 Matrix{Float64}:
0.0   1.0    0.0    0.0
0.0   0.5    0.2    0.3
0.25  0.125  0.625  0.0
0.0   0.667  0.333  0.0

initials = initial_distribution(dna"CCTCCCGGACCCTGGGCTCGGGAC")
initials
1×4 Matrix{Float64}:
0.0869565  0.434783  0.347826  0.130435

sequence = dna"CCTG"
sequence
4nt DNA Sequence:
CCTG

sequenceprobability(sequence, tpm, initials)
0.0217
"""
function sequenceprobability(sequence::LongNucOrView{4}, tpm::Matrix{Float64}, initials=Vector{Float64})
    
    nucleotideindexes = Dict(
        DNA_A => 1,
        DNA_C => 2,
        DNA_G => 3,
        DNA_T => 4
    )

    dinueclotideindexes = Dict(
        dna"AA" => [1,1],
        dna"AC" => [1,2],
        dna"AG" => [1,3],
        dna"AT" => [1,4],
        dna"CA" => [2,1],
        dna"CC" => [2,2],
        dna"CG" => [2,3],
        dna"CT" => [2,4],
        dna"GA" => [3,1],
        dna"GC" => [3,2],
        dna"GG" => [3,3],
        dna"GT" => [3,4],
        dna"TA" => [4,1],
        dna"TC" => [4,2],
        dna"TG" => [4,3],
        dna"TT" => [4,4]
    )

    init = initials[nucleotideindexes[sequence[1]]]

    probability = init
    
    for t in 1:length(sequence)-1
        
        pair = LongSequence{DNAAlphabet{4}}([sequence[t], sequence[t+1]])

        i = dinueclotideindexes[pair][1]
        j = dinueclotideindexes[pair][2]

        probability *= tpm[i,j]
    end
    return probability
end

function logodds(sequence::LongSequence{DNAAlphabet{4}}, model::TransitionModel, η::Float64)
    initcoding = model.codinginits
    initnoncoding = model.noncodinginits
    coding = model.coding
    noncoding = model.noncoding
    
    pcoding = sequenceprobability(sequence, coding, initcoding)
    pnoncoding = sequenceprobability(sequence, noncoding, initnoncoding)

    odds = log(pcoding/pnoncoding)

    if odds > η
        return "coding"
    else
        return "non-coding"
    end
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
