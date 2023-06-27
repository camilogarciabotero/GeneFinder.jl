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
``` 
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
    hasprematurestop(sequence::LongNucOrView{4})::Bool

Determine whether the `sequence` of type `LongSequence{DNAAlphabet{4}}` contains a premature stop codon.

Returns a boolean indicating whether the `sequence` has more than one stop codon.
"""
function hasprematurestop(sequence::LongNucOrView{4})::Bool
    
    stopcodons = [LongDNA{4}("TAA"), LongDNA{4}("TAG"), LongDNA{4}("TGA")]  # Create a set of stop codons
    
    length(sequence) % 3 == 0 || error("The sequence is not divisible by 3")
    
    occursin(biore"T(AG|AA|GA)"dna, sequence[end-2:end]) || error("There is no stop codon at the end of the sequence")

    @inbounds for i in 1:3:length(sequence) - 4
        codon = sequence[i:i+2]
        if codon in stopcodons
            return true
        end
    end

    return false
end

"""
    dinucleotides(sequence::LongNucOrView{4})

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
```
seq = dna"AGCTAGCTAGCT"

dinucleotides(seq)

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
function dinucleotides(sequence::LongNucOrView{4}; extended_alphabet::Bool = false)
    A = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    dinucleotides = vec([LongSequence{DNAAlphabet{4}}([n1, n2]) for n1 in A, n2 in A])

    # counts = zeros(Int64, length(dinucleotides))
    counts = Array{Int64,1}(undef, 64)
    for (index, pair) in enumerate(dinucleotides)
        count = 0
        for i in 1:length(sequence)-1
            if sequence[i] == pair[1] && sequence[i+1] == pair[2]
                count += 1
            end
        end
        counts[index] = count
    end

    pairsdict = Dict{LongSequence{DNAAlphabet{4}},Int64}()
    for (index, pair) in enumerate(dinucleotides)
        pairsdict[pair] = counts[index]
    end

    return pairsdict
end

function _int_to_dna(index::Int64; extended_alphabet::Bool = false)
    A = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    return LongSequence{DNAAlphabet{4}}([A[index]])
end

function _dna_to_int(nucleotide::DNA; extended_alphabet::Bool = false)
    A = extended_alphabet ? collect(alphabet(DNA)) : [DNA_A, DNA_C, DNA_G, DNA_T]
    return findfirst(nucleotide, LongSequence{DNAAlphabet{4}}(A))
end

# function codons(sequence::LongNucOrView{4})
#     A = [DNA_A, DNA_C, DNA_G, DNA_T]
#     trinucleotides = vec([LongSequence{DNAAlphabet{4}}([n1, n2, n3]) for n1 in A, n2 in A, n3 in A])

#     # counts = zeros(Int64, length(trinucleotides))
#     counts = Array{Int64,1}(undef, 64)
#     for (index, trio) in enumerate(trinucleotides)
#         count = 0
#         @inbounds for i in 1:length(sequence)-2
#             if sequence[i:i+2] == trio
#                 count += 1
#             end
#         end
#         counts[index] = count
#     end

#     codondict = Dict{LongSequence{DNAAlphabet{4}},Int64}()
#     for (index, codon) in enumerate(trinucleotides)
#         codondict[codon] = counts[index]
#     end

#     return codondict
# end

# """
#     eachcodon(sequence::LongSequence{DNAAlphabet{4}})

# Iterate through the codons in the `sequence` of type `LongSequence{DNAAlphabet{4}}`.

# Returns an iterator yielding `Codon` objects for each codon in the `sequence`.
# """
# function eachcodon(sequence::LongSequence{DNAAlphabet{4}})
#     seqbound = length(sequence) - 2
#     return (sequence[i:i+2] for i = 1:3:seqbound)
# end

# @testitem "eachcodon test" begin
#     using BioSequences

#     seq = dna"ATGGCGTA"

#     @test collect(eachcodon(seq)) == [Codon("ATG"), Codon("GCG")]
# end


# function hasprematurestop(sequence::LongSequence{DNAAlphabet{4}})::Bool
#     stop_codon_count = 0
#     @inbounds for codon in eachcodon(sequence)
#         if codon == dna"TAA" || codon == dna"TAG" || codon == dna"TGA"
#             println("Found stop codon: ", codon)
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