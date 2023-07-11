# General purposes methods supporting main functions
"""
    fasta_to_dna(input::String)

Converts a FASTA formatted file (even if it is a multi-fasta) to an array of `LongSequence{DNAAlphabet{4}}` objects.
"""
function fasta_to_dna(input::AbstractString)::Vector{LongSequence{DNAAlphabet{4}}}
    FASTAReader(open(input)) do reader
        return [LongSequence{DNAAlphabet{4}}(sequence(record)) for record in reader]
    end
end

# function gff_to_dna(input::AbstractString)
#     GFF3.Reader(open(input)) do reader
#         return [record for record in reader]
#     end
# end

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