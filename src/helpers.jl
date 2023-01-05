include("types.jl")

function eachcodon(sequence::LongDNA)
    seqbound = length(sequence) - 2
    return(Codon(sequence[i:i+2]) for i in 1:3:seqbound)
end

function hasprematurestop(sequence::LongDNA)::Bool
    stop_codon_count = 0
    @inbounds for codon in eachcodon(sequence)
        if codon ∈ stopcodons
            stop_codon_count += 1
        end
    end
    stop_codon_count > 1
end


# function hasprematurestop(sequence::LongDNA)::Bool
#     stop_codon_count = 0
#     seqbound = length(sequence)
#     @inbounds for i in 1:3:seqbound
#         if i+2 <= seqbound
#             codon = sequence[i:i+2]
#             if codon ∈ stopcodons
#                 stop_codon_count += 1
#             end
#         end
#     end
#     stop_codon_count > 1
# end

# function hasprematurestop02(sequence::LongDNA)
#     stop_codon_count = count(stopcodons, sequence[1:end-2])
#     stop_codon_count > 1
# end

# function _create_pairs(starts::Vector, stops::Vector)
#     combination = Vector{UnitRange}()
#     for i in starts
#         for j in stops
#             if i < j && length(i:j+2) % 3 == 0
#                 push!(combination, UnitRange(i,j))
#             end
#         end
#     end
#     return combination
# end

# This function create unique pairs while having all the starts fixed
# function _create_pairs(starts::Vector, stops::Vector)
#     combination = Dict{Int, UnitRange}()
#     for i in starts
#         for j in stops
#             if i < j && !haskey(combination, i) && length(i:j+2) % 3 == 0
#                 combination[i] = UnitRange(i,j)
#             end
#         end
#     end
#     return sort(collect(values(combination)), by=x -> x.start)
# end

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


# function _reversecomplement(sequence::LongDNA)
#     # create an empty string to hold the reverse complement
#     revcomp = LongDNA{4}()
#     # iterate over the characters in the dna sequence in reverse order
#     for base in reverse(sequence)
#         push!(revcomp, complement(base))
#     end
    
#     # return the reverse complement string
#     return revcomp

# end

