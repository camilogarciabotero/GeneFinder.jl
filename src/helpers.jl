function _create_pairs(starts::Vector, stops::Vector)
    combination = Vector{UnitRange}()
    for i in starts
        for j in stops
            if i < j && length(i:j+2) % 3 == 0
                push!(combination, UnitRange(i,j))
            end
        end
    end
    return combination
end


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


# function hasprematurestop(seq::LongDNA)
#     stop_codon_count = 0
#     for i in 1:3:length(seq)
#         if i+2 <= length(seq)
#             codon = seq[i:i+2]
#             if codon == dna"TAA" || codon == dna"TAG" || codon == dna"TGA"
#                 stop_codon_count += 1
#             end
#         end
#     end
#     stop_codon_count > 1
# end