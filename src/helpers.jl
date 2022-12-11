function _create_pairs(starts::Vector, stops::Vector)
    combination = Vector{UnitRange}()
    for i in starts
        for j in stops
            if i < j
                push!(combination, UnitRange(i,j))
            end
        end
    end
    return combination
end

## This function create unique pairs while having all the starts fixed
# function _create_pairs(starts::Vector, stops::Vector)
#     combination = Dict{Int, UnitRange}()
#     for i in starts
#         for j in stops
#             if i < j && !haskey(combination, i)
#                 combination[i] = UnitRange(i,j)
#             end
#         end
#     end
#     return sort(collect(values(combination)), by=x -> x.start)
# end

# starts = [1, 3, 5, 7]
# stops = [2, 4, 6, 8]

# _create_pairs(starts,stops)

