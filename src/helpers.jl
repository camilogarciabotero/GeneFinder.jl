# function fasta_to_dna(input::String)
#     dnaseq = LongDNA{4}() # do something with record
#     FASTAReader(open(input)) do reader
#         for record in reader
#             seq = sequence(record)
#             dnaseq = LongDNA{4}(seq) # do something with record
#         end
#         return dnaseq
#     end
# end

# function fasta_to_dna(input::String)
#     dnaseq = Vector{LongDNA{4}}()
    
#     try
#         FASTAReader(open(input)) do reader
#             for record in reader
#                 seq = sequence(record)
#                 push!(dnaseq, LongDNA{4}(seq))
#             end
#         end
#     catch e
#         error("Error reading input file: $(e)")
#     end
    
#     return dnaseq
# end
Base.sort(v::Vector{<:ORF}) = sort(v, by = _orf_sort_key)

function _orf_sort_key(orf::ORF)
    return (orf.location, orf.strand)
end

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

# function locations(seq::LongDNA)
#     location = UnitRange[]
#     for i in locationiterator(seq)
#         push!(location, i)
#     end
#     return location
# end