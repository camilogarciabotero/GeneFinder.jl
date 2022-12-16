# # this will be the main functions taking all the algorithms
# using BioSequences

# include("algorithms/simplefinder.jl")



# """
# FindGene struct
# """
# mutable struct FindGene{S1,S2}
#     a::AlignedSequence{S1}
#     b::S2
# end

# function findgenes end

# function findgenes(::SimpleFinder, sequence::LongDNA) # type::GeneticCode
#     orfs = simplefinder(sequence)
#     seqs = Vector{CDS}()
    
#     for i in orfs
#         if i.strand == '-'
#             reversedsequence = reverse_complement(sequence)
#             seq = reversedsequence[i.location]
#         else
#             seq = sequence[i.location]
#         end
#         cds = CDS(i.location, i.strand, seq)
#         push!(seqs, cds)
#     end
#     return seqs
# end