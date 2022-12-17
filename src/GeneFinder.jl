module GeneFinder

using BioSequences
using TestItems: @testitem


include("algorithms/simplefinder.jl")
export ORF, CDS, Proteins, simplefinder, findcds, findproteins

# include("findgenes.jl")



end
