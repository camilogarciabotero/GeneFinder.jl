module GeneFinder

using BioSequences
using TestItems: @testitem


include("algorithms/simplefinder.jl")
export ORF, simplefinder

include("findgenes.jl")
export CDS, Proteins, findcds, findproteins

# Write your package code here.

# Prediction types

# Gene types

# .

end
