module GeneFinder

using BioSequences
using TestItems: @testitem

include("algorithms/simplefinder.jl")
export simplefinder, findcds, findproteins

include("types.jl")
export ORF, Codon, CDS, Protein, stopcodons

end