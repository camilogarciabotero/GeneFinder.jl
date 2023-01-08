module GeneFinder

using BioSequences
using TestItems: @testitem

include("algorithms/simplefinder.jl")
export simplefinder, simplefinderextended, cdsgenerator, proteingenerator

include("types.jl")
export ORF, Codon, stopcodons

include("helpers.jl")
export eachcodon, hasprematurestop

include("findgenes.jl")
export locationgenerator, orfgenerator

end