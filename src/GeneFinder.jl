module GeneFinder

using BioSequences  #standard_genetic_code
using TestItems: @testitem

include("algorithms/simplefinder.jl")
export simplefind, simplefind_extended, simplecds_generator, simpleprot_generator

include("types.jl")
export ORF, Codon, stopcodons, startcodon, extended_startcodons

include("helpers.jl")
export eachcodon, hasprematurestop

include("findgenes.jl")
export locationgenerator, locationgenerator_extended, orfgenerator, cdsgenerator, proteingenerator

end