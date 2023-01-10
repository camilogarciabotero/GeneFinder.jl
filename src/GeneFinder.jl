module GeneFinder

using BioSequences
import BioSequences: GeneticCode
using TestItems: @testitem

include("algorithms/simplefinder.jl")
export simplefind, simplefind_extended, simplecds_generator, simpleprot_generator

include("types.jl")
export ORF, Codon, CDS, Protein

include("helpers.jl")
export eachcodon, hasprematurestop

include("findgenes.jl")
export locationgenerator, locationgenerator_extended, orfgenerator, cdsgenerator, proteingenerator

include("io.jl")
export write_cds, write_proteins, write_bed

include("constants.jl")
export STOPCODONS, STARTCODON, EXTENDED_STARTCODONS

end