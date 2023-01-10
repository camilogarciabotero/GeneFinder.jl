module GeneFinder

using BioSequences
import BioSequences: GeneticCode
using TestItems: @testitem

include("types.jl")
export ORF, Codon, CDS, Protein

include("constants.jl")
export STOPCODONS, STARTCODON, EXTENDED_STARTCODONS

include("helpers.jl")
export eachcodon, hasprematurestop

include("algorithms/simplefinder.jl")
export orf_finder, cds_generator, protein_generator, get_cds, get_proteins

include("findgenes.jl")
export locationgenerator, orfgenerator, cdsgenerator, proteingenerator

include("io.jl")
export write_cds, write_proteins, write_bed

end