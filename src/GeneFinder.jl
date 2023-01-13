module GeneFinder

using BioSequences
import BioSequences: GeneticCode
using TestItems: @testitem

include("types.jl")
export ORF, Codon, CDS, Protein

include("helpers.jl")
export eachcodon, hasprematurestop

include("algorithms/simplefinder.jl")
export orf_finder, cds_generator, protein_generator, get_cds, get_proteins

include("findgenes.jl")
export locationgenerator, orfgenerator, cdsgenerator, proteingenerator

include("io.jl")
export write_cds, write_proteins, write_bed

const STOPCODONS = [Codon("TAG"), Codon("TAA"), Codon("TGA")]
const STARTCODON = ExactSearchQuery(Codon("ATG"), iscompatible)
const EXTENDED_STARTCODONS = ExactSearchQuery(dna"DTG", iscompatible)
# const EXTENDED_STARTCODONS = PWMSearchQuery([Codon("ATG"), Codon("GTG"), Codon("TTG")], 1.0)

end