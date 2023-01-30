module GeneFinder

using BioSequences
using FASTX
import BioSequences: GeneticCode
using IterTools
using TestItems: @testitem

include("types.jl")
export ORF, Codon, CDS, Protein

include("algorithms/findorfs.jl")
export locationiterator, findorfs

include("findgenes.jl")
export cdsgenerator, proteingenerator, get_cds, get_proteins

include("io.jl")
export write_cds, write_proteins, write_bed

end