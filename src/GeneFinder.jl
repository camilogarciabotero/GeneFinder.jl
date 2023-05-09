module GeneFinder

using BioSequences
import BioSequences: GeneticCode
using FASTX
using IterTools
using TestItems: @testitem

include("types.jl")
export ORF, Codon, CDS, Protein

include("algorithms/findorfs.jl")
export locationiterator, findorfs

include("findgenes.jl")
export cdsgenerator, proteingenerator, get_cds, get_proteins

include("io.jl")
export write_cds, write_proteins, write_bed, write_gff

include("helpers.jl")
export fasta_to_dna

end