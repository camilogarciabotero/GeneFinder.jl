# https://github.com/althonos/pyrodigal/blob/main/pyrodigal/_pyrodigal.pxd

using BioSequences

abstract type ProdigalPrediction{T<:GenomicFeature} end

# Remember to make the implementation similar to BioAlignments so to have multiple algortithms to call when using 
# findgenes(string::DNAseq, algorithm::PredictionAlgorithms, type::GeneticCode)
# Make GeneticCode compatible for the type of prediction (?)

# bring gene structs

# type ProdigalPrediction


# Sequence mask

# Input seq

# Connection scorer

# Nodes

# Genes, an object for classifying and sorting predictions






