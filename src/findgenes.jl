# this will be the main functions taking all the 
include("algorithms/simplefinder.jl")

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

function findcds(sequence::LongDNA)
    orfs = simplefinder(sequence)
    seqs = Vector{LongDNA}()
    for i in orfs
        if i.strand == '-'
            cds = reverse_complement(sequence[i.start:i.stop])
        else
            cds = sequence[i.start:i.stop]
        end
        push!(seqs, cds)
    end
    return seqs
end

findcds(seq)

#function findgenes(
#    sequence::LongDNA,
#    algorithm::PredictionAlgorithms,
#    type::GeneticCode)
