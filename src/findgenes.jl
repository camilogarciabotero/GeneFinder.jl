# this will be the main functions taking all the 
include("algorithms/simplefinder.jl")

# seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

function findcds(sequence::LongDNA)
    orfs = simplefinder(sequence)
    seqs = Vector{LongDNA}()
    for i in orfs
        if i.strand == '-'
            reversedsequence = reverse_complement(sequence)
            cds = reversedsequence[i.position]
        else
            cds = sequence[i.position]
        end
        push!(seqs, cds)
    end
    return seqs
end

function findproteins(sequence::LongDNA)
    cds = findcds(sequence)
    proteins = Vector{LongAA}()
    for i in cds
        protein = translate(i)
        push!(proteins, protein)
    end
    return proteins
end

#function findgenes(
#    sequence::LongDNA,
#    algorithm::PredictionAlgorithms,
#    type::GeneticCode)