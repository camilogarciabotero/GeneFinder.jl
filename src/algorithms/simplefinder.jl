using BioSequences
struct ORF
    start::Int64
    stop::Int64 
    strand::Char
end

const startcodon = dna"ATG"

function simplefinder(sequence::LongDNA)
    stopcodons = [dna"TAG", dna"TAA", dna"TGA"]
    orfs = Vector{ORF}()
    for strand in ['+', '-']
        if strand == '-'
            seq = reverse_complement(sequence)
        else
            seq = sequence
        end
        i = 1
        for i in 1:length(seq) - 3 
            if seq[i:i+2] == startcodon
                orf = nothing
                j = i
                while j < length(seq) - 3 && seq[j:j+2] ∉ stopcodons #biore"TAG|TGA|TAA"d # !isstop(dna[j:j+2])
                    orf = ORF(i,  j + 5, strand)
                    j+=3
                end
                if j < length(seq) - 3 && seq[j:j+2] ∈ stopcodons
                    push!(orfs, orf)
                end
            end
        end
    end
    return orfs
end

# seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# @time findorf(anotherseq)

# for i in findorf(seq)
#     println(seq[i.start:i.stop])
# end



# function isstop(dna::LongDNA)
#     return dna == dna"TAA" || dna == dna"TGA" || dna == dna"TAG"
# end

# string = "ATCA"

# ls = LongSequence{DNAAlphabet{2}}(string)