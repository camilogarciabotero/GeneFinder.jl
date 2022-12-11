using BioSequences
include("../helpers.jl")

struct ORF
    start::Int64
    stop::Int64 
    strand::Char
end

function findorf(sequence::LongDNA)
    start_codon = ExactSearchQuery(dna"ATG", iscompatible)
    stop_codon = ExactSearchQuery(dna"TAR", iscompatible) # still misson TGA
    # stop_codon = ExactSearchQuery(dna"TRR", iscompatible)

    orfs = Vector{ORF}() # pre-allocate the vector with the appropriate size and then fill it in the loop. 

    seq = sequence
    for strand in ['+', '-']
        if strand == '-'
            seq = reverse_complement(sequence)
        end

        # starting_start_idx = map(x -> x.start, findall(start_codon, seq))
        starting_start_idx = [start_idx.start for start_idx in findall(start_codon, seq)]

        # stoping_stop_idx = map(x -> x.stop, findall(stop_codon, seq))
        stoping_stop_idx = [stop_idx.stop for stop_idx in findall(stop_codon, seq)]

        # combine all starts with the multiple stops creating an Array of all posible UnitRanges
        combinations = _create_pairs(starting_start_idx, stoping_stop_idx)

        for i in combinations
            if length(seq[i]) % 3 == 0 && length(findall(AA_Term, translate(seq[i]))) == 1
                # translation = translate(seq[i])
                orf = ORF(i.start, i.stop, strand) # seq[i]
                push!(orfs, orf)
            end
        end
    end
    return orfs
end

seq = dna"ATGCATGCATGCATGCTAGCTAGCTAGCTAGCTAGTAA"

# @time findorf(seq)

# for i in findorf(seq)
#     println(translate(i.cds))
# end


# start_idxs = findall(start_codon, seq)
# stop_idxs = findall(stop_codon, seq)
# starting_start_idx = collect(map(x -> x.start, start_idxs))
# stoping_stop_idx = collect(map(x -> x.stop, stop_idxs))

# struct CDS
#     orf::ORF
#     sequence::LongDNA
# end

function findcds(sequence::LongDNA)
    orfs = findorf(sequence)
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



# function findproteins(sequence::LongDNA)
    
# end