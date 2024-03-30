
export naivefinderscoring

# function _locationiterator(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}};
#     alternative_start::Bool = false
# ) where {N}
#     regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
#     # regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG([N]{3})*T(AG|AA|GA)?"dna # an attempt to make it non PCRE non-determinsitic
#     finder(x) = findfirst(regorf, sequence, first(x) + 1) # + 3
#     itr = takewhile(!isnothing, iterated(finder, findfirst(regorf, sequence)))
#     return itr
# end

function naivefinderscoring(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    seqlen = length(sequence)
    framedict = Dict(0 => 3, 1 => 1, 2 => 2)
    orfs = Vector{ORF}()

    for strand in ('+', '-')
        seq = strand == '-' ? reverse_complement(sequence) : sequence

        @inbounds for location in @views _locationiterator(seq; alternative_start)
            if length(location) >= min_len
                frame = strand == '+' ? framedict[location.start % 3] : framedict[(seqlen - location.stop + 1) % 3]
                start = strand == '+' ? location.start : seqlen - location.stop + 1
                stop = start + length(location) - 1
                push!(orfs, ORF(start:stop, strand, frame, dnaseqprobability(sequence[start:stop], ECOLICDS)))
            end
        end
    end
    return sort(orfs)
end