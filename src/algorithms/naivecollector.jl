export NaiveCollector

# abstract type GeneFinderMethod end
struct NaiveCollector <: GeneFinderMethod end

function NaiveCollector(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    overlap::Bool = false,
    kwargs...
) where {N}
    regorf = alternative_start ? biore"NTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    revseq = reverse_complement(seq)
    seqlen = length(seq)
    seqname = _varname(seq)
    
    function createorfs(x, strand)
        while length(x.captured[1]:x.captured[3]) > minlen
            if strand == STRAND_POS
                start = x.captured[1]
                stop = x.captured[3] + 1
                # frame = framedict[x.captured[1] % 3]
            else
                start = seqlen - x.captured[3] 
                stop = seqlen - x.captured[1] + 1
                # frame = framedict[(seqlen - x.captured[3]) % 3]
            end
            frm = start % 3 == 0 ? 3 : start % 3
            oseq = _orfseq(seq, start, stop, strand)
            fts = NamedTuple()
            return ORF{N,NaiveCollector}(seqname, start, stop, strand, frm, oseq, fts)
        end
    end

    orfs = Vector{ORF{N,NaiveCollector}}()
    for strand in (STRAND_POS, STRAND_NEG)
        s = strand == STRAND_NEG ? revseq : seq
        matches = eachmatch(regorf, s, overlap)
        strandedorfs = filter(!isnothing, [createorfs(x, strand) for x in matches])
        # filter!(x -> !hasprematurestop(sequence(x)), strandedorfs)
        append!(orfs, strandedorfs)
        # println(strandedorfs)
    end

    return sort!(orfs)#orfs#sort!(orfs)
end