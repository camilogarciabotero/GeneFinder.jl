export NaiveCollector

# abstract type GeneFinderMethod end
struct NaiveCollector <: GeneFinderMethod end

function NaiveCollector(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    scheme::Union{Nothing, Function} = nothing,
    overlap::Bool = false,
    kwargs...
) where {N}
    
    regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    framedict = Dict(0 => 3, 1 => 1, 2 => 2)
    revseq = reverse_complement(sequence)
    seqlen = length(sequence)
    seqname = get_var_name(sequence)

    function createorfs(x, strand)
        if length(x.captured[1]:x.captured[3]) < minlen
            return nothing
        end
        if strand == '+'
            start = x.captured[1]
            stop = x.captured[3] + 1
            frame = framedict[x.captured[1] % 3]
        else
            start = seqlen - x.captured[3] 
            stop = seqlen - x.captured[1] + 1
            frame = framedict[(seqlen - x.captured[3]) % 3]
        end
        seq = strand == '+' ? @view(sequence[start:stop]) : reverse_complement(@view(sequence[start:stop])) #@view(sequence[start:stop])
        scr = scheme === nothing ? 0.0 : scheme(seq; kwargs...)
        # fts = Dict(:score => score)
        fts = Features((score = scr,)) # rbs = RBS(biore"RRR"dna, 3:4, 1.0)
        return ORF{N,NaiveCollector}(seqname, start, stop, strand, frame, seq, fts, scheme)
    end

    orfs = Vector{ORF{N,NaiveCollector}}()
    for strand in ('+', '-')
        seq = strand == '-' ? revseq : sequence
        matches = eachmatch(regorf, seq, overlap)
        strandedorfs = filter(!isnothing, [createorfs(x, strand) for x in matches])
        append!(orfs, strandedorfs)
    end

    return sort!(orfs; kwargs...)
end
