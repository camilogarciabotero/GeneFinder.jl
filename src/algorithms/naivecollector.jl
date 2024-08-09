export NaiveCollector

# abstract type GeneFinderMethod end
struct NaiveCollector <: GeneFinderMethod end

"""
    NaiveCollector(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) -> Vector{ORFI{N,F}} where {N,F}

The `NaiveCollector` function searches for open reading frames (ORFs) in a DNA sequence. It takes the following arguments:

# Required Arguments

- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The nucleic sequence to search for ORFs.

# Keywords Arguments

- `alternative_start::Bool`: A flag indicating whether to consider alternative start codons. Default is `false`.
- `minlen::Int64`: The minimum length of an ORF. Default is `6`.
- `overlap::Bool`: A flag indicating whether to allow overlapping ORFs. Default is `false`.

The function returns a sorted vector of `ORFI{N,NaiveCollector}` objects, representing the identified ORFs.

!!! note
    This method finds, by default, non-overlapping ORFs in the given sequence. It is much faster than the `NaiveFinder` method.
     Althought it uses the same regular expression to find ORFs in a source sequence, 
     it levarages on the `eachmatch` function to find all the ORFs in the sequence.

!!! warning
    Using the `overlap = true` flag will increase the runtime of the function significantly, but some of the ORFs found may display
        premature stop codons.
"""
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
            else
                start = seqlen - x.captured[3] 
                stop = seqlen - x.captured[1] + 1
            end
            frm = start % 3 == 0 ? 3 : start % 3
            oseq = _orfseq(seq, start, stop, strand)
            fts = NamedTuple()
            return ORFI{N,NaiveCollector}(seqname, start, stop, strand, frm, oseq, fts)
        end
    end

    orfs = Vector{ORFI{N,NaiveCollector}}()
    for strand in (STRAND_POS, STRAND_NEG)
        s = strand == STRAND_NEG ? @view(revseq[begin:end]) : @view(seq[begin:end])
        matches = eachmatch(regorf, s, overlap)
        strandedorfs = filter(!isnothing, [createorfs(x, strand) for x in matches])
        # filter!(x -> !hasprematurestop(sequence(x)), strandedorfs)
        append!(orfs, strandedorfs)
        # println(strandedorfs)
    end

    return sort!(orfs)#orfs#sort!(orfs)
end