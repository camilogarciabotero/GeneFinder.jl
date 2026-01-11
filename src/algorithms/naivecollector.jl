export NaiveCollector

struct NaiveCollector <: GeneFinderMethod end

"""
    NaiveCollector(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) -> Vector{ORF{F}} where {N,F}

The `NaiveCollector` function searches for open reading frames (ORFs) in a DNA sequence. It takes the following arguments:

# Required Arguments

- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The nucleic sequence to search for ORFs.

# Keywords Arguments

- `alternative_start::Bool`: A flag indicating whether to consider alternative start codons. Default is `false`.
- `minlen::Int64`: The minimum length of an ORF. Default is `6`.
- `overlap::Bool`: A flag indicating whether to allow overlapping ORFs. Default is `false`.

The function returns a sorted vector of `ORF{NaiveCollector}` objects, representing the identified ORFs.

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
    regorf::BioRegex = alternative_start ? biore"NTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    revseq = reverse_complement(seq)
    seqlen = length(seq)
    seqname = _varname(seq)
    
    if seqname === nothing
        seqname = :unnamedseq
    else
        seqname = Symbol(seqname)
    end
    
    function createorfs(x, strand::Strand)
        if length(x.captured[1]:x.captured[3]) > minlen
            if strand == PSTRAND
                start = x.captured[1]
                stop = x.captured[3] + 1
            else
                start = seqlen - x.captured[3]
                stop = seqlen - x.captured[1] + 1
            end
            frm = start % 3 == 0 ? 3 : start % 3
            fts = NamedTuple()
            return ORF{NaiveCollector}(seqname, start:stop, strand, Int8(frm), fts)
        end
        return nothing
    end

    orfs = Vector{ORF{NaiveCollector}}()
    for strand in (PSTRAND, NSTRAND)
        s = strand == NSTRAND ? @view(revseq[begin:end]) : @view(seq[begin:end])
        matches = eachmatch(regorf, s, overlap)
        strandedorfs = filter(!isnothing, [createorfs(x, strand) for x in matches])
        append!(orfs, strandedorfs)
    end

    return sort!(orfs)
end