export NaiveCollector

"""
    NaiveCollector <: GeneFinderMethod

A fast ORF finding method that uses `eachmatch` for efficient pattern matching.

By default finds non-overlapping ORFs, which is faster than `NaiveFinder`.

See also: [`NaiveFinder`](@ref), [`GeneFinderMethod`](@ref)
"""
struct NaiveCollector <: GeneFinderMethod end

"""
    NaiveCollector(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> ORFCollection{NaiveCollector}

Find Open Reading Frames (ORFs) in a DNA sequence using efficient pattern matching.

Returns an `ORFCollection` containing the source sequence and all found ORFs.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.

# Keywords
- `alternative_start::Bool=false`: Whether to consider alternative start codons (GTG, TTG).
- `minlen::Int64=6`: Minimum ORF length in nucleotides.
- `overlap::Bool=false`: Whether to allow overlapping ORFs.

# Returns
- `ORFCollection{NaiveCollector}`: Collection of ORFs bundled with the source sequence.

!!! note
    This method finds non-overlapping ORFs by default, making it faster than `NaiveFinder`.
    It uses `eachmatch` for efficient pattern matching.

!!! warning
    Using `overlap=true` increases runtime significantly, and some ORFs found may contain
    premature stop codons.

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = NaiveCollector(seq)

# With overlapping ORFs
collection = NaiveCollector(seq, overlap=true)
```

See also: [`NaiveFinder`](@ref), [`ORFCollection`](@ref)
"""
function NaiveCollector(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    overlap::Bool = false,
    kwargs...
)::ORFCollection{NaiveCollector} where {N}
    regorf::BioRegex = alternative_start ? biore"NTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    revseq = reverse_complement(seq)
    seqv = @view seq[begin:end]
    seqlen = length(seqv)
    
    function createorf(x, strand::Strand)
        orflen = length(x.captured[1]:x.captured[3])
        if orflen >= minlen
            if strand == PSTRAND
                start = x.captured[1]
                stop = x.captured[3] + 1
            else
                start = seqlen - x.captured[3]
                stop = seqlen - x.captured[1] + 1
            end
            frm = mod1(start, 3)
            return ORF{NaiveCollector}(start:stop, strand, Int8(frm))
        end
        return nothing
    end

    orfvec = Vector{ORF{NaiveCollector}}()
    
    for strand in (PSTRAND, NSTRAND)
        s = strand == NSTRAND ? @view(revseq[begin:end]) : @view(seq[begin:end])
        matches = eachmatch(regorf, s, overlap)
        for m in matches
            orf = createorf(m, strand)
            orf !== nothing && push!(orfvec, orf)
        end
    end

    sort!(orfvec)
    return ORFCollection(seqv, orfvec)
end