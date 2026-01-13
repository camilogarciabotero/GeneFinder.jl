export NaiveFinder, NaiveFinderLazy

"""
    NaiveFinder <: GeneFinderMethod

A simple ORF finding method that detects all Open Reading Frames in a DNA sequence
using regular expression matching.

See also: [`NaiveFinderLazy`](@ref), [`GeneFinderMethod`](@ref)
"""
struct NaiveFinder <: GeneFinderMethod end

"""
    NaiveFinderLazy <: GeneFinderMethod

A memory-optimized variant of `NaiveFinder` that pre-allocates the ORF vector
based on estimated start codon counts.

See also: [`NaiveFinder`](@ref), [`GeneFinderMethod`](@ref)
"""
struct NaiveFinderLazy <: GeneFinderMethod end

"""
    _locationiterator(seq::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}

Create an iterator that yields ORF location ranges in a DNA sequence.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search.

# Keywords
- `alternative_start::Bool=false`: If `true`, matches NTG (ATG, GTG, TTG) as start codons;
  if `false`, only matches ATG.

# Returns
An iterator yielding `UnitRange{Int64}` objects representing ORF locations.

# Implementation Details
The function uses a regular expression to find ORFs:
- Pattern: `ATG(?:[N]{3})*?T(AG|AA|GA)` (or `NTG...` with alternative starts)
- `ATG` or `NTG`: Start codon
- `(?:[N]{3})*?`: Non-greedy match of any number of 3-nucleotide codons
- `T(AG|AA|GA)`: Stop codons (TAA, TAG, TGA)

The iterator uses `findfirst` with progressive offsets to find non-overlapping ORFs.

See also: [`NaiveFinder`](@ref)
"""
function _locationiterator(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false
) where {N}
    regorf = alternative_start ? biore"NTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    fdr(x) = let offset = first(x)
        r = findfirst(regorf, @view seq[offset+1:end])
        isnothing(r) ? nothing : (first(r) + offset):(last(r) + offset)
    end
    itr = takewhile(!isnothing, iterated(fdr, findfirst(regorf, seq)))
    return itr
end

@doc raw"""
    NaiveFinder(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> Vector{ORF{NaiveFinder}}

Find all Open Reading Frames (ORFs) in a DNA sequence using regular expression matching.

This method searches for ORFs on both the forward and reverse complement strands,
returning a sorted vector of `ORF{NaiveFinder}` objects.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.

# Keywords
- `alternative_start::Bool=false`: If `true`, uses extended start codons (ATG, GTG, TTG).
  Studies show that in *E. coli* K-12, these are used ~83%, ~14%, and ~3% respectively.
  Enabling this increases execution time by approximately 3x.
- `minlen::Int64=6`: Minimum ORF length in nucleotides. The default of 6 allows
  detection of the shortest possible ORF (`dna"ATGTGA"` → `aa"M*"`).

# Returns
- `Vector{ORF{NaiveFinder}}`: A sorted vector of ORFs found in the sequence.

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
orfs = NaiveFinder(seq)

# With alternative start codons and minimum length filter
orfs = NaiveFinder(seq; alternative_start=true, minlen=30)
```

# Scoring Note
This implementation does not include a scoring scheme by default. For coding potential
assessment, consider using a log-odds ratio score:

```math
S(x) = \sum_{i=1}^{L} \log \frac{a^{\mathscr{m}_{1}}_{x_{i-1} x_i}}{a^{\mathscr{m}_{2}}_{x_{i-1} x_i}}
```

Where the score compares the probability of the sequence under a coding model
versus a non-coding model. See [`lordr`](@ref) for more information.

See also: [`NaiveFinderLazy`](@ref), [`findorfs`](@ref), [`ORF`](@ref)
"""
function NaiveFinder(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    kwargs...
) where {N}
    seqlen = length(seq)
    orfs = Vector{ORF{NaiveFinder}}()
    
    # Handle the sequence name
    seqname = _varname(seq)
    if seqname === nothing
        seqname = :unnamedseq
    else
        seqname = Symbol(seqname)
    end

    @inbounds for strand in (PSTRAND, NSTRAND)
        s = strand == NSTRAND ? reverse_complement(seq) : seq
        @inbounds for location in @views _locationiterator(s; alternative_start)
            if length(location) >= minlen
                start = strand == PSTRAND ? location.start : seqlen - location.stop + 1
                stop = start + length(location) - 1
                frm = start % 3 == 0 ? 3 : start % 3

                push!(orfs, ORF{NaiveFinder}(seqname, start:stop, strand, Int8(frm), NamedTuple()))
            end
        end
    end
    return sort!(orfs)
end


### NAIVE FINDER LAZY IMPLEMENTATION ###

"""
    _estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}}) where {N} -> Int

Estimate the number of ORFs in a sequence for vector pre-allocation.

Counts ATG start codons (and their reverse complement CAT) using a k-mer iterator
to provide a heuristic estimate for the expected number of ORFs.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to analyze.

# Returns
- `Int`: Estimated ORF count (minimum of 10 to avoid zero allocation).

# Implementation
Uses `FwRvIterator` with 3-mers to efficiently count start codon occurrences
on both strands in a single pass.

See also: [`NaiveFinderLazy`](@ref)
"""
function _estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}}) where {N}
    # 4^3 = 64 possible 3-mers; only count ATG
    counts = 10 # Start with a base count to avoid zero allocation
    forwardtarget = Kmer{DNAAlphabet{4},3,1}(dna"ATG")
    reversetarget = Kmer{DNAAlphabet{4},3,1}(dna"GTA")  # Reverse of ATG

    for (fwdkmer, revkmer) in FwRvIterator{DNAAlphabet{4},3}(seq)
        if fwdkmer == forwardtarget || revkmer == reversetarget
            counts += 1
        end
    end

    return counts
end


"""
    _search_strand!(orfs, seq, seqname, strand, seqlen, alternative_start, minlen)

Search for ORFs on a single strand and append results to the ORF vector.

This is an internal helper function that avoids code duplication between
forward and reverse strand searching in `NaiveFinderLazy`.

# Arguments
- `orfs::Vector{ORF{NaiveFinderLazy}}`: Vector to append found ORFs to (mutated).
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search.
- `seqname::Symbol`: Identifier for the source sequence.
- `strand::Strand`: The strand being searched (`PSTRAND` or `NSTRAND`).
- `seqlen::Int`: Length of the original sequence (for coordinate transformation).
- `alternative_start::Bool`: Whether to use alternative start codons.
- `minlen::Int64`: Minimum ORF length filter.

# Coordinate Handling
- For `PSTRAND`: Coordinates are used directly.
- For `NSTRAND`: Coordinates are transformed from reverse complement positions
  back to original sequence positions.

See also: [`NaiveFinderLazy`](@ref)
"""
function _search_strand!(
    orfs::Vector{ORF{NaiveFinderLazy}},
    seq::NucleicSeqOrView{DNAAlphabet{N}},
    seqname::Symbol,
    strand::Strand,
    seqlen::Int,
    alternative_start::Bool,
    minlen::Int64
) where {N}
    @inbounds for location in _locationiterator(seq; alternative_start)
        if length(location) >= minlen
            # Adjust coordinates based on strand
            if strand === PSTRAND
                start, stop = location.start, location.stop
            else
                start = seqlen - location.stop + 1
                stop = seqlen - location.start + 1
            end
            
            frm = start % 3 == 0 ? 3 : start % 3
            push!(orfs, ORF{NaiveFinderLazy}(seqname, start:stop, strand, Int8(frm), NamedTuple()))
        end
    end
end

"""
    NaiveFinderLazy(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> Vector{ORF{NaiveFinderLazy}}

Memory-optimized ORF finder with smart pre-allocation.

Similar to `NaiveFinder`, but estimates the number of ORFs before searching
to pre-allocate the result vector, reducing memory allocations during the search.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.

# Keywords
- `alternative_start::Bool=false`: If `true`, uses extended start codons (ATG, GTG, TTG).
- `minlen::Int64=6`: Minimum ORF length in nucleotides.

# Returns
- `Vector{ORF{NaiveFinderLazy}}`: A sorted vector of ORFs found in the sequence.

# Performance
This variant is optimized for sequences where memory allocation overhead is significant.
It uses `_estimate_orf_count` to pre-allocate the result vector with `sizehint!`.

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
orfs = NaiveFinderLazy(seq)
```

See also: [`NaiveFinder`](@ref), [`findorfs`](@ref), [`ORF`](@ref)
"""
function NaiveFinderLazy(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    kwargs...
) where {N}
    seqlen = length(seq)
    
    # Estimate ORF count for smart pre-allocation
    estimated_orfs = _estimate_orf_count(seq)
    orfs = Vector{ORF{NaiveFinderLazy}}()
    sizehint!(orfs, estimated_orfs)
    
    # Handle the sequence name
    seqname = _varname(seq)
    seqname = seqname === nothing ? :unnamedseq : Symbol(seqname)

    # Search forward strand
    _search_strand!(orfs, seq, seqname, PSTRAND, seqlen, alternative_start, minlen)
    
    # Search reverse strand (compute reverse complement once)
    rev_seq = reverse_complement(seq)
    _search_strand!(orfs, rev_seq, seqname, NSTRAND, seqlen, alternative_start, minlen)
    
    sort!(orfs)
    return orfs
end

# function _estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}
#     atg_query = ExactSearchQuery(dna"ATG", iscompatible)
    
#     # Count ATG in forward and reverse strands
#     count_fwd = count(atg_query, seq)
#     count_rev = count(atg_query, reverse_complement(seq))
#     total_atg = count_fwd + count_rev
    
#     if alternative_start
#         # Add GTG and TTG counts
#         for motif in [dna"GTG", dna"TTG"]
#             query = ExactSearchQuery(motif, iscompatible)
#             total_atg += count(query, seq) + count(query, reverse_complement(seq))
#         end
#     end
    
#     # Heuristic: average ~1.5 ORFs per start codon
#     return max(100, div(total_atg * 3, 2))
# end
# const ATGQ1 = ApproximateSearchQuery(dna"ATG")
# const CATQ1 = ApproximateSearchQuery(dna"GTA")

# function _count_motif(q, seq::NucleicSeqOrView{DNAAlphabet{N}}) where {N}
#     n = 0
#     i = 1
#     while true
#         r = findnext(q, 0, seq, i)
#         r === nothing && return n
#         n += 1
#         i = first(r) + 1
#     end
# end

# function _estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}}) where {N}
#     return _count_motif(ATGQ1, seq) #+ _count_motif(GTAQ1, seq)
# end