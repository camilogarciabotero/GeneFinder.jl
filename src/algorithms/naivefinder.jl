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
    NaiveFinder(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> ORFCollection{NaiveFinder}

Find all Open Reading Frames (ORFs) in a DNA sequence using regular expression matching.

Returns an `ORFCollection` containing the source sequence and all found ORFs,
providing a clean API for sequence extraction.

# Arguments
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.

# Keywords
- `alternative_start::Bool=false`: If `true`, uses extended start codons (ATG, GTG, TTG).
- `minlen::Int64=6`: Minimum ORF length in nucleotides.

# Returns
- `ORFCollection{NaiveFinder}`: Collection of ORFs bundled with the source sequence.

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = NaiveFinder(seq)

# Iterate over ORFs
for orf in collection
    println(sequence(collection, orf))
end

# Index access
first_seq = sequence(collection, 1)
```

See also: [`NaiveFinderLazy`](@ref), [`ORFCollection`](@ref), [`sequence`](@ref)
"""
function NaiveFinder(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    kwargs...
)::ORFCollection{NaiveFinder} where {N}
    seqlen = length(seq)
    orfvec = Vector{ORF{NaiveFinder}}()

    @inbounds for strand in (PSTRAND, NSTRAND)
        s = strand == NSTRAND ? reverse_complement(seq) : seq
        @inbounds for location in @views _locationiterator(s; alternative_start)
            if length(location) >= minlen
                start = strand == PSTRAND ? location.start : seqlen - location.stop + 1
                stop = start + length(location) - 1
                frm = mod1(start, 3)

                push!(orfvec, ORF{NaiveFinder}(start:stop, strand, Int8(frm)))
            end
        end
    end
    sort!(orfvec)
    return ORFCollection(convert(LongSubSeq{DNAAlphabet{N}}, seq), orfvec)
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
function _estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}})::Int where {N}
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
    _search_strand!(orfs, seq, strand, seqlen, alternative_start, minlen)

Search for ORFs on a single strand and append results to the ORF vector.

This is an internal helper function that avoids code duplication between
forward and reverse strand searching in `NaiveFinderLazy`.

# Arguments
- `orfs::Vector{ORF{NaiveFinderLazy}}`: Vector to append found ORFs to (mutated).
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search.
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
    strand::Strand,
    seqlen::Int,
    alternative_start::Bool,
    minlen::Int64
) where {N}
    @inbounds for location in _locationiterator(seq; alternative_start)
        if length(location) >= minlen
            if strand === PSTRAND
                start, stop = location.start, location.stop
            else
                start = seqlen - location.stop + 1
                stop = seqlen - location.start + 1
            end
            
            frm = mod1(start, 3)
            push!(orfs, ORF{NaiveFinderLazy}(start:stop, strand, Int8(frm)))
        end
    end
end

"""
    NaiveFinderLazy(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> ORFCollection{NaiveFinderLazy}

Memory-optimized ORF finder with smart pre-allocation.

# Returns
- `ORFCollection{NaiveFinderLazy}`: Collection of ORFs bundled with the source sequence.
"""
function NaiveFinderLazy(
    seq::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    kwargs...
)::ORFCollection{NaiveFinderLazy} where {N}
    seqlen = length(seq)
    
    # Estimate ORF count for smart pre-allocation
    estimated_orfs = _estimate_orf_count(seq)
    orfvec = Vector{ORF{NaiveFinderLazy}}()
    sizehint!(orfvec, estimated_orfs)

    # Search forward strand
    _search_strand!(orfvec, seq, PSTRAND, seqlen, alternative_start, minlen)
    
    # Search reverse strand (compute reverse complement once)
    rev_seq = reverse_complement(seq)
    _search_strand!(orfvec, rev_seq, NSTRAND, seqlen, alternative_start, minlen)
    
    sort!(orfvec)
    return ORFCollection(convert(LongSubSeq{DNAAlphabet{N}}, seq), orfvec)
end