export RBS, orf_rbs_score, orbs, motifseq, offset
# public _rbswindows, _findrbs

"""
    RBS(motif::LongSubSeq, offset::UnitRange{Int64}, windowsymbol::Symbol, score::Int64, strand::Strand)

Represents a Ribosome Binding Site (RBS) in a DNA sequence.

# Fields
- `motif::LongSubSeq`: The nucleotide sequence of the RBS motif
- `offset::UnitRange{Int64}`: The position range of the RBS in the sequence
- `window::Symbol`: Symbol indicating the window type or region where the RBS was found
- `score::Int64`: Numerical score indicating the strength or confidence of the RBS prediction
- `strand::Strand`: Indicates whether the RBS is on the forward or reverse strand

# Constructor
    RBS(motif::LongSubSeq{A}, offset::UnitRange{Int64}, windowsymbol::Symbol, score::Int64, strand::Strand) where {A<:NucleicAcidAlphabet}

Creates a new RBS instance with the specified parameters.
"""
struct RBS
    motif::LongSubSeq # motif
    offset::UnitRange{Int64} # offset
    window::Symbol
    score::Int64
    strand::Strand

    function RBS(motif::LongSubSeq{A}, offset::UnitRange{Int64}, windowsymbol::Symbol, score::Int64, strand::Strand) where {A<:NucleicAcidAlphabet}
        return new(motif, offset, windowsymbol, score, strand)
    end
    # rbsinst = RBS(biore"RRR"dna, 3:4, 1.0)
end

## based on https://github.com/deprekate/PHANOTATE/blob/8b4e728171adc7d900ccbc59f9606bfd585e138c/phanotate_modules/functions.py#L48
const FORWARDRBSMOTIFS = Dict(
    ExactSearchQuery(dna"GGAGGA", iscompatible) => 27,
    ExactSearchQuery(dna"GGAGG", iscompatible)  => 24,
    ExactSearchQuery(dna"GAGGA", iscompatible)  => 22,
    ExactSearchQuery(dna"GGACGA", iscompatible) => 19,
    ExactSearchQuery(dna"GGATGA", iscompatible) => 19,
    ExactSearchQuery(dna"GGAAGA", iscompatible) => 19,
    ExactSearchQuery(dna"GGCGGA", iscompatible) => 19,
    ExactSearchQuery(dna"GGGGGA", iscompatible) => 19,
    ExactSearchQuery(dna"GGTGGA", iscompatible) => 19,
    ExactSearchQuery(dna"GGAG", iscompatible)   => 16,
    ExactSearchQuery(dna"GAGG", iscompatible)   => 16,
    ExactSearchQuery(dna"AGGA", iscompatible)   => 16,
    ExactSearchQuery(dna"GGTGG", iscompatible)  => 14,
    ExactSearchQuery(dna"GGGGG", iscompatible)  => 14,
    ExactSearchQuery(dna"GGCGA", iscompatible)  => 14,
    ExactSearchQuery(dna"AGG", iscompatible)    => 13,
    ExactSearchQuery(dna"GAG", iscompatible)    => 13,
    ExactSearchQuery(dna"GGA", iscompatible)    => 13,
    ExactSearchQuery(dna"GAAGA", iscompatible)  => 9,
    ExactSearchQuery(dna"GATGA", iscompatible)  => 9,
    ExactSearchQuery(dna"GACGA", iscompatible)  => 9,
)

const REVERSERBSMOTIFS = Dict(
    ExactSearchQuery(dna"TCCTCC", iscompatible) => 27,
    ExactSearchQuery(dna"CCTCC", iscompatible)  => 24,
    ExactSearchQuery(dna"TCCTC", iscompatible)  => 22,
    ExactSearchQuery(dna"TCGTCC", iscompatible) => 19,
    ExactSearchQuery(dna"TCATCC", iscompatible) => 19,
    ExactSearchQuery(dna"TTCCTC", iscompatible) => 19,
    ExactSearchQuery(dna"TCCGCC", iscompatible) => 19,
    ExactSearchQuery(dna"TCCCCT", iscompatible) => 19,
    ExactSearchQuery(dna"AGGTCC", iscompatible) => 19,
    ExactSearchQuery(dna"CTCC", iscompatible)   => 16,
    ExactSearchQuery(dna"CCTC", iscompatible)   => 16,
    ExactSearchQuery(dna"TCCT", iscompatible)   => 16,
    ExactSearchQuery(dna"GGTCC", iscompatible)  => 14,
    ExactSearchQuery(dna"CCCCC", iscompatible)  => 14,
    ExactSearchQuery(dna"AGCCG", iscompatible)  => 14,
    ExactSearchQuery(dna"CCT", iscompatible)    => 13,
    ExactSearchQuery(dna"CTC", iscompatible)    => 13,
    ExactSearchQuery(dna"TCC", iscompatible)    => 13,
    ExactSearchQuery(dna"TCTTC", iscompatible)  => 9,
    ExactSearchQuery(dna"TCATC", iscompatible)  => 9,
    ExactSearchQuery(dna"TCGTC", iscompatible)  => 9,
)

## based on the Prodigal paper:
# Spacer Range | Representative RBS Motif(s)
# ------------ | ---------------------------
# 3-4 bp       | GGA, GAG, AGG, AGxAG, AGGAG
# 5-10 bp      | GGAG, GAGG, GGxGG, AGGAGG
# 11-12 bp     | GGA, GAG, AGG, AGxAG, AGGAG
# 13-15 bp     | GGA, GAG, AGG, AGGA, AGGAG

## Window checking a:
#                 -10    -3  |-> start codon
#               ...|......|..ATG...

## Window checking b:
#           -16        -5    |-> start codon
#         ...|..........|....ATG...

## Window checking c:
#       -20      -11         |-> start codon
#       .|........|..........ATG...

## Example:
# ExactSearchQuery(dna"ATACG", iscompatible)  => 5, ## testing motif

#       -20  -15  -10  -5    |-> start codon
#        |....|....|....|....ATG...
#EX1.              GGAGGACCCCATGACACACACAACAC  
#                  |----|:RBS(dna"GGAGGA", 1:5, :A, 27, STRAND_POS)

"""
    _rbswindows(orf::ORF{F}) where {F} -> Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}

Generate potential ribosome binding site (RBS) windows for a given open reading frame (ORF).

Returns a tuple of UnitRanges representing three possible RBS windows relative to the ORF's start/end position,
filtering out any windows with invalid positions (< 1).

For positive strand ORFs:
- Window A: -10 to -3 positions upstream of start
- Window B: -16 to -5 positions upstream of start  
- Window C: -20 to -11 positions upstream of start

For negative strand ORFs:
- Window A: +3 to +10 positions downstream of end
- Window B: +5 to +16 positions downstream of end
- Window C: +11 to +20 positions downstream of end

# Arguments
- `orf::ORF{F}`: An ORF instance containing strand and position information

# Returns
- Tuple of valid UnitRanges representing RBS windows
"""
function _rbswindows(orf::ORF{F}) where {F}
    windowa = orf.strand == PSTRAND ? (leftposition(orf)-10:leftposition(orf)-3) : (rightposition(orf)+3:rightposition(orf)+10)
    windowb = orf.strand == PSTRAND ? (leftposition(orf)-16:leftposition(orf)-5) : (rightposition(orf)+5:rightposition(orf)+16)
    windowc = orf.strand == PSTRAND ? (leftposition(orf)-20:leftposition(orf)-11) : (rightposition(orf)+11:rightposition(orf)+20)
    windows = (windowa, windowb, windowc)
    return filter(window -> Base.first(window) >= 1 && Base.last(window) >= 1, windows)
end

"""
    _findrbs(orf::ORF{F}) where {F} -> Vector{RBS}

Search for Ribosome Binding Sites (RBS) in the given Open Reading Frame (ORF).

This function analyzes potential RBS locations in three windows upstream of the ORF start codon.
For each window, it searches for predefined RBS motifs (either forward or reverse depending on
the strand orientation) and records their locations and scores.

# Arguments
- `orf::ORF{F}`: An Open Reading Frame object containing sequence and strand information

# Returns
- `Vector{RBS}`: A vector of RBS objects, each containing:
  - The RBS sequence
  - The position range in the source sequence
  - A window symbol (`:A`, `:B`, or `:C`)
  - A score value
  - The strand orientation

# Implementation Details
- Searches in three windows upstream of the ORF
- Uses different motif sets for positive and negative strands
- Identifies all occurrences of RBS motifs in each window
- Records both the sequence and its position information

# Note
This is an internal and public function as indicated by the leading underscore.
"""
function _findrbs(orf::ORF{F}) where {F}
    wsymb = (:A, :B, :C)
    rbsvect = Vector{RBS}()
    windows = _rbswindows(orf)

    for i in 1:3
        window = windows[i]
        symbol = wsymb[i]
        wsqv = view(source(orf), window)
        motifs = strand(orf) == PSTRAND ? FORWARDRBSMOTIFS : REVERSERBSMOTIFS
        for (rbs, scr) in pairs(motifs)
            motifranges::Vector{UnitRange{Int}} = findall(rbs, wsqv)
            if !isempty(motifranges)
                for motifrange in motifranges
                    offset = (Base.first(window) + Base.first(motifrange) - 1):(Base.first(window) + Base.last(motifrange) - 1)
                    rbseq = view(wsqv, motifrange)
                    push!(rbsvect, RBS(rbseq, offset, symbol, scr, strand(orf)))
                end
            end
        end
    end
    return rbsvect
end

"""
    orf_rbs_score(orf::ORF{F}) where {F} -> Int64

Calculate a ribosome binding site (RBS) score for a given open reading frame (ORF).

This function evaluates potential ribosome binding sites in three windows upstream of
the ORF start codon by searching for known RBS motifs. It returns a total score based
on the strongest motif found in each window.

# Arguments
- `orf::ORF{F}`: An open reading frame interface object containing sequence and position information

# Returns
- `Float64`: The sum of the maximum RBS motif scores found in each window

# Details
The function:
1. Divides the upstream region into three windows (a, b, c)
2. Searches each window for matching RBS motifs
3. Keeps track of the highest scoring motif in each window
4. Returns the sum of the maximum scores across all windows

RBS motifs and their scores are predefined in `FORWARDRBSMOTIFS` for positive strand
and `REVERSERBSMOTIFS` for negative strand sequences.
"""
function orbs(orf::ORF{F}) where {F}
    wsymb = (:a, :b, :c)
    windows = _rbswindows(orf)
    maxscores = Dict(:a => 0, :b => 0, :c => 0)

    @inbounds for i in 1:3
        window = windows[i]
        symbol = wsymb[i]
        wsqv = view(source(orf), window)

        motifs = strand(orf) == PSTRAND ? FORWARDRBSMOTIFS : REVERSERBSMOTIFS
        @inbounds for (rbs, scr) in pairs(motifs)
            if occursin(rbs, wsqv)
                if scr > maxscores[symbol]
                    maxscores[symbol] = scr
                end
            end
        end
    end

    totscore = sum(values(maxscores))
    return totscore
end

const orf_rbs_score = orbs

export motifseq
function motifseq(rbs::RBS)
    return rbs.motif
end

export offset
function offset(rbs::RBS)
    return rbs.offset
end

# filter(x -> orf_rbs_score(x) > 9 && length(x) > 100, findorfs(phi))

# An idea of another implemetation of the a orf finder would be to use qstart = ExactSearchQuery(dna"ATG") 
# and then findall the start codons in the sequence and reverse sequence appended (seq * reverse_complement(seq))
# Use the memory layout of this findings to store also the frame, strand and view of the sequence. Figuring out the 
# how the strand and the location will be defined in the memory layout is the next step. 

## Another idea would be to fastly count the number of ATG in a sequence twice. This will be the approx number of ORFs
## in the sequence.Them we can use other way to store the ORFs in the memory layout.

## Another idea would be to use the findall function to find all the ATG in the sequence and then use the memory layout

# export findrbs
# function findrbs(seq::SeqOrView{A}, orf::ORFI{N,F}) where {A,N,F}
#     # s = reverse(seq)
#     score = 0
#     wsymb = (:a, :b, :c)
#     rbsvect = Vector{RBS}()
#     windows = _rbswindows(orf)

#     for i in 1:length(windows)
#         window = windows[i]
#         symbol = wsymb[i]
#         sqv = @view seq[window]
#         for (rbs, scr) in pairs(RBSMOTIFS)
#             if occursin(rbs, sqv)
#                 push!(rbsvect, RBS(rbs.seq, window, scr, symbol))
#                 # score += scr  # Use the value associated with the key
#             end
#         end
#     end
#     return rbsvect
# end