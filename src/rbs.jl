export RBS, RBSMOTIFS

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

# ExactSearchQuery(dna"ATACG", iscompatible)  => 5, ## testing motif

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
#                  |----|:27

export _rbswindows
function _rbswindows(orf::ORFI{N,F}) where {N,F}
    windowa = orf.strand == STRAND_POS ? (orf.first-10:orf.first-3) : (orf.last+3:orf.last+10)
    windowb = orf.strand == STRAND_POS ? (orf.first-16:orf.first-5) : (orf.last+5:orf.last+16)
    windowc = orf.strand == STRAND_POS ? (orf.first-20:orf.first-11) : (orf.last+11:orf.last+20)
    windows = (windowa, windowb, windowc)
    return filter(window -> first(window) >= 1 && last(window) >= 1, windows)
end

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

export _findrbs
function _findrbs(orf::ORFI{N,F}) where {N,F}
    # s = reverse(seq)
    # score = 0
    wsymb = (:A, :B, :C)
    rbsvect = Vector{RBS}()
    windows = _rbswindows(orf)

    # Iterate over the windows and their corresponding symbols
    for i in 1:3#length(windows)
        window = windows[i]
        symbol = wsymb[i]
        wsqv = view(source(orf), window)
        ## here we can apply a logic to decid whether to use the reverse motif or the forward motif
        motifs = orf.strand == STRAND_POS ? FORWARDRBSMOTIFS : REVERSERBSMOTIFS
        for (rbs, scr) in pairs(motifs)
            motifranges::Vector{UnitRange{Int}} = findall(rbs, wsqv)
            if !isempty(motifranges)
                for motifrange in motifranges
                    offset = (first(window) + first(motifrange) - 1):(first(window) + last(motifrange) - 1)
                    # rbseq = orf.strand == STRAND_POS ? view(source(orf), offset) : reverse_complement(view(source(orf), offset))
                    rbseq = view(wsqv,motifrange)
                    push!(rbsvect, RBS(rbseq, offset, symbol, scr, orf.strand))
                end
            end
            # if occursin(rbs, sqv)
            #     push!(rbsvect, RBS(rbs.seq, window, scr, symbol))
            #     # score += scr  # Use the value associated with the key
            # end
        end
    end
    return rbsvect
end

export orf_rbs_score
function orf_rbs_score(orf::ORFI{N,F}) where {N,F}
    # Initialize the score and the max scores dictionary
    wsymb = (:a, :b, :c)
    windows = _rbswindows(orf)
    maxscores = Dict(:a => 0, :b => 0, :c => 0)

    # Iterate over the windows and their corresponding symbols
    @inbounds for i in 1:3
        window = windows[i]
        symbol = wsymb[i]
        wsqv = view(source(orf), window)

        motifs = orf.strand == STRAND_POS ? FORWARDRBSMOTIFS : REVERSERBSMOTIFS
        # Check for RBS motifs in the sequence view
        @inbounds for (rbs, scr) in pairs(motifs)
            if occursin(rbs, wsqv)
                # Update the max score directly
                if scr > maxscores[symbol]
                    maxscores[symbol] = scr
                end
            end
        end
    end

    # Sum the maximum scores of each window
    totscore = sum(values(maxscores))
    
    return totscore
end

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