export RBS, RBSMOTIFS

struct RBS
    motif::SeqOrView # motif
    offset::UnitRange{Int64} # offset
    score::Int64
    window::Symbol

    function RBS(motif::SeqOrView{A}, offset::UnitRange{Int64}, score::Int64, windowsymbol::Symbol) where {A}
        return new(motif, offset, score, windowsymbol)
    end
    # rbsinst = RBS(biore"RRR"dna, 3:4, 1.0)
end

const RBSMOTIFS = Dict(
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
    # ExactSearchQuery(dna"ATACG", iscompatible)  => 5, ## testing motif
)

# export RBSMOTIFS, RBSQUERYMOTIGS
# RBSMOTIFS = (
#     dna"GGAGGA", dna"GGAGG", dna"GAGGA", dna"GGACGA", dna"GGATGA", dna"GGAAGA", dna"GGCGGA", dna"GGGGGA", dna"GGTGGA",
#     dna"GGAG", dna"GAGG", dna"AGGA", dna"GGTGG", dna"GGGGG", dna"GGCGG", dna"AGG", dna"GAG", dna"GGA", 
#     dna"GAAGA", dna"GATGA", dna"GACGA"
# )

# RBSQUERYMOTIFS = ExactSearchQuery.(RBSMOTIFS, iscompatible)


# RBSQUERYMOTIGS = ExactSearchQuery.(RBSMOTIFS, iscompatible)
    
# ExactSearchQuery(dna"ATACG", iscompatible)  => 5, ## testing motif

#       -20  -15  -10  -5    |-> start codon
#        |....|....|....|....ATG...
#EX1.              GGAGGACCCCATGACACACACAACAC  
#                  |----|:27

## Window checking a:
#           -16        -5    |-> start codon
#         ...|..........|....ATG...

## Window checking b:
#                 -10    -3  |-> start codon
#               ...|......|..ATG...

## Window checking c:
#         -18    -11         |-> start codon
#       ...|......|..........ATG...

export _rbswindows
function _rbswindows(seq::SeqOrView{A}, orf::ORFI{N,F}) where {A,N,F}

    windowa = orf.strand == STRAND_POS ? (orf.first-16:orf.first-5) : (orf.last+5:orf.last+16) 
    windowb = orf.strand == STRAND_POS ? (orf.first-10:orf.first-3) : (orf.last+3:orf.last+10)
    windowc = orf.strand == STRAND_POS ? (orf.first-18:orf.first-11) : (orf.last+11:orf.last+18)
                                                                    # This is not searching the correct RBS in the reverse complemente strand
    # windows = (a = windowa, b = windowb, c = windowc)
    windows = (windowa, windowb, windowc)
    return windows#filter(window -> first(window) >= 1, windows)
end

export findrbs
function findrbs(seq::SeqOrView{A}, orf::ORFI{N,F}) where {A,N,F}
    # s = reverse(seq)
    score = 0
    wsymb = (:a, :b, :c)
    rbsvect = Vector{RBS}()
    windows = _rbswindows(seq, orf)

    for i in 1:length(windows)
        window = windows[i]
        symbol = wsymb[i]
        sqv = @view seq[window]
        for (rbs, scr) in pairs(RBSMOTIFS)
            if occursin(rbs, sqv)
                push!(rbsvect, RBS(rbs.seq, window, scr, symbol))
                # score += scr  # Use the value associated with the key
            end
        end
    end
    return rbsvect
end

function orf_rbs_score(seq::SeqOrView{A}, orf::ORFI{N,F}) where {A,N,F}
    # Initialize the score and the max scores dictionary
    wsymb = (:a, :b, :c)
    windows = _rbswindows(seq, orf)
    maxscores = Dict(:a => 0, :b => 0, :c => 0)

    # Iterate over the windows and their corresponding symbols
    @inbounds for i in 1:length(windows)
        window = windows[i]
        symbol = wsymb[i]
        sqv = @view seq[window]

        # Check for RBS motifs in the sequence view
        @inbounds for (rbs, scr) in pairs(RBSMOTIFS)
            if occursin(rbs, sqv)
                # Update the max score directly
                if scr > maxscores[symbol]
                    maxscores[symbol] = scr
                end
            end
        end
    end

    # Sum the maximum scores of each window
    total_score = sum(values(maxscores))
    
    return total_score
end


### Another implementation


export RBSMOTIFS2
RBSMOTIFS2 = (
    RBS(dna"GGAGGA", 5:10, 27, :A),
    RBS(dna"GGAGGA", 3:8, 26, :B),
    RBS(dna"GGAGGA", 11:16, 25, :C),
    RBS(dna"GGAGG", 5:9, 24, :D),
    RBS(dna"GGAGG", 3:8, 23, :E),
    RBS(dna"GAGGA", 5:9, 22, :F),
    RBS(dna"GAGGA", 3:8, 21, :G),
    RBS(dna"GAGGA", 11:16, 20, :H),
    RBS(dna"GGACGA", 5:10, 19, :I),
    RBS(dna"GGATGA", 5:10, 19, :I),
    RBS(dna"GGAAGA", 5:10, 19, :I),
    RBS(dna"GGCGGA", 5:10, 19, :I),
    RBS(dna"GGGGGA", 5:10, 19, :I),
    RBS(dna"GGTGGA", 5:10, 19, :I),
    RBS(dna"GGAAGA", 3:9, 18, :J),
    RBS(dna"GGTGGA", 3:9, 18, :J),
    RBS(dna"GGAAGA", 11:17, 17, :K),
    RBS(dna"GGTGGA", 11:17, 17, :K),
    RBS(dna"GGAG", 5:9, 16, :L),
    RBS(dna"GAGG", 5:9, 16, :L),
    RBS(dna"AGGA", 5:9, 15, :M),
    RBS(dna"GGTGG", 5:9, 14, :N),
    RBS(dna"GGGGG", 5:9, 14, :N),
    RBS(dna"GGCNGG", 5:9, 14, :N),
    RBS(dna"AGG", 5:8, 13, :O),
    RBS(dna"GAG", 5:8, 13, :O),
    RBS(dna"GGA", 5:8, 13, :O),
    RBS(dna"AGGA", 11:15, 12, :P),
    RBS(dna"AGGA", 3:7, 11, :Q),
    RBS(dna"GAGGA", 13:19, 10, :R),
    RBS(dna"GAAGA", 5:10, 9, :S),
    RBS(dna"GATGA", 5:10, 9, :S),
    RBS(dna"GACGA", 5:10, 9, :S),
    RBS(dna"GGTGG", 3:8, 8, :T),
    RBS(dna"GGTGG", 11:16, 7, :U),
    RBS(dna"AGG", 11:14, 6, :V),
    RBS(dna"GAAGA", 3:8, 5, :W),
    RBS(dna"GAAGA", 11:16, 4, :X),
    RBS(dna"AGGA", 13:17, 3, :Y),
    RBS(dna"AGG", 13:16, 2, :Z),
    RBS(dna"GGAAGA", 13:19, 2, :Z),
    RBS(dna"GGTGG", 13:18, 2, :Z),
    RBS(dna"AGG", 3:6, 1, :Z)
)

export motifseq
function motifseq(rbs::RBS)
    return rbs.motif
end

export RBSMOTIFS2QUERY
RBSMOTIFS2QUERY = ExactSearchQuery.(motifseq.(RBSMOTIFS2))

export offset
function offset(rbs::RBS)
    return rbs.offset
end

function motifscore(seq::SeqOrView{A}) where {A}
    score = 0
    for rbs in RBSMOTIFS2
        if occursin(rbs.motif, seq)
            score += rbs.score
        end
    end
    return score
end


export orf_rbs_score2
function orf_rbs_score2(seq::SeqOrView{A}, orf::ORFI{N,F}) where {A,N,F}
    for i in 1:length(RBSMOTIFS2)
        rbs = RBSMOTIFS2[i]
        rbsq = RBSMOTIFS2QUERY[i]
        println(rbs)
        return findprev(rbsq, seq, orf.first)
    end
end


# An idea of another implemetation of the a orf finder would be to use qstart = ExactSearchQuery(dna"ATG") 
# and then findall the start codons in the sequence and reverse sequence appended (seq * reverse_complement(seq))
# Use the memory layout of this findings to store also the frame, strand and view of the sequence. Figuring out the 
# how the strand and the location will be defined in the memory layout is the next step. 

## Another idea would be to fastly count the number of ATG in a sequence twice. This will be the approx number of ORFs
## in the sequence.Them we can use other way to store the ORFs in the memory layout.