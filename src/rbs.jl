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
function _rbswindows(orf::ORFI{N,F}) where {N,F}
    windowa = orf.strand == STRAND_POS ? (orf.first-16:orf.first-5) : (orf.last+5:orf.last+16)
    windowb = orf.strand == STRAND_POS ? (orf.first-10:orf.first-3) : (orf.last+3:orf.last+10)
    windowc = orf.strand == STRAND_POS ? (orf.first-18:orf.first-11) : (orf.last+11:orf.last+18)
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
    score = 0
    wsymb = (:A, :B, :C)
    rbsvect = Vector{RBS}()
    windows = _rbswindows(orf)

    for i in 1:length(windows)
        window = windows[i]
        symbol = wsymb[i]
        sqv = @view source(orf)[window]
        for (rbs, scr) in pairs(RBSMOTIFS)
            if occursin(rbs, sqv)
                push!(rbsvect, RBS(rbs.seq, window, scr, symbol))
                # score += scr  # Use the value associated with the key
            end
        end
    end
    return rbsvect
end

# export orf_rbs_score
# function orf_rbs_score(seq::SeqOrView{A}, orf::ORFI{N,F}) where {A,N,F}
#     # Initialize the score and the max scores dictionary
#     wsymb = (:a, :b, :c)
#     windows = _rbswindows(orf)
#     maxscores = Dict(:a => 0, :b => 0, :c => 0)

#     # Iterate over the windows and their corresponding symbols
#     @inbounds for i in 1:length(windows)
#         window = windows[i]
#         symbol = wsymb[i]
#         sqv = @view seq[window]

#         # Check for RBS motifs in the sequence view
#         @inbounds for (rbs, scr) in pairs(RBSMOTIFS)
#             if occursin(rbs, sqv)
#                 # Update the max score directly
#                 if scr > maxscores[symbol]
#                     maxscores[symbol] = scr
#                 end
#             end
#         end
#     end

#     # Sum the maximum scores of each window
#     total_score = sum(values(maxscores))
    
#     return total_score
# end

export _orf_rbs_score
function _orf_rbs_score(orf::ORFI{N,F}) where {N,F}
    # Initialize the score and the max scores dictionary
    wsymb = (:a, :b, :c)
    windows = _rbswindows(orf)
    maxscores = Dict(:a => 0, :b => 0, :c => 0)

    # Iterate over the windows and their corresponding symbols
    @inbounds for i in 1:length(windows)
        window = windows[i]
        symbol = wsymb[i]
        sqv = @view source(orf)[window]

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

export motifseq
function motifseq(rbs::RBS)
    return rbs.motif
end

export offset
function offset(rbs::RBS)
    return rbs.offset
end


# An idea of another implemetation of the a orf finder would be to use qstart = ExactSearchQuery(dna"ATG") 
# and then findall the start codons in the sequence and reverse sequence appended (seq * reverse_complement(seq))
# Use the memory layout of this findings to store also the frame, strand and view of the sequence. Figuring out the 
# how the strand and the location will be defined in the memory layout is the next step. 

## Another idea would be to fastly count the number of ATG in a sequence twice. This will be the approx number of ORFs
## in the sequence.Them we can use other way to store the ORFs in the memory layout.