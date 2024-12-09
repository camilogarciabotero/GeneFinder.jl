### Another implementation

# export RBSMOTIFS, RBSQUERYMOTIGS
# RBSMOTIFS = (
#     dna"GGAGGA", dna"GGAGG", dna"GAGGA", dna"GGACGA", dna"GGATGA", dna"GGAAGA", dna"GGCGGA", dna"GGGGGA", dna"GGTGGA",
#     dna"GGAG", dna"GAGG", dna"AGGA", dna"GGTGG", dna"GGGGG", dna"GGCGG", dna"AGG", dna"GAG", dna"GGA", 
#     dna"GAAGA", dna"GATGA", dna"GACGA"
# )

# RBSQUERYMOTIFS = ExactSearchQuery.(RBSMOTIFS, iscompatible)
# RBSQUERYMOTIGS = ExactSearchQuery.(RBSMOTIFS, iscompatible)
    

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

export _motifscore
function _motifscore(seq::SeqOrView{A}) where {A}
    score = 0
    for rbs in RBSMOTIFS2
        if occursin(rbs.motif, seq)
            score += rbs.score
        end
    end
    return score
end

export RBSMOTIFS2QUERY
RBSMOTIFS2QUERY = ExactSearchQuery.(motifseq.(RBSMOTIFS2))

export _findrbs2
function _findrbs2(orf::ORFI{N,F}) where {N,F}
    rbsvect = Vector{RBS}()#Vector{RBS}()
    for i in 1:length(RBSMOTIFS2)
        rbs = RBSMOTIFS2[i]
        rbsq = RBSMOTIFS2QUERY[i]
        match = findprev(rbsq, source(orf), i) #needs to deal with strand
        if match != nothing
         push!(rbsvect, RBS(rbs.motif, -(match, rbs.offset), rbs.score, rbs.window))
         println("Match found at: ", match)
        end
    end
    return rbsvect
end
