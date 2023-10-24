const NUCLEICINDEXES = Dict(DNA_A => 1, DNA_C => 2, DNA_G => 3, DNA_T => 4)

const DINUCLEICINDEXES = Dict(
    LongDNA{4}("AA") => [1, 1],
    LongDNA{4}("AC") => [1, 2],
    LongDNA{4}("AG") => [1, 3],
    LongDNA{4}("AT") => [1, 4],
    LongDNA{4}("CA") => [2, 1],
    LongDNA{4}("CC") => [2, 2],
    LongDNA{4}("CG") => [2, 3],
    LongDNA{4}("CT") => [2, 4],
    LongDNA{4}("GA") => [3, 1],
    LongDNA{4}("GC") => [3, 2],
    LongDNA{4}("GG") => [3, 3],
    LongDNA{4}("GT") => [3, 4],
    LongDNA{4}("TA") => [4, 1],
    LongDNA{4}("TC") => [4, 2],
    LongDNA{4}("TG") => [4, 3],
    LongDNA{4}("TT") => [4, 4],
)