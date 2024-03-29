module GeneFinderTests

using Test
using BioSequences
using FASTX
using GeneFinder
using Aqua

include("findorfstest.jl")
include("iotest.jl")
include("getindextest.jl")
include("aquatest.jl")

end
