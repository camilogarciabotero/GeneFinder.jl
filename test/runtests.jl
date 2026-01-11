module GeneFinderTests

using Test
using BioSequences
using FASTX
using GeneFinder
using BioMarkovChains
using Aqua
using TestItems
using TestItemRunner

@run_package_tests

include("findorfstest.jl")
include("aquatest.jl")
# include("iotest.jl")

end

