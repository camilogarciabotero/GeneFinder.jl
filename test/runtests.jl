module GeneFinderTests

using Test
using BioSequences
using FASTX
using GeneFinder
using Aqua

include("findorfstest.jl")
include("aquatest.jl")

end
# @testset "GeneFinder.jl" begin
#     # Write your tests here.
# end

# using TestItemRunner: @run_package_tests
# @run_package_tests

# using Aqua
# using GeneFinder
# Aqua.test_all(GeneFinder)
