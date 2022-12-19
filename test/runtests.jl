# using Test

# @testset "GeneFinder.jl" begin
#     # Write your tests here.
# end

using TestItemRunner: @run_package_tests
@run_package_tests

# using Aqua
# using GeneFinder
# Aqua.test_all(GeneFinder)