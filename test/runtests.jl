module GeneFinderTests

using Test
using BioSequences
using FASTX
using GeneFinder

include("findorfstest.jl")

using Aqua: 
    test_ambiguities,
    test_persistent_tasks,
    test_piracies,
    test_stale_deps,
    test_unbound_args,
    test_undefined_exports

test_ambiguities(GeneFinder)
test_persistent_tasks(GeneFinder)
test_piracies(GeneFinder)
test_stale_deps(GeneFinder)
test_unbound_args(GeneFinder)
test_undefined_exports(GeneFinder)

end
# @testset "GeneFinder.jl" begin
#     # Write your tests here.
# end

# using TestItemRunner: @run_package_tests
# @run_package_tests

# using Aqua
# using GeneFinder
# Aqua.test_all(GeneFinder)
