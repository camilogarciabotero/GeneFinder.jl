module GeneFinderTests

using Test
using BioSequences
using FASTX
# using GeneFinder
using Aqua
using TestItems
using TestItemRunner

@run_package_tests

include("findorfstest.jl")
# include("iotest.jl")
include("getindextest.jl")
# include("aquatest.jl")



@testitem "getorfs dna" default_imports=false begin

    using BioSequences: @dna_str, DNAAlphabet
    using GeneFinder: NaiveFinder, getorfs, findorfs, ORF
    using Test: @test

    seq01 = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
    orfseqs = getorfs(seq01, DNAAlphabet{4}(), NaiveFinder())

    # @test length(orfseqs) == 5
    # @test orfseqs[1] == dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"
end


end
