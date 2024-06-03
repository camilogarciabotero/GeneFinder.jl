#!/usr/bin/env julia

fasta = ARGS[1]

using BioSequences: NucleicSeqOrView, DNAAlphabet, LongDNA
using GeneFinder: findorfs, fasta2bioseq, lors, NaiveFinder
using BenchmarkTools: @btime

seq = fasta2bioseq(fasta)[1]
time = @btime findorfs(seq, NaiveFinder(), minlen=64, scheme=lors);

println("Time: $time")