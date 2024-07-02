using PrecompileTools: @setup_workload, @compile_workload

@setup_workload let
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    # using BioSequences, GeneFinder
    seq = randdnaseq(99)
    @compile_workload let
        _locationiterator(seq; alternative_start=false)
        findorfs(seq)
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        # findorfs(seq)
    end
end