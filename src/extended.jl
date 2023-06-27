# Methods from main packages that expand their fuctions to this package structs
import BioSequences: translate, standard_genetic_code
import Base: show, length, iterate

function translate(
    ntseq::LongSubSeq{DNAAlphabet{4}};
    code::GeneticCode = standard_genetic_code,
    allow_ambiguous_codons = true,
    alternative_start = false,
)
    ntseq = copy(ntseq)
    translate(ntseq; code, allow_ambiguous_codons, alternative_start)
end

Base.sort(v::Vector{<:ORF}) = sort(v, by = _orf_sort_key)

function _orf_sort_key(orf::ORF)
    return (orf.location, orf.strand)
end

function Base.show(io::IO, tcm::TCM)
    nucleotides = sort(collect(keys(tcm.order)))

    # Print type
    println(io, "GeneFinder.TCM{$(typeof(tcm.order)), $(typeof(tcm.counts)):")

    # Print header
    max_digits = maximum([length(string(maximum(tcm.counts[:, i]))) for i in 1:size(tcm.counts, 2)])
    header_str = "   " * join([rpad(n, max_digits+1) for n in nucleotides], "")
    println(io, header_str)

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = tcm.counts[i, :]
        row_str = join([rpad(string(x), max_digits+1) for x in row], "")
        println(io, "$nucleotide1  $row_str")
    end
end

function Base.show(io::IO, tpm::TPM)
    nucleotides = sort(collect(keys(tpm.order)))

    # Print type
    println(io, "GeneFinder.tpm{$(typeof(tpm.order)), $(typeof(tpm.probabilities)):")

    # Print header
    max_digits = maximum([length(string(maximum(round.(tpm.probabilities[:, i], digits = 3)))) for i in 1:size(tpm.probabilities, 2)])
    header_str = "   " * join([rpad(n, max_digits+1) for n in nucleotides], "")
    println(io, header_str)

    # Print rows
    for (i, nucleotide1) in enumerate(nucleotides)
        row = round.(tpm.probabilities[i, :], digits = 3)
        row_str = join([rpad(string(x), max_digits+1) for x in row], "")
        println(io, "$nucleotide1  $row_str")
    end
end

function Base.show(io::IO, model::TransitionModel)
    # Print the type name
    println(io, "TransitionModel:")

    # Print the size of the transition probability matrix
    tpm_size = "  - Transition Probability Matrix (Size: $(size(model.tpm, 1)) × $(size(model.tpm, 2))):"
    println(io, tpm_size)
    for row in 1:size(model.tpm, 1)
        print(io, "    ")
        for col in 1:size(model.tpm, 2)
            print(io, round(model.tpm[row, col], digits=3))
            print(io, "\t")
        end
        println(io)
    end

    # Print the size of the initials matrix
    initials_size = "  - Initials (Size: $(size(model.initials, 1)) × $(size(model.initials, 2))):"
    println(io, initials_size)
    for row in 1:size(model.initials, 1)
        print(io, "    ")
        for col in 1:size(model.initials, 2)
            print(io, round(model.initials[row, col], digits=3))
            print(io, "\t")
        end
        println(io)
    end

    # Print the value of 'n'
    order = "  - order: $(model.n)"
    println(io, order)
end

Base.length(::TransitionModel) = 1

# function Base.iterate(model::TransitionModel, state=nothing)
#     if state === nothing
#         state = 1
#     else
#         state += 1
#     end

#     if state <= model.n
#         return (model.tpm[:, state], state), state
#     else
#         return nothing
#     end
# end

# function Base.iterate(model::TransitionModel, state=nothing)
#     if state === nothing
#         state = ones(Int, model.n)
#     else
#         state .= state .+ 1
#     end

#     valid_states = findall(state .<= model.n)
#     if !isempty(valid_states)
#         return TransitionModel((model.tpm[:, valid_states], state[valid_states]), state[valid_states])
#     else
#         return nothing
#     end
# end


# function Base.iterate(iter::LongNucOrViewIterator{4, T}, state=nothing) where T
#     if state === nothing
#         state = model.initials[nucleotideindexes[iter.iterable[1]]]
#     end

#     if iter.idx > 1
#         pair = LongDNA{4}([iter.iterable[iter.idx-1], iter.iterable[iter.idx]])
#         i, j = dinueclotideindexes[pair]
#         state *= model.tpm[i, j]
#     end

#     if iter.idx < length(iter.iterable)
#         iter.idx += 1
#         return (iter.iterable, state), state
#     else
#         return nothing
#     end
# end

# iterable = iterate(sequence)
# if iterable === nothing
#     return model.initials[nucleotideindexes[sequence]]
# else
#     _, probability = iterable
#     return probability
# end
# end

# function Base.iterate(model::TransitionModel, state=nothing)
#     if state === nothing
#         state = ones(Int, model.n)
#     else
#         state .= state .+ 1
#     end

#     valid_states = findall(state .<= model.n)
#     if !isempty(valid_states)
#         return model.tpm, model.initials
#     else
#         return nothing
#     end
# end
