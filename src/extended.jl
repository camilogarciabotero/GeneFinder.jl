# Methods from main packages that expand their fuctions to this package structs
import BioSequences: translate, standard_genetic_code
import Base: show, length, iterate, sort

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