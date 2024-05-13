export get_orfs_dna, get_orfs_aa
export getorfs

#### get_orfs_* methods ####

"""
    get_orfs_dna(sequence::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}

This function takes a `NucleicSeqOrView{DNAAlphabet{N}}` sequence and identifies the open reading frames (ORFs) using the `findorfs()` function. The function then extracts the DNA sequence of each ORF and stores it in a `Vector{LongSubSeq{DNAAlphabet{4}}}`.

# Arguments

- `sequence`: The input sequence as a `NucleicSeqOrView{DNAAlphabet{N}}`

# Keyword Arguments

- `alternative_start::Bool=false`: If set to `true`, the function considers alternative start codons when searching for ORFs. This increases the execution time by approximately 3x.
- `min_len::Int64=6`: The minimum length of the allowed ORF. By default, it allows ORFs that can encode at least one amino acid (e.g., `aa"M*"`).

"""
function get_orfs_dna(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    findermethod::Function = naivefinder,
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(sequence; findermethod, alternative_start, min_len)
    seqs = Vector{LongSubSeq{DNAAlphabet{N}}}(undef, length(orfs)) #Vector{}(undef, length(orfs)) # todo correct the output type
    @inbounds for (i, orf) in enumerate(orfs)
        seqs[i] = sequence[orf]
    end
    return seqs
end

"""
    get_orfs_aa(sequence::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N}

This function takes a `NucleicSeqOrView{DNAAlphabet{N}}` sequence and identifies the open reading frames (ORFs) using the `findorfs()` function. The function then translates each ORF into an amino acid sequence and stores it in a `Vector{LongSubSeq{AminoAcidAlphabet}}`.

# Arguments

- `sequence`: The input sequence as a `NucleicSeqOrView{DNAAlphabet{N}}`

# Keyword Arguments

- `alternative_start::Bool=false`: If set to `true`, the function considers alternative start codons when searching for ORFs. This increases the execution time by approximately 3x.
- `min_len::Int64=6`: The minimum length of the allowed ORF. By default, it allows ORFs that can encode at least one amino acid (e.g., `aa"M*"`).

"""
function get_orfs_aa(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    findermethod::Function = naivefinder,
    alternative_start::Bool = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(sequence; findermethod, alternative_start, min_len)
    aas = Vector{LongSubSeq{AminoAcidAlphabet}}(undef, length(orfs))
    @inbounds for (i, orf) in enumerate(orfs)
        aas[i] = translate(sequence[orf]; code)
    end
    return aas
end

#### record_orfs_* methods ####

"""
    record_orfs_fna(sequence::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N}

Record Open Reading Frames (ORFs) in a nucleic acid sequence.

# Arguments
- `sequence::NucleicSeqOrView{DNAAlphabet{N}}`: The nucleic acid sequence to search for ORFs.
- `alternative_start::Bool=false`: If `true`, consider alternative start codons for ORFs.
- `min_len::Int=6`: The minimum length of an ORF to be recorded.

# Returns
An array of `FASTARecord` objects representing the identified ORFs.

# Description
This function searches for Open Reading Frames (ORFs) in a given nucleic acid sequence. An ORF is a sequence of DNA that starts with a start codon and ends with a stop codon, without any other stop codons in between. By default, only the standard start codon (ATG) is considered, but if `alternative_start` is set to `true`, alternative start codons are also considered. The minimum length of an ORF to be recorded can be specified using the `min_len` argument.
"""
# function record_orfs_fna(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}}; 
#     alternative_start::Bool = false, 
#     min_len::Int64 = 6
# ) where {N}
#     orfs = findorfs(sequence; alternative_start, min_len)
#     norfs = length(orfs)
#     padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
#     records = FASTARecord[]
#     @inbounds for (index, orf) in enumerate(orfs)
#         id = string(lpad(string(index), padding, "0"))
#         header = "ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)"
#         record = FASTARecord(header, sequence[orf])
#         push!(records, record)
#     end
#     return records
# end

"""
    record_orfs_faa(sequence::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N}

This function takes a nucleic acid sequence and finds all open reading frames (ORFs) in the sequence. 
An ORF is a sequence of DNA that starts with a start codon and ends with a stop codon. 
The function returns a list of FASTA records, where each record represents an ORF in the sequence.

# Arguments
- `sequence`: The nucleic acid sequence to search for ORFs.
- `alternative_start`: A boolean indicating whether to consider alternative start codons. Default is `false`.
- `code`: The genetic code to use for translation. Default is the NCBI translation table 1.
- `min_len`: The minimum length of an ORF. ORFs shorter than this length will be ignored. Default is 6.

# Returns
- A list of FASTA records representing the ORFs found in the sequence.
"""
# function record_orfs_faa(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}};
#     alternative_start::Bool = false, 
#     code::GeneticCode = ncbi_trans_table[1],
#     min_len::Int64 = 6
# ) where {N}
#     orfs = findorfs(sequence; alternative_start, min_len)
#     norfs = length(orfs)
#     padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
#     records = FASTARecord[]
#     @inbounds for (index, orf) in enumerate(orfs)
#         id = string(lpad(string(index), padding, "0"))
#         header = "ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)"
#         record = FASTARecord(header, translate(sequence[orf]; code))
#         push!(records, record)
#     end
#     return records
# end

# function getorfs(
#     sequence::NucleicSeqOrView{DNAAlphabet{N}},
#     outseqtype::A,
#     method::M;
#     alternative_start::Bool = false,
#     min_len::Int64 = 6,
#     code::GeneticCode = ncbi_trans_table[1],
# ) where {N, A<:Alphabet, M<:GeneFinderMethod}
#     orfs = findorfs(sequence, method; alternative_start, min_len)
#     # seqs = Vector{LongSubSeq{DNAAlphabet{N}}}(undef, length(orfs)) #Vector{}(undef, length(orfs)) # todo correct the output type
    
#     if outseqtype == DNAAlphabet{N}()
#         seqs = Vector{LongSubSeq{DNAAlphabet{N}}}(undef, length(orfs)) #Vector{}(undef, length(orfs)) # todo correct the output type
#         @inbounds for (i, orf) in enumerate(orfs)
#             seqs[i] = sequence[orf]
#         end
#     elseif  outseqtype == AminoAcidAlphabet()
#         seqs = Vector{LongSubSeq{AminoAcidAlphabet}}(undef, length(orfs))
#         @inbounds for (i, orf) in enumerate(orfs)
#             seqs[i] = translate(sequence[orf]; code)
#         end
#     end

#     return seqs
# end


function getorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}},
    ::DNAAlphabet{N},
    method::M;
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N,M<:GeneFinderMethod}
    orfs = findorfs(sequence, method; alternative_start, min_len)
    seqs = Vector{LongSubSeq{DNAAlphabet{N}}}(undef, length(orfs)) #Vector{}(undef, length(orfs)) # todo correct the output type
    
    @inbounds for (i, orf) in enumerate(orfs)
        seqs[i] = sequence[orf]
    end

    return seqs
end


function getorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}},
    ::AminoAcidAlphabet,
    method::M;
    alternative_start::Bool = false,
    min_len::Int64 = 6,
    code::GeneticCode = ncbi_trans_table[1]
) where {N,M<:GeneFinderMethod}
    orfs = findorfs(sequence, method; alternative_start, min_len)
    seqs = Vector{LongSubSeq{AminoAcidAlphabet}}(undef, length(orfs))
    
    @inbounds for (i, orf) in enumerate(orfs)
        seqs[i] = translate(sequence[orf]; code)
    end

    return seqs
end