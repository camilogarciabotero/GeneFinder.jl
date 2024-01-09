"""
    locationiterator(sequence::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}

This is an iterator function that uses regular expressions to search the entire ORF (instead of start and stop codons) in a `LongSequence{DNAAlphabet{4}}` sequence.
    It uses an anonymous function that will find the first regularly expressed ORF. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.

!!! note
    As a note of the implementation we want to expand on how the ORFs are found:

    The expression `(?:[N]{3})*?` serves as the boundary between the start and stop codons. Within this expression, the character class `[N]{3}` captures exactly three occurrences of any character (representing nucleotides using IUPAC codes). This portion functions as the regular codon matches. Since it is enclosed within a non-capturing group `(?:)` and followed by `*?`, it allows for the matching of intermediate codons, but with a preference for the smallest number of repetitions. 
    
    In summary, the regular expression `ATG(?:[N]{3})*?T(AG|AA|GA)` identifies patterns that start with "ATG," followed by any number of three-character codons (represented by "N" in the IUPAC code), and ends with a stop codon "TAG," "TAA," or "TGA." This pattern is commonly used to identify potential protein-coding regions within genetic sequences.

    See more about the discussion [here](https://discourse.julialang.org/t/how-to-improve-a-generator-to-be-more-memory-efficient-when-it-is-collected/92932/8?u=camilogarciabotero)

"""
function locationiterator(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false
) where {N}
    regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
    # regorf = alternative_start ? biore"DTG(?:[N]{3})*?T(AG|AA|GA)"dna : biore"ATG([N]{3})*T(AG|AA|GA)?"dna # an attempt to make it non PCRE non-determinsitic
    finder(x) = findfirst(regorf, sequence, first(x) + 1) # + 3
    itr = takewhile(!isnothing, iterated(finder, findfirst(regorf, sequence)))
    return itr
end

const FRAMEDICT = Dict(0 => 3, 1 => 1, 2 => 2) # sets the int conversion for frame assignment later

"""
    findorfs(sequence::NucleicAlphabet{DNAAlphabet{N}}; alternative_start::Bool=false, min_len::Int64=6)::Vector{ORF} where {N}

A simple implementation that finds ORFs in a DNA sequence.

The `findorfs` function takes a LongSequence{DNAAlphabet{4}} sequence and returns a Vector{ORF} containing the ORFs found in the sequence. 
    It searches entire regularly expressed CDS, adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.
        Extending the starting codons with the `alternative_start = true` will search for ATG, GTG, and TTG.
    Some studies have shown that in *E. coli* (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.
!!! note
    This function has not ORFs scoring scheme. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFs.
    
# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function findorfs(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false, 
    min_len::Int64 = 6
) where {N}
    orfs = Vector{ORF}(undef, 0)
    reversedseq = reverse_complement(sequence)
    seqlen = length(sequence)

    for strand in ('+', '-')
        seq = strand == '-' ? reversedseq : sequence

        @inbounds for location in @views locationiterator(seq; alternative_start)
            if length(location) >= min_len
                frame = strand == '+' ? FRAMEDICT[location.start % 3] : FRAMEDICT[(seqlen - location.stop + 1) % 3]
                push!(orfs, ORF(strand == '+' ? location : (seqlen - location.stop + 1):(seqlen - location.start + 1), strand, frame))
            end
        end
    end
    return sort(orfs)
end

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
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(sequence; alternative_start, min_len)
    seqs = Vector{LongSubSeq{DNAAlphabet{4}}}(undef, length(orfs))
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
    alternative_start::Bool = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(sequence; alternative_start, min_len)
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
function record_orfs_fna(
    sequence::NucleicSeqOrView{DNAAlphabet{N}}; 
    alternative_start::Bool = false, 
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(sequence; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    records = FASTARecord[]
    @inbounds for (index, orf) in enumerate(orfs)
        id = string(lpad(string(index), padding, "0"))
        header = "ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)"
        record = FASTARecord(header, sequence[orf])
        push!(records, record)
    end
    return records
end

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
function record_orfs_faa(
    sequence::NucleicSeqOrView{DNAAlphabet{N}};
    alternative_start::Bool = false, 
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(sequence; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    records = FASTARecord[]
    @inbounds for (index, orf) in enumerate(orfs)
        id = string(lpad(string(index), padding, "0"))
        header = "ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)"
        record = FASTARecord(header, translate(sequence[orf]; code))
        push!(records, record)
    end
    return records
end