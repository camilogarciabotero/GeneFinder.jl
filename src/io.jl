export write_orfs_fna, write_orfs_faa, write_orfs_bed, write_orfs_gff

"""
    write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
    write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)

Write BED data to a file.

# Arguments
- `input`: The input DNA sequence NucSeq or a view.
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
- `finder`: The algorithm used to find ORFs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.

# Keywords
- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
- `minlen::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.
"""
function write_orfs_bed(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    @inbounds for orf in orfs
        println(output, orf.first, "\t", orf.last, "\t", orf.strand, "\t", orf.frame)
    end
end

function write_orfs_bed(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    open(output, "w") do f
        @inbounds for orf in orfs
            write(f, "$(orf.first)\t$(orf.last)\t$(orf.strand)\t$(orf.frame)\n")
        end
    end
end

"""
    write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
    write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)

Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.


# Arguments
- `input::NucleicAcidAlphabet{DNAAlphabet{N}}`: The input DNA sequence.
- `output::IO`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
- `finder::F`: The algorithm used to find ORFs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.

# Keywords

- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
- `minlen::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.

# Examples

```julia
filename = "output.fna"

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

open(filename, "w") do file
     write_orfs_fna(seq, file, NaiveFinder())
end
```
"""
function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false, 
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">", orf.groupname, " id=", id, " start=", orf.first, " stop=", orf.last, " strand=", orf.strand, " frame=", orf.frame, " score=", orf.features.fts.score)
        println(output, input[orf]) #
    end
end

function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">$(orf.groupname) id=$(id) start=$(orf.first) stop=$(orf.last) strand=$(orf.strand) frame=$(orf.frame) score=$(orf.features.fts.score)\n$(input[orf]))\n") # i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
        end
    end
end

"""
	write_orfs_faa(input::NucleicSeqOrView{DNAAlphabet{4}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
	write_orfs_faa(input::NucleicSeqOrView{DNAAlphabet{4}}, output::String, finder::F; kwargs...)

Write the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Arguments
- `input`: The input DNA sequence NucSeq or a view.
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
- `finder`: The algorithm used to find ORFs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `minlen::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.

# Examples 

```julia
filename = "output.faa"

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

open(filename, "w") do file
     write_orfs_faa(seq, file)
end
```
"""
function write_orfs_faa(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    code::GeneticCode = ncbi_trans_table[1]
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">", orf.groupname, " id=", id, " start=", orf.first, " stop=", orf.last, " strand=", orf.strand, " frame=", orf.frame, " score=", orf.features.fts.score)
        println(output, translate(input[orf]; code))
    end
end

function write_orfs_faa(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F}=NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    code::GeneticCode = ncbi_trans_table[1]
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">$(orf.groupname) id=$(id) start=$(orf.first) stop=$(orf.last) strand=$(orf.strand) frame=$(orf.frame) score=$(orf.features.fts.score)\n$(translate(input[orf]; code))\n")
        end
    end
end

"""
    write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
    write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)

Write GFF data to a file.

# Arguments
- `input`: The input DNA sequence NucSeq or a view.
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
- `finder`: The algorithm used to find ORFs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `minlen::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.

"""
function write_orfs_gff(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream,IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    println(output, "##gff-version 3\n##sequence-region Chr 1 $(length(input))")
    for (i, orf) in enumerate(orfs)
        id = string("ORF", lpad(string(i), padding, "0"))
        println(output, "Chr\t.\tORF\t", orf.first, "\t", orf.last, "\t.\t", orf.strand, "\t.\tID=", id, ";Name=", id, ";Frame=", orf.frame, ";Score=", orf.features.fts.score)
    end
end

function write_orfs_gff(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N, F<:GeneFinderMethod}
    # seq = fasta2bioseq(input)[1]
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(input))\n") 
        for (i, orf) in enumerate(orfs)
            id = string("ORF", lpad(string(i), padding, "0"))
            write(
                f,
                "Chr\t.\tORF\t$(orf.first)\t$(orf.last)\t.\t$(orf.strand)\t.\tID=$(id);Name=$(id);Frame=$(orf.frame);Score=$(orf.features.fts.score)\n"
            )
        end
    end
end


# Chr needs to be changed for fasta (?) or the name of seq (?) to fit well on IGV 