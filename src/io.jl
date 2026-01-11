export write_orfs_fna, write_orfs_faa, write_orfs_bed, write_orfs_gff, visualize, visualize_orfs

"""
    write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
    write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)

Write BED data to a file.

# Arguments
- `input`: The input DNA sequence NucSeq or a view.
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
- `finder`: The algorithm used to find ORFIs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.

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
    orfs = findorfs(input; finder, alternative_start, minlen)
    @inbounds for orf in orfs
        println(output, leftposition(orf), "\t", rightposition(orf), "\t", strand(orf), "\t", frame(orf))
    end
end

function write_orfs_bed(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    orfs = findorfs(input; finder, alternative_start, minlen)
    open(output, "w") do f
        @inbounds for orf in orfs
            write(f, "$(leftposition(orf))\t$(rightposition(orf))\t$(strand(orf))\t$(frame(orf))\n")
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
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        fts = isempty(features(orf)) ? "" : join(features(orf), ",")
        println(output, ">", seqid(orf), " id=", id, " start=", leftposition(orf), " stop=", rightposition(orf), " strand=", strand(orf), " frame=", frame(orf), " features=[", fts, "]")
        println(output, sequence(orf))
    end
end

function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            fts = isempty(features(orf)) ? "" : join(features(orf), ",")
            write(f, ">$(seqid(orf)) id=$(id) start=$(leftposition(orf)) stop=$(rightposition(orf)) strand=$(strand(orf)) frame=$(frame(orf)) features=[$(fts)]\n$(sequence(orf))\n")
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
- `finder`: The algorithm used to find ORFIs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `minlen::Int64=6`:  Length of the allowed ORFI. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFIs.

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
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        fts = isempty(features(orf)) ? "" : join(features(orf), ",")
        println(output, ">", seqid(orf), " id=", id, " start=", leftposition(orf), " stop=", rightposition(orf), " strand=", strand(orf), " frame=", frame(orf), " features=[", fts, "]")
        println(output, translate(sequence(orf); code))
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
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            fts = isempty(features(orf)) ? "" : join(features(orf), ",")
            write(f, ">$(seqid(orf)) id=$(id) start=$(leftposition(orf)) stop=$(rightposition(orf)) strand=$(strand(orf)) frame=$(frame(orf)) features=[$(fts)]\n$(translate(sequence(orf); code))\n")
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
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    println(output, "##gff-version 3\n##sequence-region Chr 1 $(length(input))")
    for (i, orf) in enumerate(orfs)
        id = string("ORFI", lpad(string(i), padding, "0"))
        fts = isempty(features(orf)) ? "" : join(features(orf), ",")
        println(output, "Chr\t.\tORFI\t", leftposition(orf), "\t", rightposition(orf), "\t.\t", strand(orf), "\t.\tID=", id, ";Name=", id, ";Frame=", frame(orf), ";Features=[", fts, "]")
    end
end

function write_orfs_gff(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N, F<:GeneFinderMethod}
    orfs = findorfs(input; finder, alternative_start, minlen)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(input))\n") 
        for (i, orf) in enumerate(orfs)
            id = string("ORFI", lpad(string(i), padding, "0"))
            fts = isempty(features(orf)) ? "" : join(features(orf), ",")
            write(f, "Chr\t.\tORFI\t$(leftposition(orf))\t$(rightposition(orf))\t.\t$(strand(orf))\t.\tID=$(id);Name=$(id);Frame=$(frame(orf));Features=[$(fts)]\n")
        end
    end
end

# Chr needs to be changed for fasta (?) or the name of seq (?) to fit well on IGV