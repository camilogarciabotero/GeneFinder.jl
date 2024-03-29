export write_orfs_fna, write_orfs_faa, write_orfs_bed, write_orfs_gff

"""
    write_orfs_bed(input::NucleicAcidAlphabet{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}
    write_orfs_bed(input::NucleicAcidAlphabet{DNAAlphabet{N}}, output::String; kwargs...) where {N}
    write_orfs_bed(input::String, output::Union{IOStream, IOBuffer}; kwargs...)
    write_orfs_bed(input::String, output::String; kwargs...)

Write BED data to a file.

# Arguments
- `input::NucleicAcidAlphabet{DNAAlphabet{N}}`: The input DNA sequence.
- `output::String`: The output file path.
- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
- `min_len::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.
"""
function write_orfs_bed(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer}; 
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    @inbounds for orf in orfs
        println(output, orf.location.start, "\t", orf.location.stop, "\t", orf.strand, "\t", orf.frame)
    end
end

function write_orfs_bed(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    open(output, "w") do f
        @inbounds for orf in orfs
            write(f, "$(orf.location.start)\t$(orf.location.stop)\t$(orf.strand)\t$(orf.frame)\n")
        end
    end
end

function write_orfs_bed(
    input::String,
    output::Union{IOStream, IOBuffer}; 
    alternative_start::Bool = false,
    min_len::Int64 = 6
)
    input = fasta_to_dna(input)[1] # rewrite input to be a DNA sequence
    orfs = findorfs(input; alternative_start, min_len)
    @inbounds for orf in orfs
        println(output, orf.location.start, "\t", orf.location.stop, "\t", orf.strand, "\t", orf.frame)
    end
end

function write_orfs_bed(
    input::String,
    output::String;
    alternative_start::Bool = false,
    min_len::Int64 = 6
)
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    open(output, "w") do f
        @inbounds for orf in orfs
            write(f, "$(orf.location.start)\t$(orf.location.stop)\t$(orf.strand)\t$(orf.frame)\n")
        end
    end
end

"""
    write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}
    write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}
    write_orfs_fna(input::String, output::Union{IOStream, IOBuffer}; kwargs...)
    write_orfs_fna(input::String, output::String; kwargs...)

Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
- `min_len::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.
# Examples

```julia
filename = "output.fna"

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

open(filename, "w") do file
     write_orfs_fna(seq, file)
end
```
"""
function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    alternative_start::Bool = false, 
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF", id, " id=", id, " start=", orf.location.start, " stop=", orf.location.stop, " strand=", orf.strand, " frame=", orf.frame)
        println(output, input[orf])
    end
end

function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}}, 
    output::String; 
    alternative_start::Bool = false,
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(input[orf])\n") # i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
        end
    end
end

function write_orfs_fna(
    input::String,
    output::Union{IOStream, IOBuffer};
    alternative_start::Bool = false, 
    min_len::Int64 = 6
)
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF", id, " id=", id, " start=", orf.location.start, " stop=", orf.location.stop, " strand=", orf.strand, " frame=", orf.frame)
        println(output, input[orf])
    end
end

function write_orfs_fna(
    input::String, 
    output::String; 
    alternative_start::Bool = false,
    min_len::Int64 = 6
) 
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(input[orf])\n") # i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
        end
    end
end

"""
	write_orfs_faa(input::LongSequence{DNAAlphabet{4}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}
	write_orfs_faa(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...) where {N}
    write_orfs_faa(input::String, output::Union{IOStream, IOBuffer}; kwargs...)
    write_orfs_faa(input::String, output::String; kwargs...)

Write the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
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
    alternative_start::Bool = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF", id, " id=", id, " start=", orf.location.start, " stop=", orf.location.stop, " strand=", orf.strand, " frame=", orf.frame)
        println(output, translate(input[orf]; code))
    end
end

function write_orfs_faa(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    alternative_start::Bool = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6,
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(translate(input[orf]; code))\n")
        end
    end
end

function write_orfs_faa(
    input::String,
    output::Union{IOStream, IOBuffer};
    alternative_start::Bool = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6
) 
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF", id, " id=", id, " start=", orf.location.start, " stop=", orf.location.stop, " strand=", orf.strand, " frame=", orf.frame)
        println(output, translate(input[orf]; code))
    end
end

function write_orfs_faa(
    input::String,
    output::String;
    alternative_start::Bool = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len::Int64 = 6,
)
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(translate(input[orf]; code))\n")
        end
    end
end

"""
    write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}
    write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}
    write_orfs_gff(input::String, output::Union{IOStream, IOBuffer}; kwargs...)
    write_orfs_gff(input::String, output::String; kwargs...)

Write GFF data to a file.

# Arguments
- `input`: The input DNA sequence.
- `output`: The output file to write the GFF data to.
- `alternative_start::Bool=false`: If true, extended start codons will be considered during the search, increasing the execution time.
- `min_len::Int64=6`: The minimum length of the allowed ORF. The default value allows for possible encoding proteins with the `aa"M*"` sequence.

"""
function write_orfs_gff(
    input::NucleicSeqOrView{A}, 
    output::Union{IOStream, IOBuffer}; 
    alternative_start::Bool = false, 
    min_len::Int64 = 6
) where {A}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    println(output, "##gff-version 3\n##sequence-region Chr 1 $(length(input))")
    for (i, orf) in enumerate(orfs)
        id = string("ORF", lpad(string(i), padding, "0"))
        println(output, "Chr\t.\tORF\t", orf.location.start, "\t", orf.location.stop, "\t.\t", orf.strand, "\t.\tID=", id, ";Name=", id, ";Frame=", orf.frame)
        end
end

function write_orfs_gff(
    input::NucleicSeqOrView{A},
    output::String; 
    alternative_start::Bool = false, 
    min_len::Int64 = 6
) where {A}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(input))\n") 
        for (i, orf) in enumerate(orfs)
            id = string("ORF", lpad(string(i), padding, "0"))
            write(
                f,
                "Chr\t.\tORF\t$(orf.location.start)\t$(orf.location.stop)\t.\t$(orf.strand)\t.\tID=$(id);Name=$(id);Frame=$(orf.frame)\n",
            )
        end
    end
end

function write_orfs_gff(
    input::String, 
    output::Union{IOStream, IOBuffer}; 
    alternative_start::Bool = false, 
    min_len::Int64 = 6
)
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    println(output, "##gff-version 3\n##sequence-region Chr 1 $(length(input))")
    for (i, orf) in enumerate(orfs)
        id = string("ORF", lpad(string(i), padding, "0"))
        println(output, "Chr\t.\tORF\t", orf.location.start, "\t", orf.location.stop, "\t.\t", orf.strand, "\t.\tID=", id, ";Name=", id, ";Frame=", orf.frame)
    end
end

function write_orfs_gff(
    input::String,
    output::String; 
    alternative_start::Bool = false, 
    min_len::Int64 = 6
) 
    input = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(input))\n") 
        for (i, orf) in enumerate(orfs)
            id = string("ORF", lpad(string(i), padding, "0"))
            write(
                f,
                "Chr\t.\tORF\t$(orf.location.start)\t$(orf.location.stop)\t.\t$(orf.strand)\t.\tID=$(id);Name=$(id);Frame=$(orf.frame)\n",
            )
        end
    end
end

# Chr needs to be changed for fasta (?) or the name of seq (?) to fit well on IGV 