"""
	write_orfs_bed(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...)
	write_orfs_bed(input::String, output::String; kwargs...)

Write BED data to a file.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_orfs_bed(input::LongSequence{DNAAlphabet{4}}, output::String; alternative_start = false, min_len = 6)
    open(output, "w") do f
        @simd for i in findorfs(input; alternative_start, min_len)
            write(f, "$(i.location.start)\t$(i.location.stop)\t$(i.strand)\t$(i.frame)\n")
        end
    end
end

function write_orfs_bed(input::String, output::String; alternative_start = false, min_len = 6)
    dnaseq = fasta_to_dna(input)[1]

    open(output, "w") do f
        @simd for i in findorfs(dnaseq; alternative_start, min_len)
            write(f, "$(i.location.start)\t$(i.location.stop)\t$(i.strand)\t$(i.frame)\n")
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
    alternative_start = false, 
    min_len = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(input[orf])") # i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
    end
end

function write_orfs_fna(
    input::String,
    output::Union{IOStream, IOBuffer};
    alternative_start = false, 
    min_len = 6
)
    dnaseq = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(dnaseq[orf])") # i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
    end
end

function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}}, 
    output::String; 
    alternative_start = false,
    min_len = 6
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
    output::String; 
    alternative_start = false,
    min_len = 6
)
    dnaseq = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(dnaseq[orf])\n") # i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
        end
    end
end

"""
	write_orfs_faa(input::LongSequence{DNAAlphabet{4}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}
	write_orfs_faa(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...) where {N}
	write_orfs_faa(input::String, output::Union{IOStream, IOBuffer}; kwargs...)
	write_orfs_faa(input::String, output::Union{IOStream, IOBuffer}; kwargs...)


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
    alternative_start = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(i), padding, "0"))
        println(output, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(translate(input[orf]; code))")
    end
end

function write_orfs_faa(
    input::String,
    output::Union{IOStream, IOBuffer};
    alternative_start = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len = 6
)
    dnaseq = fasta_to_dna(input)[1]
    orfs = findorfs(dnaseq; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(orfs)
        id = string(lpad(string(index), padding, "0"))
        println(output, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(translate(dnaseq[orf]; code))")
    end
end

function write_orfs_faa(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    alternative_start = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len = 6,
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
    output::String;
    alternative_start = false,
    min_len = 6,
)
    dnaseq = fasta_to_dna(input)[1]
    orfs = findorfs(dnaseq; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(orfs)
            id = string(lpad(string(i), padding, "0"))
            write(f, ">ORF$(id) id=$(id) start=$(orf.location.start) stop=$(orf.location.stop) strand=$(orf.strand) frame=$(orf.frame)\n$(translate(dnaseq[orf]; code))\n")
        end
    end
end

"""
	write_orfs_gff(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...)
	write_orfs_gff(input::String, output::String; kwargs...)

Write GFF data to a file.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_orfs_gff(input::LongSequence{DNAAlphabet{4}}, output::String; alternative_start = false, min_len = 6)
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

function write_orfs_gff(input::String, output::String; alternative_start = false, min_len = 6)
    dnaseq = fasta_to_dna(input)[1]

    open(output, "w") do f
        write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(dnaseq))\n")
        for (index, i) in enumerate(findorfs(dnaseq; alternative_start, min_len))
            id = string("ORF", lpad(string(index), 5, "0"))
            write(
                f,
                "Chr\t.\tORF\t$(i.location.start)\t$(i.location.stop)\t.\t$(i.strand)\t.\tID=$(id);Name=$(id);Frame=$(i.frame)\n",
            )
        end
    end
end

# Chr needs to be changed for fasta (?) or the name of seq (?) to fit well on IGV 