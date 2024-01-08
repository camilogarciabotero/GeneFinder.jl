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
    write_orfs_dna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}
    write_orfs_dna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}
    write_orfs_dna(input::String, output::Union{IOStream, IOBuffer}; kwargs...)
    write_orfs_dna(input::String, output::String; kwargs...)

Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
- `min_len::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.
# Examples

```julia
filename = "output.fna"

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

open(filename, "w") do file
     write_orfs_dna(seq, file)
end
```
"""
function write_orfs_dna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    alternative_start = false, 
    min_len = 6
) where {N}
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (index, i) in enumerate(orfs)
        id = string(lpad(string(index), padding, "0"))
        # sequence = i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
        println(output, ">ORF$(id) location=$(i.location) strand=$(i.strand) frame=$(i.frame)\n$(i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location]))")
    end
end

function write_orfs_dna(
    input::String,
    output::Union{IOStream, IOBuffer};
    alternative_start = false, 
    min_len = 6
)
    dnaseq = fasta_to_dna(input)[1]
    orfs = findorfs(input; alternative_start, min_len)
    norfs = length(orfs)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (index, i) in enumerate(orfs)
        id = string(lpad(string(index), padding, "0"))
        # sequence = i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location])
        println(output, ">ORF$(id) location=$(i.location) strand=$(i.strand) frame=$(i.frame)\n$(i.strand == '+' ? dnaseq[i.location] : reverse_complement(@view dnaseq[i.location]))")
    end
end

function write_orfs_dna(
    input::NucleicSeqOrView{DNAAlphabet{N}}, 
    output::String; 
    alternative_start = false,
    min_len = 6
) where {N}
    open(output, "w") do f
        for i in findorfs(input; alternative_start, min_len)
            write(f, ">location=$(i.location) strand=$(i.strand) frame=$(i.frame)\n$(i.strand == '+' ? input[i.location] : reverse_complement(@view input[i.location]))\n")
        end
    end
end

function write_orfs_dna(
    input::String, 
    output::String; 
    alternative_start = false,
    min_len = 6
)
    dnaseq = fasta_to_dna(input)[1]

    open(output, "w") do f
        for i in findorfs(dnaseq; alternative_start, min_len)
            # sequence = i.strand == '+' ? dnaseq[i.location] : reverse_complement(@view dnaseq[i.location])
            write(f, ">location=$(i.location) strand=$(i.strand) frame=$(i.frame)\n$(i.strand == '+' ? dnaseq[i.location] : reverse_complement(@view dnaseq[i.location]))\n")
        end
    end
end

"""
	write_orfs_aa(input::LongSequence{DNAAlphabet{4}}, output::String; kwargs...)

Write the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_orfs_aa(
    input::LongSequence{DNAAlphabet{N}},
    output::String;
    alternative_start = false,
    code::GeneticCode = ncbi_trans_table[1],
    min_len = 6,
) where {N}
    open(output, "w") do f
        for i in findorfs(input; alternative_start, min_len)
            translation = i.strand == '+' ? translate(input[i.location]; code) : translate(reverse_complement(@view input[i.location]); code)
            write(f, ">location=$(i.location) strand=$(i.strand) frame=$(i.frame)\n$(translation)\n")
        end
    end
end

function write_orfs_aa(
    input::String,
    output::String;
    alternative_start = false,
    min_len = 6,
)
    dnaseq = fasta_to_dna(input)[1]

    open(output, "w") do f
        for i in findorfs(dnaseq; alternative_start, min_len)
            translation = i.strand == '+' ? translate(dnaseq[i.location]) : translate(reverse_complement(@view dnaseq[i.location]))
            write(f, ">location=$(i.location) strand=$(i.strand) frame=$(i.frame)\n$(translation)\n")
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

function write_orfs_gff(input::LongSequence{DNAAlphabet{4}}, output::String; alternative_start = false, min_len = 6)
    open(output, "w") do f
        write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(input))\n") 
        for (index, i) in enumerate(findorfs(input; alternative_start, min_len))
            id = string("ORF", lpad(string(index), 5, "0"))
            write(
                f,
                "Chr\t.\tORF\t$(i.location.start)\t$(i.location.stop)\t.\t$(i.strand)\t.\tID=$(id);Name=$(id);Frame=$(i.frame)\n",
            )
        end
    end
end

# Chr needs to be changed for fasta (?) or the name of seq (?) to fit well on IGV 