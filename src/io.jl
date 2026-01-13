export write_orfs_fna, write_orfs_faa, write_orfs_bed, write_orfs_gff

"""
    write_orfs_bed(input, output; kwargs...)

Write ORF coordinates in BED format.

BED (Browser Extensible Data) format is a tab-delimited text format commonly used
for genomic annotations. Each line contains: start, stop, strand, and frame.

# Arguments
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).

# Keywords
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
- `minlen::Int64=6`: Minimum ORF length in nucleotides.

# Output Format
```
start	stop	strand	frame
35	79	+	2
120	180	-	3
```

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file
write_orfs_bed(seq, "output.bed"; finder=NaiveFinder)

# Write to IO stream
open("output.bed", "w") do io
    write_orfs_bed(seq, io; finder=NaiveFinder)
end
```

See also: [`write_orfs_gff`](@ref), [`write_orfs_fna`](@ref), [`findorfs`](@ref)
"""
function write_orfs_bed(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    @inbounds for orf in collection
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
    collection = findorfs(input; finder, alternative_start, minlen)
    open(output, "w") do f
        @inbounds for orf in collection
            write(f, "$(leftposition(orf))\t$(rightposition(orf))\t$(strand(orf))\t$(frame(orf))\n")
        end
    end
end

"""
    write_orfs_fna(input, output; kwargs...)

Write ORF nucleotide sequences in FASTA format.

Outputs the DNA sequences of all detected ORFs, with headers containing
metadata (ID, coordinates, strand, frame, and features).

# Arguments
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).

# Keywords
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
- `minlen::Int64=6`: Minimum ORF length in nucleotides.

# Output Format
```
>ORF01 id=01 start=35 stop=79 strand=+ frame=2 features=[]
ATGCATGCATGCATGCTAG
>ORF02 id=02 start=120 stop=180 strand=- frame=3 features=[]
ATGCTAGCTAGCTAGCTAA
```

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file
write_orfs_fna(seq, "output.fna"; finder=NaiveFinder)

# Write to IO stream
open("output.fna", "w") do io
    write_orfs_fna(seq, io; finder=NaiveFinder, minlen=9)
end
```

See also: [`write_orfs_faa`](@ref), [`write_orfs_gff`](@ref), [`findorfs`](@ref)
"""
function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false, 
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    norfs = length(collection)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(collection)
        id = lpad(string(i), padding, "0")
        fts = isempty(features(orf)) ? "" : join(values(features(orf)), ",")
        strandchar = strand(orf) == PSTRAND ? '+' : '-'
        println(output, ">ORF", id, " id=", id, " start=", leftposition(orf), " stop=", rightposition(orf), " strand=", strandchar, " frame=", frame(orf), " features=[", fts, "]")
        println(output, sequence(collection, orf))
    end
end

function write_orfs_fna(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    norfs = length(collection)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(collection)
            id = lpad(string(i), padding, "0")
            fts = isempty(features(orf)) ? "" : join(values(features(orf)), ",")
            strandchar = strand(orf) == PSTRAND ? '+' : '-'
            write(f, ">ORF$(id) id=$(id) start=$(leftposition(orf)) stop=$(rightposition(orf)) strand=$(strandchar) frame=$(frame(orf)) features=[$(fts)]\n$(sequence(collection, orf))\n")
        end
    end
end

"""
    write_orfs_faa(input, output; kwargs...)

Write translated ORF sequences in FASTA format (amino acids).

Outputs the protein sequences of all detected ORFs, translating each ORF
using the specified genetic code. Headers contain the same metadata as `write_orfs_fna`.

# Arguments
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).

# Keywords
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
- `code::GeneticCode=ncbi_trans_table[1]`: Genetic code for translation (default: standard code).
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
- `minlen::Int64=6`: Minimum ORF length in nucleotides.

# Output Format
```
>ORF01 id=01 start=35 stop=79 strand=+ frame=2 features=[]
MHACA*
>ORF02 id=02 start=120 stop=180 strand=- frame=3 features=[]
MLALA*
```

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file with standard genetic code
write_orfs_faa(seq, "output.faa"; finder=NaiveFinder)

# Use bacterial genetic code (table 11)
write_orfs_faa(seq, "output.faa"; finder=NaiveFinder, code=ncbi_trans_table[11])

# Write to IO stream
open("output.faa", "w") do io
    write_orfs_faa(seq, io; finder=NaiveFinder)
end
```

See also: [`write_orfs_fna`](@ref), [`write_orfs_gff`](@ref), [`findorfs`](@ref), [`BioSequences.translate`](@extref)
"""
function write_orfs_faa(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream, IOBuffer};
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    code::GeneticCode = ncbi_trans_table[1]
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    norfs = length(collection)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    @inbounds for (i, orf) in enumerate(collection)
        id = lpad(string(i), padding, "0")
        fts = isempty(features(orf)) ? "" : join(values(features(orf)), ",")
        strandchar = strand(orf) == PSTRAND ? '+' : '-'
        println(output, ">ORF", id, " id=", id, " start=", leftposition(orf), " stop=", rightposition(orf), " strand=", strandchar, " frame=", frame(orf), " features=[", fts, "]")
        println(output, translate(sequence(collection, orf); code))
    end
end

function write_orfs_faa(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    alternative_start::Bool = false,
    minlen::Int64 = 6,
    code::GeneticCode = ncbi_trans_table[1]
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    norfs = length(collection)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        for (i, orf) in enumerate(collection)
            id = lpad(string(i), padding, "0")
            fts = isempty(features(orf)) ? "" : join(values(features(orf)), ",")
            strandchar = strand(orf) == PSTRAND ? '+' : '-'
            write(f, ">ORF$(id) id=$(id) start=$(leftposition(orf)) stop=$(rightposition(orf)) strand=$(strandchar) frame=$(frame(orf)) features=[$(fts)]\n$(translate(sequence(collection, orf); code))\n")
        end
    end
end

"""
    write_orfs_gff(input, output; kwargs...)

Write ORF annotations in GFF3 (General Feature Format version 3).

GFF3 is a standard format for genomic annotations, compatible with genome browsers
like IGV, JBrowse, and the UCSC Genome Browser.

# Arguments
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).

# Keywords
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
- `seqname::String="Chr"`: Sequence/chromosome name for the first GFF column.
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
- `minlen::Int64=6`: Minimum ORF length in nucleotides.

# Output Format
```
##gff-version 3
##sequence-region Chr 1 1000
Chr	.	ORF	35	79	.	+	.	ID=ORF01;Name=ORF01;Frame=2;Features=[]
Chr	.	ORF	120	180	.	-	.	ID=ORF02;Name=ORF02;Frame=3;Features=[]
```

# Example
```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file
write_orfs_gff(seq, "output.gff"; finder=NaiveFinder)

# Custom chromosome name for genome browser visualization
write_orfs_gff(seq, "output.gff"; finder=NaiveFinder, seqname="scaffold_1")

# Write to IO stream
open("output.gff", "w") do io
    write_orfs_gff(seq, io; finder=NaiveFinder)
end
```

See also: [`write_orfs_bed`](@ref), [`write_orfs_fna`](@ref), [`findorfs`](@ref)
"""
function write_orfs_gff(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::Union{IOStream,IOBuffer};
    finder::Type{F} = NaiveFinder,
    seqname::String = "Chr",
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    norfs = length(collection)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    println(output, "##gff-version 3")
    println(output, "##sequence-region ", seqname, " 1 ", length(input))
    @inbounds for (i, orf) in enumerate(collection)
        id = lpad(string(i), padding, "0")
        fts = isempty(features(orf)) ? "" : join(values(features(orf)), ",")
        strandchar = strand(orf) == PSTRAND ? '+' : '-'
        println(output, seqname, "\t.\tORF\t", leftposition(orf), "\t", rightposition(orf), "\t.\t", strandchar, "\t.\tID=ORF", id, ";Name=ORF", id, ";Frame=", frame(orf), ";Features=[", fts, "]")
    end
end

function write_orfs_gff(
    input::NucleicSeqOrView{DNAAlphabet{N}},
    output::String;
    finder::Type{F} = NaiveFinder,
    seqname::String = "Chr",
    alternative_start::Bool = false,
    minlen::Int64 = 6
) where {N,F<:GeneFinderMethod}
    collection = findorfs(input; finder, alternative_start, minlen)
    norfs = length(collection)
    padding = norfs < 10 ? length(string(norfs)) + 1 : length(string(norfs))
    open(output, "w") do f
        write(f, "##gff-version 3\n")
        write(f, "##sequence-region $(seqname) 1 $(length(input))\n")
        for (i, orf) in enumerate(collection)
            id = lpad(string(i), padding, "0")
            fts = isempty(features(orf)) ? "" : join(values(features(orf)), ",")
            strandchar = strand(orf) == PSTRAND ? '+' : '-'
            write(f, "$(seqname)\t.\tORF\t$(leftposition(orf))\t$(rightposition(orf))\t.\t$(strandchar)\t.\tID=ORF$(id);Name=ORF$(id);Frame=$(frame(orf));Features=[$(fts)]\n")
        end
    end
end