"""
	write_bed(input::LongDNA, output::String; kwargs...)
	write_bed(input::String, output::String; kwargs...)

Write BED data to a file.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_bed(input::LongDNA, output::String; alternative_start=false, min_len=6)
   open(output, "w") do f
      @simd for i in findorfs(input; alternative_start, min_len)
         write(f, "$(i.location.start)\t$(i.location.stop)\t$(i.strand)\n")
      end
   end
end

function write_bed(input::String, output::String; alternative_start=false, min_len=6)
   rdr = FASTA.Reader(open(input))
   record = first(rdr)
   seq = sequence(record)
   dnaseq = LongDNA{4}(seq)

   open(output, "w") do f
      @simd for i in findorfs(dnaseq; alternative_start, min_len)
         write(f, "$(i.location.start)\t$(i.location.stop)\t$(i.strand)\n")
      end
   end
end

"""
	write_cds(input::LongDNA, output::String; kwargs...)
	write_cds(input::String, output::String; kwargs...)

Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `alternative_start`: A boolean value indicating whether alternative start codons should be used when identifying CDSs. Default is `false`.
- `min_len`: An integer representing the minimum length that a CDS must have in order to be included in the output file. Default is `6`.
"""
function write_cds(input::LongDNA, output::String; alternative_start=false, min_len=6)
   open(output, "w") do f
      for i in cdsgenerator(input; alternative_start, min_len)
         write(f, ">location=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

function write_cds(input::String, output::String; alternative_start=false, min_len=6)
   rdr = FASTA.Reader(open(input))
   record = first(rdr)
   seq = sequence(record)
   dnaseq = LongDNA{4}(seq)

   open(output, "w") do f
      for i in cdsgenerator(dnaseq; alternative_start, min_len)
         write(f, ">location=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

"""
	write_proteins(input::LongDNA, output::String; kwargs...)

Write the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_proteins(input::LongDNA ,output::String; alternative_start=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len=6)
   open(output, "w") do f
      for i in proteingenerator(input; alternative_start, code, min_len)
         write(f, ">location=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

function write_proteins(input::String, output::String; alternative_start=false, min_len=6)
   rdr = FASTA.Reader(open(input))
   record = first(rdr)
   seq = sequence(record)
   dnaseq = LongDNA{4}(seq)

   open(output, "w") do f
      for i in proteingenerator(dnaseq; alternative_start, min_len)
         write(f, ">location=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

"""
	write_gff(input::LongDNA, output::String; kwargs...)
	write_gff(input::String, output::String; kwargs...)

Write GFF data to a file.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_gff(input::String, output::String; alternative_start=false, min_len=6)
   rdr = FASTA.Reader(open(input))
   record = first(rdr)
   seq = sequence(record)
   dnaseq = LongDNA{4}(seq)

   open(output, "w") do f
      write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(dnaseq))\n")
      for (index, i) in enumerate(findorfs(dnaseq; alternative_start, min_len))
         id = string("ORF", lpad(string(index), 5, "0"))
         write(f, "Chr\tGeneFinder\tORF\t$(i.location.start)\t$(i.location.stop)\t.\t$(i.strand)\t.\tid=$(id);name=$(id)\n")
      end
   end
end

function write_gff(input::LongDNA, output::String; alternative_start=false, min_len=6)
   open(output, "w") do f
      write(f, "##gff-version 3\n##sequence-region Chr 1 $(length(input))\n")
      for (index, i) in enumerate(findorfs(input; alternative_start, min_len))
         id = string("ORF", lpad(string(index), 5, "0"))
         write(f, "Chr\tGeneFinder\tORF\t$(i.location.start)\t$(i.location.stop)\t.\t$(i.strand)\t.\tid=$(id);name=$(id)\n")
      end
   end
end