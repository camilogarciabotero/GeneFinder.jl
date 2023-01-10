"""
   write_bed(file::String, seq::LongDNA; kwargs...)

Write BED data to a file.

# Keywords

- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_bed(file::String, seq::LongDNA; alternative_start=false, min_len=6)
   open(file, "w") do f
      @simd for i in orf_finder(seq; alternative_start, min_len)
         write(f, "$(i.location.start)\t$(i.location.stop)\t$(i.strand)\n")
      end
   end
end

"""
   write_cds(file::String, seq::LongDNA; kwargs...)

Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `alternative_start`: A boolean value indicating whether alternative start codons should be used when identifying CDSs. Default is `false`.
- `min_len`: An integer representing the minimum length that a CDS must have in order to be included in the output file. Default is `6`.
"""
function write_cds(file::String, seq::LongDNA; alternative_start=false, min_len=6)
   open(file, "w") do f
      for i in cds_generator(seq; alternative_start, min_len)
         write(f, ">locus=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

"""
   write_proteins(file::String, seq::LongDNA; kwargs...)

Write the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

# Keywords

- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
- `min_len::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
"""
function write_proteins(file::String, seq::LongDNA; alternative_start=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len=6)
   open(file, "w") do f
      for i in cds_generator(seq; alternative_start, code, min_len)
         write(f, ">locus=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

# FASTA.Writer(open("some_file.fna", "w")) do writer
#     write(writer, record) # a FASTA.Record
# end