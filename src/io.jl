"""
   write_bed(file::String, seq::LongDNA; alternative_start::Bool=false)

Write BED data to a file.

Parameters:
- file: string; the file name to which the BED data will be written
- seq: an instance of `LongDNA`, representing a long DNA sequence
- alternative_start: boolean (optional); default is `false`

"""
function write_bed(file::String, seq::LongDNA; alternative_start=false, min_len=6)
   open(file, "w") do f
      @simd for i in simplefind(seq; alternative_start, min_len)
         write(f, "$(i.location.start)\t$(i.location.stop)\t$(i.strand)\n")
      end
   end
end

"""
   write_cds(file::String, seq::LongDNA; alternative_start=false, min_len = 6)

Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.

Parameters:
- `file`: A string representing the file path and name where the CDSs should be written.
- `seq`: A `LongDNA` object representing the DNA sequence from which the CDSs should be extracted.

Keyword Arguments:
- `alternative_start`: A boolean value indicating whether alternative start codons should be used when identifying CDSs. Default is `false`.
- `min_len`: An integer representing the minimum length that a CDS must have in order to be included in the output file. Default is `6`.
"""
function write_cds(file::String, seq::LongDNA; alternative_start=false, min_len=6)
   open(file, "w") do f
      for i in simplecds_generator(seq; alternative_start, min_len)
         write(f, ">locus=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

"""
   write_proteins(file::String, seq::LongDNA; alternative_start = false, code::GeneticCode = BioSequences.standard_genetic_code, min_len = 6)

Write a file containing the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

Parameters:
- `file`: A string representing the file path and name where the protein sequences should be written.
- `seq`: A `LongDNA` object representing the DNA sequence from which the CDSs and protein sequences should be extracted.

Keyword Arguments:
- `alternative_start`: A boolean value indicating whether alternative start codons should be used when identifying CDSs and translating them into protein sequences. Default is `false`.
- `code`: A `GeneticCode` object representing the genetic code that should be used to translate the CDSs into protein sequences. Default is the standard genetic code.
- `min_len`: An integer representing the minimum length that a protein sequence must have in order to be included in the output file. Default is `6`.
"""
function write_proteins(file::String, seq::LongDNA; alternative_start=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len=6)
   open(file, "w") do f
      for i in simpleprot_generator(seq; alternative_start, code, min_len)
         write(f, ">locus=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

# FASTA.Writer(open("some_file.fna", "w")) do writer
#     write(writer, record) # a FASTA.Record
# end