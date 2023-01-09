
function write_cds(file::String, seq::LongDNA; alternative_start=false, min_len = 6)
    open(file, "w") do f
       for i in simplecds_generator(seq; alternative_start, min_len)
          write(f, ">locus=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
       end
    end
end

function write_proteins(file::String, seq::LongDNA; alternative_start = false, code::GeneticCode = BioSequences.standard_genetic_code, min_len = 6)
   open(file, "w") do f
      for i in simpleprot_generator(seq; alternative_start, code, min_len)
         write(f, ">locus=$(i.orf.location) strand=$(i.orf.strand)\n$(i.sequence)\n")
      end
   end
end

# FASTA.Writer(open("some_file.fna", "w")) do writer
#     write(writer, record) # a FASTA.Record
# end