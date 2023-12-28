# const re = Automa.RegExp

# start = re"ATG"
# codons = re"[ACGT]+"
# stop = re"(TAA|TAG|TGA)"

# orf = re.cat(start, codons, stop)

# machine = compile(orf)

# @eval function orfparser(data)
#     seqs = String[]
#     $(generate_code(machine))
#     return seqs
# end

# orfparser("ATG")