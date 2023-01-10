const STOPCODONS = [Codon("TAG"), Codon("TAA"), Codon("TGA")]
const STARTCODON = ExactSearchQuery(Codon("ATG"), iscompatible)
const EXTENDED_STARTCODONS = PWMSearchQuery([Codon("ATG"), Codon("GTG"), Codon("TTG")], 1.0)