## Writting ORFs into bioinformatic formats

This package facilitates the creation of `FASTA`, `BED`, and `GFF` files, specifically extracting Open Reading Frame (ORF) information from `BioSequence` instances, particularly those of type `NucleicSeqOrView{A} where A`, and then writing the information into the desired format.

Functionality:

The package provides four distinct functions for writing files in different formats:

| Function          | Description                                            |
|-------------------|--------------------------------------------------------|
| `write_orfs_fna`    | Writes nucleotide sequences in FASTA format.     |
| `write_orfs_faa`    | Writes amino acid sequences in FASTA format.  |
| `write_orfs_bed`    | Outputs information in BED format.                           |
| `write_orfs_gff`    | Generates files in GFF format.                              |


All these functions support processing both `BioSequence` instances and external `FASTA` files. In the case of a `BioSequence` instace into external files, simply provide the path to the `FASTA` file using a `String` to the path. To demonstrate the use of the `write_*` methods with a `BioSequence`, consider the following example:

```julia
using BioSequences, GeneFinder

# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)
seq = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"
```

Once a `BioSequence` object has been instantiated, the `write_orfs_fna` function proves useful for generating a `FASTA` file containing the nucleotide sequences of the ORFs. Notably, the `write_orfs*` methods support either an `IOStream` or an `IOBuffer` as an output argument, allowing flexibility in directing the output either to a file or a buffer. In the following example, we demonstrate writing the output directly to a file.

```julia
outfile = "LFLS01000089.fna"

open(outfile, "w") do io
    write_orfs_fna(seq, io, NaiveFinder())
end
```

```bash
cat LFLS01000089.fna

>seq id=01 start=29 stop=40 strand=+ frame=2 features=[]
ATGCAACCCTGA
>seq id=02 start=137 stop=145 strand=+ frame=2 features=[]
ATGCGCTGA
>seq id=03 start=164 stop=184 strand=+ frame=2 features=[]
ATGCGTCGAATGGCACGGTGA
>seq id=04 start=173 stop=184 strand=+ frame=2 features=[]
ATGGCACGGTGA
>seq id=05 start=236 stop=241 strand=+ frame=2 features=[]
ATGTGA
>seq id=06 start=248 stop=268 strand=+ frame=2 features=[]
ATGTGTCCAACGGCAGTCTGA
>seq id=07 start=362 stop=373 strand=+ frame=2 features=[]
ATGCAACCCTGA
>seq id=08 start=470 stop=496 strand=+ frame=2 features=[]
ATGCACTGGCTGGTCCTGTCAATCTGA
>seq id=09 start=551 stop=574 strand=+ frame=2 features=[]
ATGTCACCGCACAAGGCAATGTGA
>seq id=10 start=569 stop=574 strand=+ frame=2 features=[]
ATGTGA
>seq id=11 start=581 stop=601 strand=+ frame=2 features=[]
ATGTGTCCAACGGCAGCCTGA
>seq id=12 start=695 stop=706 strand=+ frame=2 features=[]
ATGCAACCCTGA
```

This could also be done to writting a `FASTA` file with the nucleotide sequences of the ORFs using the `write_orfs_fna` function. Similarly for the `BED` and `GFF` files using the `write_orfs_bed` and `write_orfs_gff` functions respectively.