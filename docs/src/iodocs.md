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
    write_orfs_fna(seq, io)
end
```

```bash
cat LFLS01000089.fna

>ORF01 id=01 start=29 stop=40 strand=+ frame=2
ATGCAACCCTGA
>ORF02 id=02 start=137 stop=145 strand=+ frame=2
ATGCGCTGA
>ORF03 id=03 start=164 stop=184 strand=+ frame=2
ATGCGTCGAATGGCACGGTGA
>ORF04 id=04 start=173 stop=184 strand=+ frame=2
ATGGCACGGTGA
>ORF05 id=05 start=236 stop=241 strand=+ frame=2
ATGTGA
>ORF06 id=06 start=248 stop=268 strand=+ frame=2
ATGTGTCCAACGGCAGTCTGA
>ORF07 id=07 start=362 stop=373 strand=+ frame=2
ATGCAACCCTGA
>ORF08 id=08 start=470 stop=496 strand=+ frame=2
ATGCACTGGCTGGTCCTGTCAATCTGA
>ORF09 id=09 start=551 stop=574 strand=+ frame=2
ATGTCACCGCACAAGGCAATGTGA
>ORF10 id=10 start=569 stop=574 strand=+ frame=2
ATGTGA
>ORF11 id=11 start=581 stop=601 strand=+ frame=2
ATGTGTCCAACGGCAGCCTGA
>ORF12 id=12 start=695 stop=706 strand=+ frame=2
ATGCAACCCTGA
```

### Combining `FASTX` for reading and writing fastas

We can now combine the `FASTX` package with the function `write_orfs_faa` to write a `FASTA` file with the protein sequences of the translated ORFs obtained from an external `FASTA` file. 

```julia
infile = "test/data/NC_001884.fasta"
outfile = "test/data/NC_001884-orfs.faa"

open(inputfile) do io
    write_orfs_faa(infile, outfile)
end
```

```bash
head test/data/NC_001884-orfs.faa

>ORF0001 id=0001 start=41 stop=145 strand=- frame=2
MTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
>ORF0002 id=0002 start=41 stop=172 strand=- frame=2
MVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*
>ORF0003 id=0003 start=41 stop=454 strand=- frame=2
MSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDL
TATNSFH*
>ORF0004 id=0004 start=41 stop=472 strand=- frame=2
MKTKKQMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRF
NFIFDLTATNSFH*
>ORF0005 id=0005 start=41 stop=505 strand=- frame=2
MLSKYEDDNSNMKTKKQMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIP
AQFEITPILRFNFIFDLTATNSFH*
```

This could also be done to writting a `FASTA` file with the nucleotide sequences of the ORFs using the `write_orfs_fna` function. Similarly for the `BED` and `GFF` files using the `write_orfs_bed` and `write_orfs_gff` functions respectively.