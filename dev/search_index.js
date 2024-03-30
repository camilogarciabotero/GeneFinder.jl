var documenterSearchIndex = {"docs":
[{"location":"iodocs/#Writting-ORF-information-into-bioinformatic-formats","page":"Wrtiting ORFs in files","title":"Writting ORF information into bioinformatic formats","text":"","category":"section"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"This package facilitates the creation of FASTA, BED, and GFF files, specifically extracting Open Reading Frame (ORF) information from BioSequence instances, particularly those of type NucleicSeqOrView{A} where A, and then writing the information into the desired format.","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"Functionality:","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"The package provides four distinct functions for writing files in different formats:","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"Function Description\nwrite_orfs_fna Writes nucleotide sequences in FASTA format.\nwrite_orfs_faa Writes amino acid sequences in FASTA format.\nwrite_orfs_bed Outputs information in BED format.\nwrite_orfs_gff Generates files in GFF format.","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"All these functions support processing both BioSequence instances and external FASTA files. In the case of a BioSequence instace into external files, simply provide the path to the FASTA file using a String to the path. To demonstrate the use of the write_* methods with a BioSequence, consider the following example:","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"using BioSequences, GeneFinder\n\n# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)\nseq = dna\"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC\"","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"Once a BioSequence object has been instantiated, the write_orfs_fna function proves useful for generating a FASTA file containing the nucleotide sequences of the ORFs. Notably, the write_orfs* methods support either an IOStream or an IOBuffer as an output argument, allowing flexibility in directing the output either to a file or a buffer. In the following example, we demonstrate writing the output directly to a file.","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"outfile = \"LFLS01000089.fna\"\n\nopen(outfile, \"w\") do io\n    write_orfs_fna(seq, io)\nend","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"cat LFLS01000089.fna\n\n>ORF01 id=01 start=29 stop=40 strand=+ frame=2\nATGCAACCCTGA\n>ORF02 id=02 start=137 stop=145 strand=+ frame=2\nATGCGCTGA\n>ORF03 id=03 start=164 stop=184 strand=+ frame=2\nATGCGTCGAATGGCACGGTGA\n>ORF04 id=04 start=173 stop=184 strand=+ frame=2\nATGGCACGGTGA\n>ORF05 id=05 start=236 stop=241 strand=+ frame=2\nATGTGA\n>ORF06 id=06 start=248 stop=268 strand=+ frame=2\nATGTGTCCAACGGCAGTCTGA\n>ORF07 id=07 start=362 stop=373 strand=+ frame=2\nATGCAACCCTGA\n>ORF08 id=08 start=470 stop=496 strand=+ frame=2\nATGCACTGGCTGGTCCTGTCAATCTGA\n>ORF09 id=09 start=551 stop=574 strand=+ frame=2\nATGTCACCGCACAAGGCAATGTGA\n>ORF10 id=10 start=569 stop=574 strand=+ frame=2\nATGTGA\n>ORF11 id=11 start=581 stop=601 strand=+ frame=2\nATGTGTCCAACGGCAGCCTGA\n>ORF12 id=12 start=695 stop=706 strand=+ frame=2\nATGCAACCCTGA","category":"page"},{"location":"iodocs/#Combining-FASTX-for-reading-and-writing-fastas","page":"Wrtiting ORFs in files","title":"Combining FASTX for reading and writing fastas","text":"","category":"section"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"We can now combine the FASTX package with the function write_orfs_faa to write a FASTA file with the protein sequences of the translated ORFs obtained from an external FASTA file. ","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"infile = \"test/data/NC_001884.fasta\"\noutfile = \"test/data/NC_001884-orfs.faa\"\n\nopen(inputfile) do io\n    write_orfs_faa(infile, outfile)\nend","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"head test/data/NC_001884-orfs.faa\n\n>ORF0001 id=0001 start=41 stop=145 strand=- frame=2\nMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*\n>ORF0002 id=0002 start=41 stop=172 strand=- frame=2\nMVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDLTATNSFH*\n>ORF0003 id=0003 start=41 stop=454 strand=- frame=2\nMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRFNFIFDL\nTATNSFH*\n>ORF0004 id=0004 start=41 stop=472 strand=- frame=2\nMKTKKQMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIPAQFEITPILRF\nNFIFDLTATNSFH*\n>ORF0005 id=0005 start=41 stop=505 strand=- frame=2\nMLSKYEDDNSNMKTKKQMSEHLSQKEKELKNKENFIFDKYESGIYSDELFLKRKAALDEEFKELQNAKNELNGLQDTQSEIDSNTVRNNINKIIDQYHIESSSEKKNELLRMVLKDVIVNMTQKRKGPIP\nAQFEITPILRFNFIFDLTATNSFH*","category":"page"},{"location":"iodocs/","page":"Wrtiting ORFs in files","title":"Wrtiting ORFs in files","text":"This could also be done to writting a FASTA file with the nucleotide sequences of the ORFs using the write_orfs_fna function. Similarly for the BED and GFF files using the write_orfs_bed and write_orfs_gff functions respectively.","category":"page"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = GeneFinder\nDocTestSetup = quote\n    using GeneFinder\nend","category":"page"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [GeneFinder]","category":"page"},{"location":"api/#GeneFinder.ORF","page":"API","title":"GeneFinder.ORF","text":"struct ORF\n    location::UnitRange{Int64}\n    strand::Char\n    frame::Int\n    score::Union{Nothing, Float64}\nend\n\nThe ORF struct represents an open reading frame in a DNA sequence. It has two fields: \n\nlocation: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence\nstrand:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.\nframe: is an Int type indicating the reading frame of the ORF. The frame is the position of the first nucleotide of the codon that starts the ORF, relative to the start of the sequence. It can be 1, 2, or 3.\nscore: is a Union{Nothing, Float64} type indicating the score of the ORF. It can be a Float64 or nothing.\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder._locationiterator-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N","page":"API","title":"GeneFinder._locationiterator","text":"locationiterator(sequence::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}\n\nThis is an iterator function that uses regular expressions to search the entire ORF (instead of start and stop codons) in a LongSequence{DNAAlphabet{4}} sequence.     It uses an anonymous function that will find the first regularly expressed ORF. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.\n\nnote: Note\nAs a note of the implementation we want to expand on how the ORFs are found:The expression (?:[N]{3})*? serves as the boundary between the start and stop codons. Within this expression, the character class [N]{3} captures exactly three occurrences of any character (representing nucleotides using IUPAC codes). This portion functions as the regular codon matches. Since it is enclosed within a non-capturing group (?:) and followed by *?, it allows for the matching of intermediate codons, but with a preference for the smallest number of repetitions. In summary, the regular expression ATG(?:[N]{3})*?T(AG|AA|GA) identifies patterns that start with \"ATG,\" followed by any number of three-character codons (represented by \"N\" in the IUPAC code), and ends with a stop codon \"TAG,\" \"TAA,\" or \"TGA.\" This pattern is commonly used to identify potential protein-coding regions within genetic sequences.See more about the discussion here\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N","page":"API","title":"GeneFinder.findorfs","text":"findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; findermethod=naivefinder, alternative_start=false, min_len=6) where {N}\n\nThis is the main interface method for finding open reading frames (ORFs) in a DNA sequence.\n\nIt takes the following arguments:\n\nsequence: The nucleic acid sequence to search for ORFs.\nfindermethod: The algorithm used to find ORFs. Default is naivefinder.\nalternative_start: A boolean indicating whether to consider alternative start codons. Default is false.\nmin_len: The minimum length of an ORF. Default is 6.\n\nReturns a vector of ORF objects representing the found ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.get_orfs_aa-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N","page":"API","title":"GeneFinder.get_orfs_aa","text":"get_orfs_aa(sequence::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N}\n\nThis function takes a NucleicSeqOrView{DNAAlphabet{N}} sequence and identifies the open reading frames (ORFs) using the findorfs() function. The function then translates each ORF into an amino acid sequence and stores it in a Vector{LongSubSeq{AminoAcidAlphabet}}.\n\nArguments\n\nsequence: The input sequence as a NucleicSeqOrView{DNAAlphabet{N}}\n\nKeyword Arguments\n\nalternative_start::Bool=false: If set to true, the function considers alternative start codons when searching for ORFs. This increases the execution time by approximately 3x.\nmin_len::Int64=6: The minimum length of the allowed ORF. By default, it allows ORFs that can encode at least one amino acid (e.g., aa\"M*\").\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.get_orfs_dna-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N","page":"API","title":"GeneFinder.get_orfs_dna","text":"get_orfs_dna(sequence::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}\n\nThis function takes a NucleicSeqOrView{DNAAlphabet{N}} sequence and identifies the open reading frames (ORFs) using the findorfs() function. The function then extracts the DNA sequence of each ORF and stores it in a Vector{LongSubSeq{DNAAlphabet{4}}}.\n\nArguments\n\nsequence: The input sequence as a NucleicSeqOrView{DNAAlphabet{N}}\n\nKeyword Arguments\n\nalternative_start::Bool=false: If set to true, the function considers alternative start codons when searching for ORFs. This increases the execution time by approximately 3x.\nmin_len::Int64=6: The minimum length of the allowed ORF. By default, it allows ORFs that can encode at least one amino acid (e.g., aa\"M*\").\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.naivefinder-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N","page":"API","title":"GeneFinder.naivefinder","text":"naivefinder(sequence::NucleicAlphabet{DNAAlphabet{N}}; alternativestart::Bool=false, minlen::Int64=6)::Vector{ORF} where {N}\n\nA simple implementation that finds ORFs in a DNA sequence.\n\nThe naivefinder function takes a LongSequence{DNAAlphabet{4}} sequence and returns a Vector{ORF} containing the ORFs found in the sequence.      It searches entire regularly expressed CDS, adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.         Extending the starting codons with the alternative_start = true will search for ATG, GTG, and TTG.     Some studies have shown that in E. coli (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.\n\nnote: Note\nThis function has not ORFs scoring scheme. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_orfs_bed-Union{Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where N","page":"API","title":"GeneFinder.write_orfs_bed","text":"write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}\nwrite_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}\nwrite_orfs_bed(input::String, output::Union{IOStream, IOBuffer}; kwargs...)\nwrite_orfs_bed(input::String, output::String; kwargs...)\n\nWrite BED data to a file.\n\nArguments\n\ninput::NucleicAcidAlphabet{DNAAlphabet{N}}: The input DNA sequence.\noutput::String: The output file path.\nalternative_start::Bool=false: If true, alternative start codons will be used when identifying CDSs. Default is false.\nmin_len::Int64=6: The minimum length that a CDS must have in order to be included in the output file. Default is 6.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_orfs_faa-Union{Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where N","page":"API","title":"GeneFinder.write_orfs_faa","text":"write_orfs_faa(input::NucleicSeqOrView{DNAAlphabet{4}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}\nwrite_orfs_faa(input::NucleicSeqOrView{DNAAlphabet{4}}, output::String; kwargs...) where {N}\nwrite_orfs_faa(input::String, output::Union{IOStream, IOBuffer}; kwargs...)\nwrite_orfs_faa(input::String, output::String; kwargs...)\n\nWrite the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.\n\nKeywords\n\ncode::GeneticCode=BioSequences.standard_genetic_code: The genetic code by which codons will be translated. See BioSequences.ncbi_trans_table for more info. \nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\nExamples\n\nfilename = \"output.faa\"\n\nseq = dna\"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA\"\n\nopen(filename, \"w\") do file\n     write_orfs_faa(seq, file)\nend\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_orfs_fna-Union{Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where N","page":"API","title":"GeneFinder.write_orfs_fna","text":"write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}\nwrite_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}\nwrite_orfs_fna(input::String, output::Union{IOStream, IOBuffer}; kwargs...)\nwrite_orfs_fna(input::String, output::String; kwargs...)\n\nWrite a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.\n\nKeywords\n\nalternative_start::Bool=false: If true, alternative start codons will be used when identifying CDSs. Default is false.\nmin_len::Int64=6: The minimum length that a CDS must have in order to be included in the output file. Default is 6.\n\nExamples\n\nfilename = \"output.fna\"\n\nseq = dna\"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA\"\n\nopen(filename, \"w\") do file\n     write_orfs_fna(seq, file)\nend\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_orfs_gff-Union{Tuple{A}, Tuple{Union{BioSequences.LongSequence{A}, BioSequences.LongSubSeq{A}}, Union{IOStream, IOBuffer}}} where A","page":"API","title":"GeneFinder.write_orfs_gff","text":"write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}; kwargs...) where {N}\nwrite_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String; kwargs...) where {N}\nwrite_orfs_gff(input::String, output::Union{IOStream, IOBuffer}; kwargs...)\nwrite_orfs_gff(input::String, output::String; kwargs...)\n\nWrite GFF data to a file.\n\nArguments\n\ninput: The input DNA sequence.\noutput: The output file to write the GFF data to.\nalternative_start::Bool=false: If true, extended start codons will be considered during the search, increasing the execution time.\nmin_len::Int64=6: The minimum length of the allowed ORF. The default value allows for possible encoding proteins with the aa\"M*\" sequence.\n\n\n\n\n\n","category":"method"},{"location":"simplecodingrule/#The-*log-odds-ratio*-decision-rule","page":"-","title":"The log-odds ratio decision rule","text":"","category":"section"},{"location":"simplecodingrule/","page":"-","title":"-","text":"The sequence probability given a transition probability model (eq. 2) could be used as the source of a sequence classification based on a decision rule to classify whether a sequence correspond to a model or another. Now, imagine we got two DNA sequence transition models, a CDS model and a No-CDS model. The log-odds ratio decision rule could be establish as:","category":"page"},{"location":"simplecodingrule/","page":"-","title":"-","text":"beginalign\nS(X) = log fracP_C(X_1=i_1 ldots X_T=i_T)P_N(X_1=i_1 ldots X_T=i_T)  begincases  eta  Rightarrow textcoding   eta  Rightarrow textnoncoding endcases\nendalign","category":"page"},{"location":"simplecodingrule/","page":"-","title":"-","text":"Where the P_C is the probability of the sequence given a CDS model, P_N is the probability of the sequence given a No-CDS model, the decision rule is finally based on whether the ratio is greater or lesser than a given threshold η of significance level.","category":"page"},{"location":"simplecodingrule/","page":"-","title":"-","text":"In the GeneFinder we have implemented this rule and a couple of basic transition probability models of CDS and No-CDS of E. coli from Axelson-Fisk (2015) work. To check whether a random sequence could be coding based on these decision we use the predicate iscoding with the ECOLICDS and ECOLINOCDS models:","category":"page"},{"location":"simplecodingrule/","page":"-","title":"-","text":"using GeneFinder, BioSequences\nrandseq = get_orfs_dna(randdnaseq(99))[1] # this will retrieved a random coding ORF\n\niscoding(randseq, ECOLICDS, ECOLINOCDS)","category":"page"},{"location":"simplecodingrule/","page":"-","title":"-","text":"true","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"engine: knitr cache: true –-","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"<p align=\"center\">   <img src=\"../assets/logo.svg\" height=\"150\"><br/>   <i>A Gene Finder framework for Julia.</i><br/><br/>   <a href=\"https://www.repostatus.org/#wip\">     <img src=\"https://www.repostatus.org/badges/latest/wip.svg\">   </a>   <a href=\"https://codecov.io/gh/camilogarciabotero/GeneFinder.jl\">     <img src=\"https://img.shields.io/codecov/c/github/camilogarciabotero/GeneFinder.jl?logo=codecov&logoColor=white\">   </a>   <a href=\"https://camilogarciabotero.github.io/GeneFinder.jl/dev/\">     <img src=\"https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white\">   </a>   <a href=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl\">     <img src=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main\">   <a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE\">     <img src=\"https://img.shields.io/badge/license-MIT-green.svg\">   </a> </p>","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"","category":"page"},{"location":"oldindex/#Overview","page":"-","title":"Overview","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"oldindex/#Installation","page":"-","title":"Installation","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"add GeneFinder","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"oldindex/#Algorithms","page":"-","title":"Algorithms","text":"","category":"section"},{"location":"oldindex/#Coding-genes-(CDS-ORFs)","page":"-","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"[x] Simple finder\n[ ] EasyGene\n[ ] GLIMMER\n[ ] Prodigal - Pyrodigal\n[ ] PHANOTATE\n[ ] k-mer based gene finders (?)\n[ ] Augustus (?)","category":"page"},{"location":"oldindex/#Non-coding-genes-(RNA)","page":"-","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"[ ] Infernal\n[ ] tRNAscan","category":"page"},{"location":"oldindex/#Other-features","page":"-","title":"Other features","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"[ ] parallelism SIMD ?\n[ ] memory management (?)\n[ ] specialized types\n[x] Gene\n[x] ORF\n[x] CDS\n[ ] EukaryoticGene (?)\n[ ] ProkaryoticGene (?)\n[ ] Codon\n[ ] Intron\n[ ] Exon\n[ ] GFF –> See other packages\n[ ] FASTX –> See I/O in other packages","category":"page"},{"location":"oldindex/#Compatibilities","page":"-","title":"Compatibilities","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"Must interact with or extend:","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl","category":"page"},{"location":"oldindex/#Contributing","page":"-","title":"Contributing","text":"","category":"section"},{"location":"oldindex/#Citing","page":"-","title":"Citing","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"Logo: gene analysis by Vector Points from the Noun Project","category":"page"},{"location":"roadmap/#Roadmap","page":"Roadmap","title":"Roadmap","text":"","category":"section"},{"location":"roadmap/#Coding-genes-(CDS-ORFs)","page":"Roadmap","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"☒ Finding ORFs\n☐ EasyGene\n☐ GLIMMER\n☐ Prodigal - Pyrodigal\n☐ PHANOTATE\n☐ k-mer based gene finders (?)\n☐ Augustus (?)","category":"page"},{"location":"roadmap/#Non-coding-genes-(RNA)","page":"Roadmap","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"☐ Infernal\n☐ tRNAscan","category":"page"},{"location":"roadmap/#Other-features","page":"Roadmap","title":"Other features","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"☐ parallelism SIMD ?\n☐ memory management (?)\n☐ incorporate Ribosime Binding Sites (RBS)\n☐ incorporate Programmed Reading Frame Shifting (PRFS)\n☐ specialized types\n☒ Gene\n☒ ORF\n☒ Codon\n☒ CDS\n☐ EukaryoticGene (?)\n☐ ProkaryoticGene (?)\n☐ Intron\n☐ Exon\n☐ GFF –\\> See other packages\n☐ FASTX –\\> See I/O in other packages","category":"page"},{"location":"roadmap/#Compatibilities","page":"Roadmap","title":"Compatibilities","text":"","category":"section"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"Must interact with or extend:","category":"page"},{"location":"roadmap/","page":"Roadmap","title":"Roadmap","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl\nGraphs.jl","category":"page"},{"location":"simplefinder/#Finding-complete-and-overlapped-ORFs","page":"Finding ORFs","title":"Finding complete and overlapped ORFs","text":"","category":"section"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"The first implemented function is findorfs a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"using BioSequences, GeneFinder\n\n# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)\nseq = dna\"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC\"","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"Now lest us find the ORFs","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"findorfs(seq)\n\n12-element Vector{ORF}:\n ORF(29:40, '+', 2, 0.0)\n ORF(137:145, '+', 2, 0.0)\n ORF(164:184, '+', 2, 0.0)\n ORF(173:184, '+', 2, 0.0)\n ORF(236:241, '+', 2, 0.0)\n ORF(248:268, '+', 2, 0.0)\n ORF(362:373, '+', 2, 0.0)\n ORF(470:496, '+', 2, 0.0)\n ORF(551:574, '+', 2, 0.0)\n ORF(569:574, '+', 2, 0.0)\n ORF(581:601, '+', 2, 0.0)\n ORF(695:706, '+', 2, 0.0)","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"Two other functions (get_orfs_dna and get_orfs_aa) are implemented to get the ORFs in DNA and amino acid sequences, respectively. They use the findorfs function to first get the ORFs and then get the correspondance array of BioSequence objects.","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"get_orfs_dna(seq)\n\n12-element Vector{LongSubSeq{DNAAlphabet{4}}}:\n ATGCAACCCTGA\n ATGCGCTGA\n ATGCGTCGAATGGCACGGTGA\n ATGGCACGGTGA\n ATGTGA\n ATGTGTCCAACGGCAGTCTGA\n ATGCAACCCTGA\n ATGCACTGGCTGGTCCTGTCAATCTGA\n ATGTCACCGCACAAGGCAATGTGA\n ATGTGA\n ATGTGTCCAACGGCAGCCTGA\n ATGCAACCCTGA","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"get_orfs_aa(seq)\n\n12-element Vector{LongSubSeq{AminoAcidAlphabet}}:\n MQP*\n MR*\n MRRMAR*\n MAR*\n M*\n MCPTAV*\n MQP*\n MHWLVLSI*\n MSPHKAM*\n M*\n MCPTAA*\n MQP*","category":"page"},{"location":"simplefinder/#The-ORF-type","page":"Finding ORFs","title":"The ORF type","text":"","category":"section"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"For convenience, the ORF type is more stringent in preventing the creation of incompatible instances. As a result, attempting to create an instance with incompatible parameters will result in an error. For instance, the following code snippet will trigger an error:","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"ORF(1:10, '+', 4)\n\nERROR: AssertionError: Invalid frame value. Frame must be 1, 2, or 3.\nStacktrace:\n [1] ORF(location::UnitRange{Int64}, strand::Char, frame::Int64)\n   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:47\n [2] top-level scope\n   @ REPL[25]:1","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"Similar behavior will be encountered when the strand is neither + nor -. This precautionary measure helps prevent the creation of invalid ORFs, ensuring greater stability and enabling the extension of its interface. For example, after creating a specific ORF, users can seamlessly iterate over a sequence of interest and verify whether the ORF is contained within the sequence.","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"orf = ORF(137:145, '+', 2)\nseq[orf]\n\n9nt DNA Sequence:\nATGCGCTGA","category":"page"},{"location":"simplefinder/","page":"Finding ORFs","title":"Finding ORFs","text":"warning: Warning\nIt is still possible to create an ORF and pass it to a sequence that does not necessarily contain an actual open reading frame. This will be addressed in future versions of the package. But the benefit of having it is that it will retrieve the corresponding subsequence of the sequence in a convinient way (5' to 3') regardless of the strand.","category":"page"},{"location":"","page":"Home","title":"Home","text":"\n<p align=\"center\">\n  <img src=\"assets/logo.svg\" height=\"150\"><br/>\n  <i>A Gene Finder framework for Julia.</i>\n</p>","category":"page"},{"location":"","page":"Home","title":"Home","text":"\n<div style=\"text-align: center;\">\n\n<a href=\"https://camilogarciabotero.github.io/GeneFinder.jl/dev/\">\n  <img src=\"https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white\" alt=\"Documentation\">\n</a>\n<a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/releases/latest\">\n  <img src=\"https://img.shields.io/github/release/camilogarciabotero/GeneFinder.jl.svg\" alt=\"Release\">\n</a>\n<a href=\"https://doi.org/10.5281/zenodo.7519184\">\n  <img src=\"https://zenodo.org/badge/DOI/10.5281/zenodo.7519184.svg\" alt=\"DOI\">\n</a>\n<a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml\">\n<br>\n  <img src=\"https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg\" alt=\"GitHub Actions\">\n</a>\n<a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE\">\n  <img src=\"https://img.shields.io/badge/license-MIT-green.svg\" alt=\"License\">\n</a>\n<a href=\"https://www.repostatus.org/#wip\">\n  <img src=\"https://www.repostatus.org/badges/latest/wip.svg\" alt=\"Repo Status\">\n</a>\n<a href=\"https://pkgs.genieframework.com?packages=GeneFinder\">\n  <img src=\"https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/GeneFinder&label=downloads\" alt=\"Downloads\">\n</a>\n<a href=\"https://github.com/JuliaTesting/Aqua.jl\">\n  <img src=\"https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg\" alt=\"Aqua QA\">\n</a>\n\n</div>\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add GeneFinder\n","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"@misc{GeneFinder.jl,\n\tauthor  = {Camilo García},\n\ttitle   = {GeneFinder.jl},\n\turl     = {https://github.com/camilogarciabotero/GeneFinder.jl},\n\tversion = {v0.2.0},\n\tyear    = {2022},\n\tmonth   = {11}\n}","category":"page"}]
}
