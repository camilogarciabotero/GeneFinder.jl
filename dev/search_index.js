var documenterSearchIndex = {"docs":
[{"location":"_index/","page":"-","title":"-","text":"<p align=\"center\"> <img src=\"../assets/logo.svg\" height=\"150\"><br/> <i>A Gene Finder framework for Julia.</i><br/><br/> <a href=\"https://www.repostatus.org/#wip\"> <img src=\"https://www.repostatus.org/badges/latest/wip.svg\"> </a> <a href=\"https://codecov.io/gh/camilogarciabotero/GeneFinder.jl\"> <img src=\"https://img.shields.io/codecov/c/github/camilogarciabotero/GeneFinder.jl?logo=codecov&logoColor=white\"> </a> <a href=\"https://camilogarciabotero.github.io/GeneFinder.jl/dev/\"> <img src=\"https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white\"> </a> <a href=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl\"> <img src=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main\"> <a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE\"> <img src=\"https://img.shields.io/badge/license-MIT-green.svg\"> </a> </p>","category":"page"},{"location":"_index/","page":"-","title":"-","text":"","category":"page"},{"location":"_index/#Overview","page":"-","title":"Overview","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"_index/","page":"-","title":"-","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"_index/#Installation","page":"-","title":"Installation","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"_index/","page":"-","title":"-","text":"add GeneFinder","category":"page"},{"location":"_index/","page":"-","title":"-","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"_index/#Algorithms","page":"-","title":"Algorithms","text":"","category":"section"},{"location":"_index/#Coding-genes-(CDS-ORFs)","page":"-","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"☒ Simple   finder\n☐ EasyGene\n☐ GLIMMER\n☐ Prodigal - Pyrodigal\n☐ PHANOTATE\n☐ k-mer based gene finders (?)\n☐ Augustus (?)","category":"page"},{"location":"_index/#Non-coding-genes-(RNA)","page":"-","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"☐ Infernal\n☐ tRNAscan","category":"page"},{"location":"_index/#Other-features","page":"-","title":"Other features","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"☐ parallelism SIMD ?\n☐ memory management (?)\n☐ specialized types\n☒ Gene\n☒ ORF\n☒ CDS\n☐ EukaryoticGene (?)\n☐ ProkaryoticGene (?)\n☐ Codon\n☐ Intron\n☐ Exon\n☐ GFF –\\> See other packages\n☐ FASTX –\\> See I/O in other packages","category":"page"},{"location":"_index/#Compatibilities","page":"-","title":"Compatibilities","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"Must interact with or extend:","category":"page"},{"location":"_index/","page":"-","title":"-","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl","category":"page"},{"location":"_index/#Contributing","page":"-","title":"Contributing","text":"","category":"section"},{"location":"_index/#Citing","page":"-","title":"Citing","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"See CITATION.bib for the relevant reference(s).","category":"page"},{"location":"_index/","page":"-","title":"-","text":"Logo: gene analysis by Vector Points from the Noun Project","category":"page"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = GeneFinder\nDocTestSetup = quote\n    using GeneFinder\nend","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [GeneFinder]","category":"page"},{"location":"api/#GeneFinder.CDS","page":"API","title":"GeneFinder.CDS","text":"struct CDS\n    orf::ORF\n    sequence::LongSubSeq{DNAAlphabet{4}}\nend\n\nThe CDS struct represents a coding sequence in a DNA sequence. It has three fields:\n\norf: is the basic composible type (location::UnitRange{Int}, strand::Char) displaying the location of the ORF and the associate strand: forward ('+') or reverse ('-')\nsequence: a LongDNA sequence representing the actual sequence of the CDS\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.GeneFeatures","page":"API","title":"GeneFinder.GeneFeatures","text":"struct GeneFeatures\n    seqname::String\n    start::Int64\n    stop::Int64\n    score::Float64\n    strand::Char\n    frame::'.'\n    attribute::\nend\n\nThis is the main Gene struct, based on the fields that could be found in a GFF3, still needs to be defined correctly,     The idea is correct the frame and attributes that will have something like a possible list (id=Char;name=;locus_tag).     The write and get functions should have a dedicated method for this struct.\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.ORF","page":"API","title":"GeneFinder.ORF","text":"struct ORF\n    location::UnitRange{Int64}\n    strand::Char\nend\n\nThe ORF struct represents an open reading frame in a DNA sequence. It has two fields: \n\nlocation: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence\nstrand:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.Protein","page":"API","title":"GeneFinder.Protein","text":"struct Protein\n    sequence::LongSequence\n    orf::ORF\nend\n\nSimilarly to the CDS struct, the Protein struct represents a encoded protein sequence in a DNA sequence.      It has three fields:\n\norf: is the basic composible type (location::UnitRange{Int}, strand::Char) of the sequence\nsequence: a LongSequence sequence representing the actual translated sequence of the CDS\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.cdsgenerator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.cdsgenerator","text":"cdsgenerator(sequence::LongDNA; kwargs...)\ncdsgenerator(sequence::String; kwargs...)\n\nA function to generete CDSs sequence out of a DNA sequence.\n\nThe cdsgenerator is a generator function that takes a LongDNA sequence and returns an iterator over the given sequence,     containing the coding sequences (CDSs) found in the sequence and the ORF.      It uses the findorfs function to find open reading frames (ORFs) in the sequence,      and then it extracts the actual CDS sequence from each ORF, returining both.      The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.eachcodon-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.eachcodon","text":"eachcodon(sequence::LongDNA)\n\nIterate through the codons in the sequence of type LongDNA.\n\nReturns an iterator yielding Codon objects for each codon in the sequence.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.fasta_to_dna-Tuple{String}","page":"API","title":"GeneFinder.fasta_to_dna","text":"fasta_to_dna(input::String)\n\nConverts a FASTA formatted file (even if it is a multi-fasta) to an array of LongDNA objects.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.findorfs-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.findorfs","text":"findorfs(sequence::LongDNA; kwargs...)::Vector{ORF}\nfindorfs(sequence::String; kwargs...)::Vector{ORF}\n\nA simple implementation that finds ORFs in a DNA sequence.\n\nThe findorfs function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence.      It searches entire regularly expressed CDS, adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.         Extending the starting codons with the alternative_start = true will search for ATG, GTG, and TTG.     Some studies have shown that in E. coli (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.\n\nnote: Note\nThis function has not ORFs scoring scheme. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.getcds-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.getcds","text":"getcds(input::LongDNA, output::String; kwargs...)\ngetcds(input::String, output::String; kwargs...)\n\nThis function will take a LongDNA or String sequence and by means of the findorfs() function will push LongSubSeq{DNAAlphabet{4}} into a Vector{}\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.getproteins-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.getproteins","text":"getproteins(input::LongDNA, output::String; kwargs...)\ngetproteins(input::String, output::String; kwargs...)\n\nSimilar to getcds() function, it will take a LongDNA or String sequence and by means of the findorfs() and the translate() function will push LongAAs into a Vector\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.hasprematurestop-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.hasprematurestop","text":"hasprematurestop(sequence::LongDNA)::Bool\n\nDetermine whether the sequence of type LongDNA contains a premature stop codon.\n\nReturns a boolean indicating whether the sequence has more than one stop codon.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.locationiterator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.locationiterator","text":"locationiterator(sequence::LongDNA; alternative_start::Bool=false)\n\nThis is an iterator function that uses regular expressions to search the entire CDS (instead of start and stop codons) in a LongDNA sequence.     It uses an anonymous function that will find the first regularly expressed CDS. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.proteingenerator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.proteingenerator","text":"proteingenerator(sequence::LongDNA; kwargs...)\nproteingenerator(sequence::String; kwargs...)\n\nAs its name suggest this generator function iterates over the sequence to find proteins directly from a DNA sequence.      The proteingenerator function takes a LongDNA sequence and returns a Vector{CDS} containing the      coding sequences (CDSs) found in the sequence and the associated ORF.\n\nKeywords\n\ncode::GeneticCode=BioSequences.standard_genetic_code: The genetic code by which codons will be translated. See BioSequences.ncbi_trans_table for more info. \nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_bed-Tuple{BioSequences.LongDNA, String}","page":"API","title":"GeneFinder.write_bed","text":"write_bed(input::LongDNA, output::String; kwargs...)\nwrite_bed(input::String, output::String; kwargs...)\n\nWrite BED data to a file.\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_cds-Tuple{BioSequences.LongDNA, String}","page":"API","title":"GeneFinder.write_cds","text":"write_cds(input::LongDNA, output::String; kwargs...)\nwrite_cds(input::String, output::String; kwargs...)\n\nWrite a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.\n\nKeywords\n\nalternative_start: A boolean value indicating whether alternative start codons should be used when identifying CDSs. Default is false.\nmin_len: An integer representing the minimum length that a CDS must have in order to be included in the output file. Default is 6.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_gff-Tuple{String, String}","page":"API","title":"GeneFinder.write_gff","text":"write_gff(input::LongDNA, output::String; kwargs...)\nwrite_gff(input::String, output::String; kwargs...)\n\nWrite GFF data to a file.\n\nKeywords\n\nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.write_proteins-Tuple{BioSequences.LongDNA, String}","page":"API","title":"GeneFinder.write_proteins","text":"write_proteins(input::LongDNA, output::String; kwargs...)\n\nWrite the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.\n\nKeywords\n\ncode::GeneticCode=BioSequences.standard_genetic_code: The genetic code by which codons will be translated. See BioSequences.ncbi_trans_table for more info. \nalternative_start::Bool=false: If true will pass the extended start codons to search. This will increase 3x the exec. time.\nmin_len::Int64=6:  Length of the allowed ORF. Default value allow aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"oldindex/","page":"-","title":"-","text":"","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"engine: knitr cache: true –-","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"<p align=\"center\">   <img src=\"../assets/logo.svg\" height=\"150\"><br/>   <i>A Gene Finder framework for Julia.</i><br/><br/>   <a href=\"https://www.repostatus.org/#wip\">     <img src=\"https://www.repostatus.org/badges/latest/wip.svg\">   </a>   <a href=\"https://codecov.io/gh/camilogarciabotero/GeneFinder.jl\">     <img src=\"https://img.shields.io/codecov/c/github/camilogarciabotero/GeneFinder.jl?logo=codecov&logoColor=white\">   </a>   <a href=\"https://camilogarciabotero.github.io/GeneFinder.jl/dev/\">     <img src=\"https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white\">   </a>   <a href=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl\">     <img src=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main\">   <a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE\">     <img src=\"https://img.shields.io/badge/license-MIT-green.svg\">   </a> </p>","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"","category":"page"},{"location":"oldindex/#Overview","page":"-","title":"Overview","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"oldindex/#Installation","page":"-","title":"Installation","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"add GeneFinder","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"oldindex/#Algorithms","page":"-","title":"Algorithms","text":"","category":"section"},{"location":"oldindex/#Coding-genes-(CDS-ORFs)","page":"-","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"[x] Simple finder\n[ ] EasyGene\n[ ] GLIMMER\n[ ] Prodigal - Pyrodigal\n[ ] PHANOTATE\n[ ] k-mer based gene finders (?)\n[ ] Augustus (?)","category":"page"},{"location":"oldindex/#Non-coding-genes-(RNA)","page":"-","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"[ ] Infernal\n[ ] tRNAscan","category":"page"},{"location":"oldindex/#Other-features","page":"-","title":"Other features","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"[ ] parallelism SIMD ?\n[ ] memory management (?)\n[ ] specialized types\n[x] Gene\n[x] ORF\n[x] CDS\n[ ] EukaryoticGene (?)\n[ ] ProkaryoticGene (?)\n[ ] Codon\n[ ] Intron\n[ ] Exon\n[ ] GFF –> See other packages\n[ ] FASTX –> See I/O in other packages","category":"page"},{"location":"oldindex/#Compatibilities","page":"-","title":"Compatibilities","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"Must interact with or extend:","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl","category":"page"},{"location":"oldindex/#Contributing","page":"-","title":"Contributing","text":"","category":"section"},{"location":"oldindex/#Citing","page":"-","title":"Citing","text":"","category":"section"},{"location":"oldindex/","page":"-","title":"-","text":"See CITATION.bib for the relevant reference(s).","category":"page"},{"location":"oldindex/","page":"-","title":"-","text":"Logo: gene analysis by Vector Points from the Noun Project","category":"page"},{"location":"simplefinder/#A-simple-algorithm","page":"A first algorithm","title":"A simple algorithm","text":"","category":"section"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"The first implemented function is findorfs a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"using BioSequences, GeneFinder\n\n# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)\nseq = dna\"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC\"","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"726nt DNA Sequence:\nAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC","category":"page"},{"location":"simplefinder/#Finding-all-ORFs","page":"A first algorithm","title":"Finding all ORFs","text":"","category":"section"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"findorfs(seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"12-element Vector{ORF}:\n ORF(29:40, '+')\n ORF(137:145, '+')\n ORF(164:184, '+')\n ORF(173:184, '+')\n ORF(236:241, '+')\n ORF(248:268, '+')\n ORF(362:373, '+')\n ORF(470:496, '+')\n ORF(551:574, '+')\n ORF(569:574, '+')\n ORF(581:601, '+')\n ORF(695:706, '+')","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"Two other functions (getcds and getproteins) pass the sequence to findorfs take the ORFs and act as generators of the sequence, so this way the can be collected in the REPL as an standard output or written into a file more conviniently using the FASTX IO system:","category":"page"},{"location":"simplefinder/#Generting-cds-and-proteins-with-its-ORF","page":"A first algorithm","title":"Generting cds and proteins with its ORF","text":"","category":"section"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"getcds(seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"12-element Vector{LongSequence{DNAAlphabet{4}}}:\n ATGCAACCCTGA\n ATGCGCTGA\n ATGCGTCGAATGGCACGGTGA\n ATGGCACGGTGA\n ATGTGA\n ATGTGTCCAACGGCAGTCTGA\n ATGCAACCCTGA\n ATGCACTGGCTGGTCCTGTCAATCTGA\n ATGTCACCGCACAAGGCAATGTGA\n ATGTGA\n ATGTGTCCAACGGCAGCCTGA\n ATGCAACCCTGA","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"getproteins(seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"12-element Vector{LongAA}:\n MQP*\n MR*\n MRRMAR*\n MAR*\n M*\n MCPTAV*\n MQP*\n MHWLVLSI*\n MSPHKAM*\n M*\n MCPTAA*\n MQP*","category":"page"},{"location":"simplefinder/#Combining-FASTX-for-reading-and-writing-a-fasta-record","page":"A first algorithm","title":"Combining FASTX for reading and writing a fasta record","text":"","category":"section"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"using FASTX\n\nwrite_proteins(\"../../test/data/NC_001884.fasta\", \"proteins.fasta\")","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"head proteins.fasta\n>location=75:113 strand=+\nMKLNLRIGVISN*\n>location=144:215 strand=+\nMLTITSFKTILNSSFFFSELDSM*\n>location=210:215 strand=+\nM*\n>location=237:374 strand=+\nMLFLTVLLSISDCVSCNPLSSFFAFWSSLNSSSNAAFLFKKSSSL*\n>location=337:402 strand=+\nMQLFSSKKVHHCKCHFHIYRR*","category":"page"},{"location":"simplefinder/#Writting-cds-and-proteins-fastas","page":"A first algorithm","title":"Writting cds and proteins fastas","text":"","category":"section"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"write_cds(\"cds.fasta\", seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"cat cds.fasta\n>location=29:40 strand=+\nATGCAACCCTGA\n>location=137:145 strand=+\nATGCGCTGA\n>location=164:184 strand=+\nATGCGTCGAATGGCACGGTGA\n>location=173:184 strand=+\nATGGCACGGTGA\n>location=236:241 strand=+\nATGTGA\n>location=248:268 strand=+\nATGTGTCCAACGGCAGTCTGA\n>location=362:373 strand=+\nATGCAACCCTGA\n>location=470:496 strand=+\nATGCACTGGCTGGTCCTGTCAATCTGA\n>location=551:574 strand=+\nATGTCACCGCACAAGGCAATGTGA\n>location=569:574 strand=+\nATGTGA\n>location=581:601 strand=+\nATGTGTCCAACGGCAGCCTGA\n>location=695:706 strand=+\nATGCAACCCTGA","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"write_proteins(\"proteins.fasta\", seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"cat proteins.fasta\n>location=29:40 strand=+\nMQP*\n>location=137:145 strand=+\nMR*\n>location=164:184 strand=+\nMRRMAR*\n>location=173:184 strand=+\nMAR*\n>location=236:241 strand=+\nM*\n>location=248:268 strand=+\nMCPTAV*\n>location=362:373 strand=+\nMQP*\n>location=470:496 strand=+\nMHWLVLSI*\n>location=551:574 strand=+\nMSPHKAM*\n>location=569:574 strand=+\nM*\n>location=581:601 strand=+\nMCPTAA*\n>location=695:706 strand=+\nMQP*","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add GeneFinder","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Algorithms","page":"Home","title":"Algorithms","text":"","category":"section"},{"location":"#Coding-genes-(CDS-ORFs)","page":"Home","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"☒ Simple finder\n☐ EasyGene\n☐ GLIMMER\n☐ Prodigal - Pyrodigal\n☐ PHANOTATE\n☐ k-mer based gene finders (?)\n☐ Augustus (?)","category":"page"},{"location":"#Non-coding-genes-(RNA)","page":"Home","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"☐ Infernal\n☐ tRNAscan","category":"page"},{"location":"#Other-features","page":"Home","title":"Other features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"☐ parallelism SIMD ?\n☐ memory management (?)\n☐ incorporate Ribosime Binding Sites (RBS)\n☐ incorporate Programmed Reading Frame Shifting (PRFS)\n☐ specialized types\n☒ Gene\n☒ ORF\n☒ Codon\n☒ CDS\n☐ EukaryoticGene (?)\n☐ ProkaryoticGene (?)\n☐ Intron\n☐ Exon\n☐ GFF –\\> See other packages\n☐ FASTX –\\> See I/O in other packages","category":"page"},{"location":"#Compatibilities","page":"Home","title":"Compatibilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Must interact with or extend:","category":"page"},{"location":"","page":"Home","title":"Home","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl\nGraphs..jl","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See CITATION.bib for the relevant reference(s).","category":"page"}]
}
