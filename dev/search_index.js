var documenterSearchIndex = {"docs":
[{"location":"_index/","page":"-","title":"-","text":"<p align=\"center\"> <img src=\"../assets/logo.svg\" height=\"150\"><br/> <i>A Gene Finder framework for Julia.</i><br/><br/> <a href=\"https://www.repostatus.org/#wip\"> <img src=\"https://www.repostatus.org/badges/latest/wip.svg\"> </a> <a href=\"https://codecov.io/gh/camilogarciabotero/GeneFinder.jl\"> <img src=\"https://img.shields.io/codecov/c/github/camilogarciabotero/GeneFinder.jl?logo=codecov&logoColor=white\"> </a> <a href=\"https://camilogarciabotero.github.io/GeneFinder.jl/dev/\"> <img src=\"https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white\"> </a> <a href=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl\"> <img src=\"https://travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main\"> <a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE\"> <img src=\"https://img.shields.io/badge/license-MIT-green.svg\"> </a> </p>","category":"page"},{"location":"_index/","page":"-","title":"-","text":"","category":"page"},{"location":"_index/#Overview","page":"-","title":"Overview","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"_index/","page":"-","title":"-","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"_index/#Installation","page":"-","title":"Installation","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"_index/","page":"-","title":"-","text":"add GeneFinder","category":"page"},{"location":"_index/","page":"-","title":"-","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"_index/#Algorithms","page":"-","title":"Algorithms","text":"","category":"section"},{"location":"_index/#Coding-genes-(CDS-ORFs)","page":"-","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"☒ Simple   finder\n☐ EasyGene\n☐ GLIMMER\n☐ Prodigal - Pyrodigal\n☐ PHANOTATE\n☐ k-mer based gene finders (?)\n☐ Augustus (?)","category":"page"},{"location":"_index/#Non-coding-genes-(RNA)","page":"-","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"☐ Infernal\n☐ tRNAscan","category":"page"},{"location":"_index/#Other-features","page":"-","title":"Other features","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"☐ parallelism SIMD ?\n☐ memory management (?)\n☐ specialized types\n☒ Gene\n☒ ORF\n☒ CDS\n☐ EukaryoticGene (?)\n☐ ProkaryoticGene (?)\n☐ Codon\n☐ Intron\n☐ Exon\n☐ GFF –\\> See other packages\n☐ FASTX –\\> See I/O in other packages","category":"page"},{"location":"_index/#Compatibilities","page":"-","title":"Compatibilities","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"Must interact with or extend:","category":"page"},{"location":"_index/","page":"-","title":"-","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl","category":"page"},{"location":"_index/#Contributing","page":"-","title":"Contributing","text":"","category":"section"},{"location":"_index/#Citing","page":"-","title":"Citing","text":"","category":"section"},{"location":"_index/","page":"-","title":"-","text":"See CITATION.bib for the relevant reference(s).","category":"page"},{"location":"_index/","page":"-","title":"-","text":"Logo: gene analysis by Vector Points from the Noun Project","category":"page"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = GeneFinder\nDocTestSetup = quote\n    using GeneFinder\nend","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [GeneFinder]","category":"page"},{"location":"api/#GeneFinder.Codon","page":"API","title":"GeneFinder.Codon","text":"Codon <: BioSequence{DNAAlphabet{2}}\n\nA Struct representing a codon, which is a subtype of BioSequence with an Alphabet of type DNAAlphabet{2}. It has a single field x of type UInt8. This was implemente in The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.ORF","page":"API","title":"GeneFinder.ORF","text":"struct ORF\n    location::UnitRange{Int64}\n    strand::Char\n\n    ORF(location, strand) = new(location, strand)\nend\n\nThe ORF struct represents an open reading frame in a DNA sequence. It has two fields: \n\nlocation: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence\nstrand:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.cdsgenerator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.cdsgenerator","text":"cdsgenerator(sequence::LongDNA)\n\nA function to generete CDSs sequence out of a DNA sequence.\n\nThe cdsgenerator is a generator function that takes a LongDNA sequence and returns an iterator over the given sequence,     containing the coding sequences (CDSs) found in the sequence and the ORF.      It uses the simplefinder function to find open reading frames (ORFs) in the sequence,      and then it extracts the actual CDS sequence from each ORF.      The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.locationgenerator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.locationgenerator","text":"locationgenerator(sequence::LongDNA)\n\nGenerate the locations of ORFs in the given DNA sequence.\n\nThis function searches the sequence for start codons, and generates ranges of indices corresponding to the locations of ORFs in the sequence. The ORFs are generated by iterating over the start codon indices and searching for the first stop codon that follows each start codon. ORFs that contain premature stop codons are filtered out using the hasprematurestop function.\n\nThe sequence argument must be a LongDNA object, which is a type of DNA sequence with a longer maximum length than the DNA type.\n\nReturns:     A generator expression that yields ranges of indices corresponding to the locations of ORFs in the sequence.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.orfgenerator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.orfgenerator","text":"orfgenerator(sequence::LongDNA)\n\nGenerate ORFs from the given DNA sequence.\n\nThis function generates ORFs from the forward and reverse complement strands of the sequence using the locationgenerator function. It generates an ORF object for each range of indices returned by locationgenerator, and includes a '+' or '-' strand label to indicate the strand from which the ORF was generated.\n\nThe sequence argument must be a LongDNA object, which is a type of DNA sequence with a longer maximum length than the DNA type.\n\nReturns:     A generator expression that yields ORF objects corresponding to the ORFs in the sequence.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.proteingenerator-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.proteingenerator","text":"proteingenerator(sequence::LongDNA)\n\nAs its name suggest this generator function that iterates over the sequence to find proteins directly from a DNA sequence.      The cdsgenerator function takes a LongDNA sequence and returns a Vector{CDS} containing the      coding sequences (CDSs) found in the sequence. \n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.simplefinder-Tuple{BioSequences.LongDNA}","page":"API","title":"GeneFinder.simplefinder","text":"simplefinder(sequence::LongDNA)\n\nThe simplest algorithm that finds ORFs in a DNA sequence.\n\nThe simplefinder function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.     This function has not ORFs size and overlapping condition contraints. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"simplefinder/#A-simple-algorithm","page":"A first algorithm","title":"A simple algorithm","text":"","category":"section"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"The first implemented function is simplefinder a very non-restrictive ORF finder function that will catch all ORFs in a dedicated structure. Note that this will catch random ORFs not necesarily genes since it has no ORFs size or overlapping condition contraints. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"using BioSequences, GeneFinder\n\n# > 180195.SAMN03785337.LFLS01000089 -> finds only 1 gene in Prodigal (from Pyrodigal tests)\nseq = dna\"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC\"","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"726nt DNA Sequence:\nAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTG…GCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"simplefinder(seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"12-element Vector{ORF}:\n ORF(29:40, '+')\n ORF(137:145, '+')\n ORF(164:184, '+')\n ORF(173:184, '+')\n ORF(236:241, '+')\n ORF(248:268, '+')\n ORF(362:373, '+')\n ORF(470:496, '+')\n ORF(551:574, '+')\n ORF(569:574, '+')\n ORF(581:601, '+')\n ORF(695:706, '+')","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"Two other functions (findcds and findproteins) pass the sequence to simplefinder take the ORFs to index search the CDS and traslate into Protein:","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"findcds(seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"12-element Vector{CDS}:\n CDS(ORF(29:40, '+'), ATGCAACCCTGA)\n CDS(ORF(137:145, '+'), ATGCGCTGA)\n CDS(ORF(164:184, '+'), ATGCGTCGAATGGCACGGTGA)\n CDS(ORF(173:184, '+'), ATGGCACGGTGA)\n CDS(ORF(236:241, '+'), ATGTGA)\n CDS(ORF(248:268, '+'), ATGTGTCCAACGGCAGTCTGA)\n CDS(ORF(362:373, '+'), ATGCAACCCTGA)\n CDS(ORF(470:496, '+'), ATGCACTGGCTGGTCCTGTCAATCTGA)\n CDS(ORF(551:574, '+'), ATGTCACCGCACAAGGCAATGTGA)\n CDS(ORF(569:574, '+'), ATGTGA)\n CDS(ORF(581:601, '+'), ATGTGTCCAACGGCAGCCTGA)\n CDS(ORF(695:706, '+'), ATGCAACCCTGA)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"findproteins(seq)","category":"page"},{"location":"simplefinder/","page":"A first algorithm","title":"A first algorithm","text":"12-element Vector{Protein}:\n Protein(ORF(29:40, '+'), MQP*)\n Protein(ORF(137:145, '+'), MR*)\n Protein(ORF(164:184, '+'), MRRMAR*)\n Protein(ORF(173:184, '+'), MAR*)\n Protein(ORF(236:241, '+'), M*)\n Protein(ORF(248:268, '+'), MCPTAV*)\n Protein(ORF(362:373, '+'), MQP*)\n Protein(ORF(470:496, '+'), MHWLVLSI*)\n Protein(ORF(551:574, '+'), MSPHKAM*)\n Protein(ORF(569:574, '+'), M*)\n Protein(ORF(581:601, '+'), MCPTAA*)\n Protein(ORF(695:706, '+'), MQP*)","category":"page"},{"location":"","page":"Home","title":"Home","text":"<p align=\"center\">\n<img src=\"docs/assets/logo.svg\" height=\"150\"><br/> <i>A Gene Finder\nframework for Julia.</i><br/><br/>\n<a href=\"https://www.repostatus.org/#wip\">\n<img src=\"https://www.repostatus.org/badges/latest/wip.svg\"> </a>\n<a href=\"https://camilogarciabotero.github.io/GeneFinder.jl/dev/\">\n<img src=\"https://img.shields.io/badge/documentation-online-blue.svg?logo=Julia&logoColor=white\">\n</a>\n<a href=\"https://app.travis-ci.com/camilogarciabotero/GeneFinder.jl\">\n<img src=\"https://app.travis-ci.com/camilogarciabotero/GeneFinder.jl.svg?branch=main\">\n<a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml\">\n<img src=\"https://github.com/camilogarciabotero/GeneFinder.jl/actions/workflows/CI.yml/badge.svg\">\n<a href=\"https://github.com/camilogarciabotero/GeneFinder.jl/blob/main/LICENSE\">\n<img src=\"https://img.shields.io/badge/license-MIT-green.svg\"> </a> </a>\n</p>","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The main goal is to create a versatile module that enables apply different implemented algorithm to DNA sequences. See, for instance, BioAlignment implementations of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install GeneFinder from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add GeneFinder","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Algorithms","page":"Home","title":"Algorithms","text":"","category":"section"},{"location":"#Coding-genes-(CDS-ORFs)","page":"Home","title":"Coding genes (CDS - ORFs)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"☒ Simple   finder\n☐ EasyGene\n☐ GLIMMER\n☐ Prodigal - Pyrodigal\n☐ PHANOTATE\n☐ k-mer based gene finders (?)\n☐ Augustus (?)","category":"page"},{"location":"#Non-coding-genes-(RNA)","page":"Home","title":"Non-coding genes (RNA)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"☐ Infernal\n☐ tRNAscan","category":"page"},{"location":"#Other-features","page":"Home","title":"Other features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"☐ parallelism SIMD ?\n☐ memory management (?)\n☐ specialized types\n☒ Gene\n☒ ORF\n☒ Codon\n☒ CDS\n☐ EukaryoticGene (?)\n☐ ProkaryoticGene (?)\n☐ Codon\n☐ Intron\n☐ Exon\n☐ GFF –\\> See other packages\n☐ FASTX –\\> See I/O in other packages","category":"page"},{"location":"#Compatibilities","page":"Home","title":"Compatibilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Must interact with or extend:","category":"page"},{"location":"","page":"Home","title":"Home","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See CITATION.bib for the relevant reference(s).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Logo: gene analysis by Vector Points from the Noun Project","category":"page"}]
}
