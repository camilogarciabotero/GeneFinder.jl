var documenterSearchIndex = {"docs":
[{"location":"api/","page":"-","title":"-","text":"CurrentModule = GeneFinder\nDocTestSetup = quote\n    using GeneFinder\nend","category":"page"},{"location":"api/","page":"-","title":"-","text":"Modules = [GeneFinder]","category":"page"},{"location":"api/#GeneFinder.CDS","page":"-","title":"GeneFinder.CDS","text":"struct CDS\n    location::UnitRange{Int64}\n    strand::Char\n    sequence::LongDNA\nend\n\nThe CDS struct represents a coding sequence in a DNA sequence. It has three fields:\n\n- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) displaying the location of the ORF and the associate strand: forward ('+') or reverse ('-')\n- `sequence`: a `LongDNA` sequence representing the actual sequence of the CDS\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.ORF","page":"-","title":"GeneFinder.ORF","text":"struct ORF\n    location::UnitRange{Int64}\n    strand::Char\nend\n\nThe ORF struct represents an open reading frame in a DNA sequence. It has two fields: \n\n- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence\n- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.Protein","page":"-","title":"GeneFinder.Protein","text":"struct Protein\n    location::UnitRange{Int64}\n    strand::Char\n    sequence::LongDNA\nend\n\nSimilarly to the CDS struct, the Protein struct represents a encoded protein sequence in a DNA sequence.      It has three fields:\n\n- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) of the sequence\n- `sequence`: a `LongAA` sequence representing the actual translated sequence of the CDS\n\n\n\n\n\n","category":"type"},{"location":"api/#GeneFinder.findcds-Tuple{BioSequences.LongDNA}","page":"-","title":"GeneFinder.findcds","text":"`findcds(sequence::LongDNA)`\n\nA function to generete CDSs sequence out of a DNA sequence.\n\nThe findcds function takes a LongDNA sequence and returns a Vector{CDS}      containing the coding sequences (CDSs) found in the sequence.      It uses the simplefinder function to find open reading frames (ORFs) in the sequence,      and then it extracts the actual CDS sequence from each ORF.      The function also searches the reverse complement of the sequence, so it finds CDSs on both strands.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.findproteins-Tuple{BioSequences.LongDNA}","page":"-","title":"GeneFinder.findproteins","text":"`findproteins(sequence::LongDNA)`\n\nAs its name suggest this function generate the possible proteins directly from a DNA sequence.      The findcds function takes a LongDNA sequence and returns a Vector{CDS} containing the      coding sequences (CDSs) found in the sequence. \n\n\n\n\n\n","category":"method"},{"location":"api/#GeneFinder.simplefinder-Tuple{BioSequences.LongDNA}","page":"-","title":"GeneFinder.simplefinder","text":"`simplefinder(sequence::LongDNA)`\n\nThe simplest algorithm that finds ORFs in a DNA sequence.\n\nThe simplefinder function takes a LongDNA sequence and returns a Vector{ORF} containing the ORFs found in the sequence. It searches the sequence for start codons (ATG) and stops when it finds a stop codon (TAG, TAA, or TGA), adding each ORF it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFs on both strands.     This function has not ORFs size and overlapping condition contraints. Thus it might consider aa\"M*\" a posible encoding protein from the resulting ORFs.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"<!– # <img src=\"../assets/logo.svg\" width=\"30%\" align=\"right\" /> GeneFinder –>","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: MIT license)","category":"page"},{"location":"","page":"Home","title":"Home","text":"<!– (Image: Build Status) –> <!– (Image: Dev) –>","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"<!– (Image: CI) –>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<!– (Image: Aqua QA) –>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<!– (Image: Unit tests) –>","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a species-agnostic, algorithm extensible, sequence-anonymous (genome, metagenomes) gene finder library for the Julia Language.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The main idea is to create versatile module that enables apply different implemented algorithm to DNA sequences. See the BioAlignment implementation of different sequence alignment algorithms (local, global, edit-distance).","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install BioSequences from the julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add GeneFinder","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Algorithms","page":"Home","title":"Algorithms","text":"","category":"section"},{"location":"#Compatibilities","page":"Home","title":"Compatibilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Must interact with or extend:","category":"page"},{"location":"","page":"Home","title":"Home","text":"GenomicAnnotations.jl\nBioSequences.jl\nSequenceVariation.jl\nGenomicFeatures.jl\nFASTX.jl\nKmers.jl","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See CITATION.bib for the relevant reference(s).","category":"page"}]
}
