


## The Main ORF type {#The-Main-ORF-type}

The main type of the package is `ORFI` which represents an Open Reading Frame Interval. It is a subtype of the `GenomicInterval` type from the `GenomicFeatures` package.
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.GeneFinderMethod' href='#GeneFinder.GeneFinderMethod'><span class="jlbinding">GeneFinder.GeneFinderMethod</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type GeneFinderMethod
```


Abstract base type for different ORF finding methods/algorithms.

Subtypes should implement the calling interface to find ORFs in a sequence.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/types.jl#L5-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.OpenReadingFrame' href='#GeneFinder.OpenReadingFrame'><span class="jlbinding">GeneFinder.OpenReadingFrame</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct ORF{F}
```


The `ORF` struct represents an Open Reading Frame (ORF) in genomics.

**Fields**
- `seqid::Symbol`: The identifier of the sequence to which the ORF belongs.
  
- `range::UnitRange{<:Int64}`: The position range of the ORF on the sequence.
  
- `strand::Strand`: The strand on which the ORF is located.
  
- `frame::Int8`: The reading frame of the ORF (1, 2, or 3).
  
- `features::NamedTuple`: The features associated with the ORF.
  

**Example**

```julia
ORF{NaiveFinder}(:seq01, 1:33, PSTRAND, Int8(1), (;score = 0.8))
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/types.jl#L57-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.Strand' href='#GeneFinder.Strand'><span class="jlbinding">GeneFinder.Strand</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Strand
```


An enumeration type representing DNA strand orientation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/types.jl#L16-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.features-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.features-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.features</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
features(i::ORF{F})
```


Extracts the features from an `ORF` object.

**Arguments**
- `i::ORF{F}`: An `ORF` object.
  

**Returns**

The features of the `ORF` object.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/types.jl#L165-L175" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.sequence-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.sequence-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.sequence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sequence(i::ORF{F})
```


Extracts the DNA sequence corresponding to the given open reading frame (ORF). Uses the source sequence referenced by the ORF&#39;s seqid.

For positive strand ORFs, returns a LongSubSeq view (avoiding unnecessary copying). For negative strand ORFs, returns the reverse complement (requires allocation).

**Arguments**
- `i::ORF{F}`: The open reading frame (ORF) for which the DNA sequence needs to be extracted.
  

**Returns**
- A LongSubSeq for positive strand ORFs, or a reverse complement sequence for negative strand ORFs.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/types.jl#L99-L113" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.source-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.source-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.source</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
source(i::ORF{F})
```


Get the source sequence associated with the given `ORF` object.

**Arguments**
- `i::ORF{F}`: The `ORF` object for which to retrieve the source sequence.
  

**Returns**

The source sequence associated with the `ORF` object.

::: warning Warning

The `source` method works if the sequence is defined in the global scope. Otherwise it will throw an error.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/types.jl#L147-L160" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Finding ORFIs {#Finding-ORFIs}

The function `findorfs` serves as a method interface as it is generic method that can handle different gene finding methods.
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{F}, Tuple{N}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{F}, Tuple{N}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.findorfs</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; ::F, kwargs...) where {N, F<:GeneFinderMethod}
```


This is the main interface method for finding open reading frames (ORFs) in a DNA sequence.

It takes the following required arguments:
- `sequence`: The nucleic acid sequence to search for ORFs.
  
- `finder`: The algorithm used to find ORFs. It can be either `NaiveFinder`, `NaiveCollector` or yet other implementations.
  

**Keyword Arguments regardless of the finder method:**
- `alternative_start::Bool`: A boolean indicating whether to consider alternative start codons. Default is `false`.
  
- `minlen::Int`: The minimum length of an ORF. Default is `6`.
  
- `scheme::Function`: The scoring scheme to use for scoring the sequence from the ORF. Default is `nothing`.
  

**Returns**

A vector of `ORF` objects representing the found ORFs.

**Example**

```julia
sequence = randdnaseq(120)

120nt DNA Sequence:
 GCCGGACAGCGAAGGCTAATAAATGCCCGTGCCAGTATC…TCTGAGTTACTGTACACCCGAAAGACGTTGTACGCATTT

findorfs(sequence, finder=NaiveFinder)

1-element Vector{ORF}:
 ORF{NaiveFinder}(77:118, '-', 2)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/findorfs.jl#L3-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Finding ORFs using BioRegex {#Finding-ORFs-using-BioRegex}
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder.NaiveFinder</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
NaiveFinder(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) -> Vector{ORF{F}} where {N,F}
```


A simple implementation that finds ORFIs in a DNA sequence.

The `NaiveFinder` method takes a LongSequence{DNAAlphabet{4}} sequence and returns a Vector{ORFIs} containing the ORFIs found in the sequence.      It searches entire regularly expressed CDS, adding each ORFI it finds to the vector. The function also searches the reverse complement of the sequence, so it finds ORFIs on both strands.         Extending the starting codons with the `alternative_start = true` will search for ATG, GTG, and TTG.     Some studies have shown that in _E. coli_ (K-12 strain), ATG, GTG and TTG are used 83 %, 14 % and 3 % respectively.

::: tip Note

This function has neither ORFIs scoring scheme by default nor length constraints. Thus it might consider `aa"M*"` a posible encoding protein from the resulting ORFIs.

:::

**Required Arguments**
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The nucleic acid sequence to search for ORFIs.
  

**Keywords Arguments**
- `alternative_start::Bool`: If true will pass the extended start codons to search. This will increase 3x the execution time. Default is `false`.
  
- `minlen::Int64=6`:  Length of the allowed ORFI. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFIs.
  

. Default value allow `aa"M*"` a posible encoding protein from the resulting ORF

::: tip Note

As the scheme is generally a scoring function that at least requires a sequence, one simple scheme is the log-odds ratio score. This score is a log-odds ratio that compares the probability of the sequence generated by a coding model to the probability of the sequence generated by a non-coding model:

$$S(x) = \sum_{i=1}^{L} \beta_{x_{i}x} = \sum_{i=1} \log \frac{a^{\mathscr{m}_{1}}_{i-1} x_i}{a^{\mathscr{m}_{2}}_{i-1} x_i}$$

If the log-odds ratio exceeds a given threshold (`η`), the sequence is considered likely to be coding. See [`lordr`](@ref) for more information about coding creteria.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/algorithms/naivefinder.jl#L39-L68" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._estimate_orf_count-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder._estimate_orf_count-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder._estimate_orf_count</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}
```


Estimate the number of ORFs based on start codon count.

Counts ATG start codons in the sequence using a k-mer iterator. Returns the count as a heuristic estimate for pre-allocating the ORF vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/algorithms/naivefinder.jl#L104-L111" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._locationiterator-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder._locationiterator-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder._locationiterator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_locationiterator(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N}
```


This is an iterator function that uses regular expressions to search the entire ORFI (instead of start and stop codons) in a `LongSequence{DNAAlphabet{4}}` sequence.     It uses an anonymous function that will find the first regularly expressed ORFI. Then using this anonymous function it creates an iterator that will apply it until there is no other CDS.

::: tip Note

As a note of the implementation we want to expand on how the ORFIs are found:

The expression `(?:[N]{3})*?` serves as the boundary between the start and stop codons.  Within this expression, the character class `[N]{3}` captures exactly three occurrences of any character (representing nucleotides using IUPAC codes).  This portion functions as the regular codon matches.  Since it is enclosed within a non-capturing group `(?:)` and followed by `*?`, it allows for the matching of intermediate codons, but with a preference for the smallest number of repetitions. 

In summary, the regular expression `ATG(?:[N]{3})*?T(AG|AA|GA)` identifies patterns that start with &quot;ATG,&quot; followed by any number of three-character codons (represented by &quot;N&quot; in the IUPAC code), and ends with a stop codon &quot;TAG,&quot; &quot;TAA,&quot; or &quot;TGA.&quot; This pattern is commonly used to identify potential protein-coding regions within genetic sequences.

See more about the discussion [here](https://discourse.julialang.org/t/how-to-improve-a-generator-to-be-more-memory-efficient-when-it-is-collected/92932/8?u=camilogarciabotero)

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/algorithms/naivefinder.jl#L6-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._search_strand!-Union{Tuple{N}, Tuple{Vector{OpenReadingFrame{NaiveFinderLazy}}, Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Symbol, Strand, Int64, Bool, Int64}} where N' href='#GeneFinder._search_strand!-Union{Tuple{N}, Tuple{Vector{OpenReadingFrame{NaiveFinderLazy}}, Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Symbol, Strand, Int64, Bool, Int64}} where N'><span class="jlbinding">GeneFinder._search_strand!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_search_strand!(orfs::Vector{ORF{NaiveFinderLazy}}, seq::NucleicSeqOrView{DNAAlphabet{N}}, seqname::Symbol, strand::Strand, seqlen::Int, alternative_start::Bool, minlen::Int64) where {N}
```


Helper function to search for ORFs in a single strand direction.

Avoids code duplication between forward and reverse strand searching.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/algorithms/naivefinder.jl#L128-L134" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.NaiveCollector-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder.NaiveCollector-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder.NaiveCollector</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
NaiveCollector(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) -> Vector{ORF{F}} where {N,F}
```


The `NaiveCollector` function searches for open reading frames (ORFs) in a DNA sequence. It takes the following arguments:

**Required Arguments**
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The nucleic sequence to search for ORFs.
  

**Keywords Arguments**
- `alternative_start::Bool`: A flag indicating whether to consider alternative start codons. Default is `false`.
  
- `minlen::Int64`: The minimum length of an ORF. Default is `6`.
  
- `overlap::Bool`: A flag indicating whether to allow overlapping ORFs. Default is `false`.
  

The function returns a sorted vector of `ORF{NaiveCollector}` objects, representing the identified ORFs.

::: tip Note

This method finds, by default, non-overlapping ORFs in the given sequence. It is much faster than the `NaiveFinder` method.  Althought it uses the same regular expression to find ORFs in a source sequence,   it levarages on the `eachmatch` function to find all the ORFs in the sequence.

:::

::: warning Warning

Using the `overlap = true` flag will increase the runtime of the function significantly, but some of the ORFs found may display     premature stop codons.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/algorithms/naivecollector.jl#L5-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Writing ORFs to files {#Writing-ORFs-to-files}
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_bed-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_bed-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_bed</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
write_orfs_bed(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)
```


Write BED data to a file.

**Arguments**
- `input`: The input DNA sequence NucSeq or a view.
  
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
  
- `finder`: The algorithm used to find ORFIs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.
  

**Keywords**
- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
  
- `minlen::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/io.jl#L3-L17" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_faa-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_faa-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_faa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_faa(input::NucleicSeqOrView{DNAAlphabet{4}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
write_orfs_faa(input::NucleicSeqOrView{DNAAlphabet{4}}, output::String, finder::F; kwargs...)
```


Write the protein sequences encoded by the coding sequences (CDSs) of a given DNA sequence to the specified file.

**Arguments**
- `input`: The input DNA sequence NucSeq or a view.
  
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
  
- `finder`: The algorithm used to find ORFIs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.
  

**Keywords**
- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
  
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
  
- `minlen::Int64=6`:  Length of the allowed ORFI. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFIs.
  

**Examples**

```julia
filename = "output.faa"

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

open(filename, "w") do file
     write_orfs_faa(seq, file)
end
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/io.jl#L112-L140" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_fna-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_fna-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_fna</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
write_orfs_fna(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)
```


Write a file containing the coding sequences (CDSs) of a given DNA sequence to the specified file.

**Arguments**
- `input::NucleicAcidAlphabet{DNAAlphabet{N}}`: The input DNA sequence.
  
- `output::IO`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
  
- `finder::F`: The algorithm used to find ORFs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.
  

**Keywords**
- `alternative_start::Bool=false`: If true, alternative start codons will be used when identifying CDSs. Default is `false`.
  
- `minlen::Int64=6`: The minimum length that a CDS must have in order to be included in the output file. Default is `6`.
  

**Examples**

```julia
filename = "output.fna"

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

open(filename, "w") do file
     write_orfs_fna(seq, file, NaiveFinder())
end
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/io.jl#L46-L74" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_gff-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_gff-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_gff</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::Union{IOStream, IOBuffer}, finder::F; kwargs...)
write_orfs_gff(input::NucleicSeqOrView{DNAAlphabet{N}}, output::String, finder::F; kwargs...)
```


Write GFF data to a file.

**Arguments**
- `input`: The input DNA sequence NucSeq or a view.
  
- `output`: The otput format, it can be a file (`String`) or a buffer (`IOStream` or `IOBuffer)
  
- `finder`: The algorithm used to find ORFs. It can be either `NaiveFinder()` or `NaiveFinderScored()`.
  

**Keywords**
- `code::GeneticCode=BioSequences.standard_genetic_code`: The genetic code by which codons will be translated. See `BioSequences.ncbi_trans_table` for more info. 
  
- `alternative_start::Bool=false`: If true will pass the extended start codons to search. This will increase 3x the exec. time.
  
- `minlen::Int64=6`:  Length of the allowed ORF. Default value allow `aa"M*"` a posible encoding protein from the resulting ORFs.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/50478f27f5ac3d8dfd1db72b5bbd5aa2edae9ecc/src/io.jl#L180-L197" target="_blank" rel="noreferrer">source</a></Badge>

</details>

