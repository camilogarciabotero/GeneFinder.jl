


## Core Types {#Core-Types}

The main types of the package are `ORF` (Open Reading Frame) and `ORFCollection`.  An `ORF` stores coordinates and metadata, while `ORFCollection` bundles ORFs with  their source sequence for clean sequence extraction.
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.GeneFinderMethod' href='#GeneFinder.GeneFinderMethod'><span class="jlbinding">GeneFinder.GeneFinderMethod</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
abstract type GeneFinderMethod
```


Abstract base type for different ORF finding methods/algorithms.

Subtypes should implement the calling interface to find ORFs in a sequence, returning an `ORFCollection`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L6-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.ORFCollection' href='#GeneFinder.ORFCollection'><span class="jlbinding">GeneFinder.ORFCollection</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct ORFCollection{F<:GeneFinderMethod, S<:NucleicSeqOrView}
```


A collection of Open Reading Frames (ORFs) bundled with a view of their source sequence.

This is the primary return type for all `GeneFinderMethod` implementations. It provides a clean API for accessing ORFs and their corresponding sequences without relying on global state.

The source is always stored as a view (`LongSubSeq`) to avoid unnecessary copying while maintaining a reference to the original sequence data.

**Type Parameters**
- `F<:GeneFinderMethod`: The gene finding algorithm used to identify these ORFs.
  
- `S<:NucleicSeqOrView`: The type of the source sequence view.
  

**Fields**
- `source::S`: A view of the source DNA sequence containing the ORFs.
  
- `orfs::Vector{ORF{F}}`: The vector of ORFs found in the source sequence.
  

**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAG"
collection = findorfs(seq, finder=NaiveFinder)

# Source is stored as a view
typeof(source(collection))  # LongSubSeq{DNAAlphabet{4}}

# Iteration
for orf in collection
    println(orf)
end

# Sequence extraction
orfseq = sequence(collection, 1)
```


See also: [`ORF`](@ref), [`sequence`](/api#GeneFinder.sequence-Tuple{ORFCollection,%20Int64}), [`source`](/api#GeneFinder.source-Tuple{ORFCollection})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L126-L166" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.OpenReadingFrame' href='#GeneFinder.OpenReadingFrame'><span class="jlbinding">GeneFinder.OpenReadingFrame</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct OpenReadingFrame{F<:GeneFinderMethod}
```


The `OpenReadingFrame` (aliased as `ORF`) struct represents an Open Reading Frame in genomics.

An ORF is a lightweight coordinate-based structure that stores the location and metadata of a potential coding region. ORFs are typically contained within an `ORFCollection`, which provides access to the source sequence.

**Type Parameter**
- `F<:GeneFinderMethod`: The gene finding algorithm used to identify this ORF.
  

**Fields**
- `range::UnitRange{Int64}`: The position range (start:stop) of the ORF on the sequence.
  
- `strand::Strand`: The strand orientation (`PSTRAND` or `NSTRAND`).
  
- `frame::Int8`: The reading frame (1, 2, or 3).
  
- `features::NamedTuple`: Additional features/metadata associated with the ORF.
  

**Example**

```julia
using BioSequences, GeneFinder

# ORFs are typically obtained from an ORFCollection
seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = findorfs(seq)

# Access individual ORF
orf = collection[1]

# Extract sequence through the collection
orfseq = sequence(collection, 1)
```


See also: [`ORFCollection`](/orftype#ORFCollection), [`sequence`](/api#GeneFinder.sequence-Tuple{ORFCollection,%20Int64}), [`features`](/api#GeneFinder.features-Union{Tuple{OpenReadingFrame{F}},%20Tuple{F}}%20where%20F)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L68-L102" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.Strand' href='#GeneFinder.Strand'><span class="jlbinding">GeneFinder.Strand</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Strand
```


An enumeration type representing DNA strand orientation.

**Values**
- `PSTRAND = 1`: Positive/forward strand (+)
  
- `NSTRAND = 2`: Negative/reverse strand (-)
  

**Example**

```julia
strand = PSTRAND  # Positive strand
strand = NSTRAND  # Negative strand
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L18-L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._extract_sequence-Union{Tuple{F}, Tuple{Union{BioSequences.LongSequence{var"#s26"}, BioSequences.LongSubSeq{var"#s26"}} where var"#s26"<:BioSequences.NucleicAcidAlphabet, OpenReadingFrame{F}}} where F' href='#GeneFinder._extract_sequence-Union{Tuple{F}, Tuple{Union{BioSequences.LongSequence{var"#s26"}, BioSequences.LongSubSeq{var"#s26"}} where var"#s26"<:BioSequences.NucleicAcidAlphabet, OpenReadingFrame{F}}} where F'><span class="jlbinding">GeneFinder._extract_sequence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_extract_sequence(seq::NucleicSeqOrView, orf::ORF) -> DNA sequence
```


Internal function to extract the DNA sequence for an ORF from a source sequence.

Handles strand orientation:
- `PSTRAND`: Returns a view (no allocation)
  
- `NSTRAND`: Returns reverse complement (allocates)
  

Includes bounds checking and codon validation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L342-L352" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._isvalidorf-Tuple{UnitRange{Int64}, Strand, Int8, NamedTuple}' href='#GeneFinder._isvalidorf-Tuple{UnitRange{Int64}, Strand, Int8, NamedTuple}'><span class="jlbinding">GeneFinder._isvalidorf</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_isvalidorf(range, strand, frame, features) -> Bool
```


Validate ORF parameters and throw descriptive errors for invalid inputs.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L38-L42" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.features-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.features-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.features</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
features(orf::ORF{F}) where {F}
```


Extract the features/metadata from an ORF.

**Arguments**
- `orf::ORF{F}`: The ORF from which to extract features.
  

**Returns**
- `NamedTuple`: The features associated with the ORF (may be empty).
  

**Example**

```julia
orf = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1), (score=0.95, gc=0.52))

feats = features(orf)  # Returns (score = 0.95, gc = 0.52)
feats.score            # Access individual feature: 0.95
```


See also: [`ORF`](@ref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L384-L404" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.finder-Union{Tuple{ORFCollection{F}}, Tuple{F}} where F' href='#GeneFinder.finder-Union{Tuple{ORFCollection{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.finder</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
finder(collection::ORFCollection{F}) where {F}
```


Get the gene finding method type used for this collection.

**Returns**
- `Type{F}`: The gene finder method type (e.g., `NaiveFinder`).
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L227-L234" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.finder-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.finder-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.finder</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
finder(orf::ORF{F}) where {F}
```


Get the gene finding method type used to identify this ORF.

**Arguments**
- `orf::ORF{F}`: The ORF to query.
  

**Returns**
- `Type{F}`: The gene finder method type (e.g., `NaiveFinder`).
  

**Example**

```julia
orf = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
finder(orf)  # Returns NaiveFinder
```


See also: [`GeneFinderMethod`](/api#GeneFinder.GeneFinderMethod), [`ORF`](@ref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L501-L519" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.frame-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.frame-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.frame</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
frame(orf::ORF{F}) where {F}
```


Get the reading frame of the ORF.

**Arguments**
- `orf::ORF{F}`: The ORF to query.
  

**Returns**
- `Int8`: The reading frame (1, 2, or 3).
  

**Example**

```julia
orf = ORF{NaiveFinder}(1:33, PSTRAND, Int8(2))
frame(orf)  # Returns 2
```


See also: [`strand`](/api#GeneFinder.strand-Union{Tuple{OpenReadingFrame{F}},%20Tuple{F}}%20where%20F), [`ORF`](@ref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L455-L473" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.leftposition-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.leftposition-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.leftposition</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
leftposition(orf::ORF{F}) where {F}
```


Get the left (start) position of the ORF range.

**Arguments**
- `orf::ORF{F}`: The ORF to query.
  

**Returns**
- `Int`: The first position of the ORF range.
  

**Example**

```julia
orf = ORF{NaiveFinder}(10:42, PSTRAND, Int8(1))
leftposition(orf)  # Returns 10
```


See also: [`rightposition`](/api#GeneFinder.rightposition-Union{Tuple{OpenReadingFrame{F}},%20Tuple{F}}%20where%20F), [`ORF`](@ref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L409-L427" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.orfvector-Tuple{ORFCollection}' href='#GeneFinder.orfvector-Tuple{ORFCollection}'><span class="jlbinding">GeneFinder.orfvector</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
orfvector(collection::ORFCollection)
```


Get the vector of ORFs from a collection.

**Arguments**
- `collection::ORFCollection`: The collection to query.
  

**Returns**
- `Vector{ORF{F}}`: The vector of ORFs.
  

**Example**

```julia
collection = findorfs(seq)
orf_vector = orfvector(collection)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L208-L224" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.rightposition-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.rightposition-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.rightposition</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rightposition(orf::ORF{F}) where {F}
```


Get the right (end) position of the ORF range.

**Arguments**
- `orf::ORF{F}`: The ORF to query.
  

**Returns**
- `Int`: The last position of the ORF range.
  

**Example**

```julia
orf = ORF{NaiveFinder}(10:42, PSTRAND, Int8(1))
rightposition(orf)  # Returns 42
```


See also: [`leftposition`](/api#GeneFinder.leftposition-Union{Tuple{OpenReadingFrame{F}},%20Tuple{F}}%20where%20F), [`ORF`](@ref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L432-L450" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.sequence-Tuple{ORFCollection, Int64}' href='#GeneFinder.sequence-Tuple{ORFCollection, Int64}'><span class="jlbinding">GeneFinder.sequence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sequence(collection::ORFCollection, i::Int)
```


Extract the DNA sequence for the ORF at index `i` in the collection.

**Arguments**
- `collection::ORFCollection`: The collection containing the ORF and source sequence.
  
- `i::Int`: The index of the ORF.
  

**Returns**
- `LongSubSeq{DNAAlphabet{4}}`: The DNA sequence corresponding to the ORF as a view.
  

**Example**

```julia
collection = findorfs(seq)
orfseq = sequence(collection, 1)  # Get sequence of first ORF
```


**Example**

For getting the sequence of all ORFs are several alternatives:

```julia
collection = findorfs(seq)
# Using a for loop with push!
orfseqs = Vector{LongSubSeq{DNAAlphabet{4}}}()
for orf in collection
    push!(orfseqs, sequence(collection, orf))
end

# Using broadcasting
orfseq = sequence.(Ref(collection), collection.orfs)

# Using list comprehension
orfseqs = [sequence(collection, orf) for orf in collection]

# Using map
orfseqs = map(orf -> sequence(collection, orf), collection)

# Using indices
orfseqs = [sequence(collection, i) for i in eachindex(collection)]

# Using map with indices
orfseqs = map(i -> sequence(collection, i), eachindex(collection))

# Using a generator expression
orfseqs = collect(sequence(collection, orf) for orf in collection)

# Using a generator with indices
orfseqs = collect(sequence(collection, i) for i in eachindex(collection))
```


See also: [`sequences`](@ref), [`ORFCollection`](/orftype#ORFCollection)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L259-L312" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.sequence-Tuple{ORFCollection, OpenReadingFrame}' href='#GeneFinder.sequence-Tuple{ORFCollection, OpenReadingFrame}'><span class="jlbinding">GeneFinder.sequence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sequence(collection::ORFCollection, orf::ORF)
```


Extract the DNA sequence for a specific ORF using the collection&#39;s source.

**Arguments**
- `collection::ORFCollection`: The collection containing the source sequence.
  
- `orf::ORF`: The ORF for which to extract the sequence.
  

**Returns**
- The DNA sequence (LongSubSeq{DNAAlphabet{4}}) corresponding to the ORF.
  

**Example**

```julia
collection = findorfs(seq)
orf = collection[1]
orfseq = sequence(collection, orf)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L319-L337" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.source-Tuple{ORFCollection}' href='#GeneFinder.source-Tuple{ORFCollection}'><span class="jlbinding">GeneFinder.source</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
source(collection::ORFCollection)
```


Get the source sequence associated with an ORF collection.

**Arguments**
- `collection::ORFCollection`: The collection to query.
  

**Returns**
- The source DNA sequence.
  

**Example**

```julia
collection = findorfs(seq)
src = source(collection)  # Returns the original sequence
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L189-L205" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.strand-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F' href='#GeneFinder.strand-Union{Tuple{OpenReadingFrame{F}}, Tuple{F}} where F'><span class="jlbinding">GeneFinder.strand</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
strand(orf::ORF{F}) where {F}
```


Get the strand orientation of the ORF.

**Arguments**
- `orf::ORF{F}`: The ORF to query.
  

**Returns**
- `Strand`: Either `PSTRAND` (positive/forward) or `NSTRAND` (negative/reverse).
  

**Example**

```julia
orf = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
strand(orf)  # Returns PSTRAND
```


See also: [`frame`](/api#GeneFinder.frame-Union{Tuple{OpenReadingFrame{F}},%20Tuple{F}}%20where%20F), [`Strand`](/api#GeneFinder.Strand), [`ORF`](@ref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/types.jl#L478-L496" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Finding ORFs {#Finding-ORFs}

The `findorfs` function serves as a unified interface for different gene finding methods. All methods return an `ORFCollection`.
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{F}, Tuple{N}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{F}, Tuple{N}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.findorfs</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
findorfs(sequence::NucleicSeqOrView{DNAAlphabet{N}}; finder::Type{F}, kwargs...) where {N, F<:GeneFinderMethod}
```


Main interface for finding Open Reading Frames (ORFs) in a DNA sequence.

Returns an `ORFCollection` containing the found ORFs bundled with the source sequence, providing a clean API for sequence extraction.

**Arguments**
- `sequence::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.
  

**Keywords**
- `finder::Type{F}=NaiveFinder`: The algorithm to use (`NaiveFinder`, `NaiveFinderLazy`, etc.).
  
- `alternative_start::Bool=false`: Whether to consider alternative start codons (GTG, TTG).
  
- `minlen::Int=6`: Minimum ORF length in nucleotides.
  

**Returns**
- `ORFCollection{F}`: Collection of ORFs bundled with the source sequence.
  

**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = findorfs(seq)

# Access ORFs sequences
orfseqs = sequence.(Ref(collection), collection.orfs)  # Get sequences of all ORFs

# Use a different finder
collection = findorfs(seq, finder=NaiveFinderLazy)
```


See also: [`NaiveFinder`](/api#NaiveFinder), [`NaiveFinderLazy`](/api#GeneFinder.NaiveFinderLazy), [`ORFCollection`](/orftype#ORFCollection)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/findorfs.jl#L3-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## ORF Finding Algorithms {#ORF-Finding-Algorithms}

### NaiveFinder {#NaiveFinder}

Uses regular expression matching to find ORFs.
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.NaiveFinder' href='#GeneFinder.NaiveFinder'><span class="jlbinding">GeneFinder.NaiveFinder</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
NaiveFinder <: GeneFinderMethod
```


A simple ORF finding method that detects all Open Reading Frames in a DNA sequence using regular expression matching.

See also: [`NaiveFinderLazy`](/api#GeneFinder.NaiveFinderLazy), [`GeneFinderMethod`](/api#GeneFinder.GeneFinderMethod)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L3-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder.NaiveFinder-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder.NaiveFinder</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
NaiveFinder(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> ORFCollection{NaiveFinder}
```


Find all Open Reading Frames (ORFs) in a DNA sequence using regular expression matching.

Returns an `ORFCollection` containing the source sequence and all found ORFs, providing a clean API for sequence extraction.

**Arguments**
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search for ORFs.
  

**Keywords**
- `alternative_start::Bool=false`: If `true`, uses extended start codons (ATG, GTG, TTG).
  
- `minlen::Int64=6`: Minimum ORF length in nucleotides.
  

**Returns**
- `ORFCollection{NaiveFinder}`: Collection of ORFs bundled with the source sequence.
  

**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = NaiveFinder(seq)

# Iterate over ORFs
for orf in collection
    println(sequence(collection, orf))
end

# Index access
first_seq = sequence(collection, 1)
```


See also: [`NaiveFinderLazy`](/api#GeneFinder.NaiveFinderLazy), [`ORFCollection`](/orftype#ORFCollection), [`sequence`](/api#GeneFinder.sequence-Tuple{ORFCollection,%20Int64})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L62-L97" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.NaiveFinderLazy' href='#GeneFinder.NaiveFinderLazy'><span class="jlbinding">GeneFinder.NaiveFinderLazy</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
NaiveFinderLazy <: GeneFinderMethod
```


A memory-optimized variant of `NaiveFinder` that pre-allocates the ORF vector based on estimated start codon counts.

See also: [`NaiveFinder`](/api#NaiveFinder), [`GeneFinderMethod`](/api#GeneFinder.GeneFinderMethod)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L13-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.NaiveFinderLazy-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder.NaiveFinderLazy-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder.NaiveFinderLazy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
NaiveFinderLazy(seq::NucleicSeqOrView{DNAAlphabet{N}}; kwargs...) where {N} -> ORFCollection{NaiveFinderLazy}
```


Memory-optimized ORF finder with smart pre-allocation.

**Returns**
- `ORFCollection{NaiveFinderLazy}`: Collection of ORFs bundled with the source sequence.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L208-L215" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._estimate_orf_count-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder._estimate_orf_count-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder._estimate_orf_count</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_estimate_orf_count(seq::NucleicSeqOrView{DNAAlphabet{N}}) where {N} -> Int
```


Estimate the number of ORFs in a sequence for vector pre-allocation.

Counts ATG start codons (and their reverse complement CAT) using a k-mer iterator to provide a heuristic estimate for the expected number of ORFs.

**Arguments**
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to analyze.
  

**Returns**
- `Int`: Estimated ORF count (minimum of 10 to avoid zero allocation).
  

**Implementation**

Uses `FwRvIterator` with 3-mers to efficiently count start codon occurrences on both strands in a single pass.

See also: [`NaiveFinderLazy`](/api#GeneFinder.NaiveFinderLazy)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L126-L145" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._locationiterator-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N' href='#GeneFinder._locationiterator-Union{Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}}, Tuple{N}} where N'><span class="jlbinding">GeneFinder._locationiterator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_locationiterator(seq::NucleicSeqOrView{DNAAlphabet{N}}; alternative_start::Bool=false) where {N}
```


Create an iterator that yields ORF location ranges in a DNA sequence.

**Arguments**
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search.
  

**Keywords**
- `alternative_start::Bool=false`: If `true`, matches NTG (ATG, GTG, TTG) as start codons; if `false`, only matches ATG.
  

**Returns**

An iterator yielding `UnitRange{Int64}` objects representing ORF locations.

**Implementation Details**

The function uses a regular expression to find ORFs:
- Pattern: `ATG(?:[N]{3})*?T(AG|AA|GA)` (or `NTG...` with alternative starts)
  
- `ATG` or `NTG`: Start codon
  
- `(?:[N]{3})*?`: Non-greedy match of any number of 3-nucleotide codons
  
- `T(AG|AA|GA)`: Stop codons (TAA, TAG, TGA)
  

The iterator uses `findfirst` with progressive offsets to find non-overlapping ORFs.

See also: [`NaiveFinder`](/api#NaiveFinder)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L23-L48" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder._search_strand!-Union{Tuple{N}, Tuple{Vector{OpenReadingFrame{NaiveFinderLazy}}, Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Strand, Int64, Bool, Int64}} where N' href='#GeneFinder._search_strand!-Union{Tuple{N}, Tuple{Vector{OpenReadingFrame{NaiveFinderLazy}}, Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Strand, Int64, Bool, Int64}} where N'><span class="jlbinding">GeneFinder._search_strand!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_search_strand!(orfs, seq, strand, seqlen, alternative_start, minlen)
```


Search for ORFs on a single strand and append results to the ORF vector.

This is an internal helper function that avoids code duplication between forward and reverse strand searching in `NaiveFinderLazy`.

**Arguments**
- `orfs::Vector{ORF{NaiveFinderLazy}}`: Vector to append found ORFs to (mutated).
  
- `seq::NucleicSeqOrView{DNAAlphabet{N}}`: The DNA sequence to search.
  
- `strand::Strand`: The strand being searched (`PSTRAND` or `NSTRAND`).
  
- `seqlen::Int`: Length of the original sequence (for coordinate transformation).
  
- `alternative_start::Bool`: Whether to use alternative start codons.
  
- `minlen::Int64`: Minimum ORF length filter.
  

**Coordinate Handling**
- For `PSTRAND`: Coordinates are used directly.
  
- For `NSTRAND`: Coordinates are transformed from reverse complement positions back to original sequence positions.
  

See also: [`NaiveFinderLazy`](/api#GeneFinder.NaiveFinderLazy)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/algorithms/naivefinder.jl#L162-L184" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Writing ORFs to Files {#Writing-ORFs-to-Files}

Export ORFs in various formats (FASTA, BED, GFF3).
<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_bed-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_bed-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_bed</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_bed(input, output; kwargs...)
```


Write ORF coordinates in BED format.

BED (Browser Extensible Data) format is a tab-delimited text format commonly used for genomic annotations. Each line contains: start, stop, strand, and frame.

**Arguments**
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
  
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).
  

**Keywords**
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
  
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
  
- `minlen::Int64=6`: Minimum ORF length in nucleotides.
  

**Output Format**

```
start	stop	strand	frame
35	79	+	2
120	180	-	3
```


**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file
write_orfs_bed(seq, "output.bed"; finder=NaiveFinder)

# Write to IO stream
open("output.bed", "w") do io
    write_orfs_bed(seq, io; finder=NaiveFinder)
end
```


See also: [`write_orfs_gff`](/api#GeneFinder.write_orfs_gff-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`write_orfs_fna`](/api#GeneFinder.write_orfs_fna-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`findorfs`](/api#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{F},%20Tuple{N}}%20where%20{N,%20F<:GeneFinderMethod})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/io.jl#L3-L43" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_faa-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_faa-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_faa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_faa(input, output; kwargs...)
```


Write translated ORF sequences in FASTA format (amino acids).

Outputs the protein sequences of all detected ORFs, translating each ORF using the specified genetic code. Headers contain the same metadata as `write_orfs_fna`.

**Arguments**
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
  
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).
  

**Keywords**
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
  
- `code::GeneticCode=ncbi_trans_table[1]`: Genetic code for translation (default: standard code).
  
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
  
- `minlen::Int64=6`: Minimum ORF length in nucleotides.
  

**Output Format**

```
>ORF01 id=01 start=35 stop=79 strand=+ frame=2 features=[]
MHACA*
>ORF02 id=02 start=120 stop=180 strand=- frame=3 features=[]
MLALA*
```


**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file with standard genetic code
write_orfs_faa(seq, "output.faa"; finder=NaiveFinder)

# Use bacterial genetic code (table 11)
write_orfs_faa(seq, "output.faa"; finder=NaiveFinder, code=ncbi_trans_table[11])

# Write to IO stream
open("output.faa", "w") do io
    write_orfs_faa(seq, io; finder=NaiveFinder)
end
```


See also: [`write_orfs_fna`](/api#GeneFinder.write_orfs_fna-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`write_orfs_gff`](/api#GeneFinder.write_orfs_gff-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`findorfs`](/api#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{F},%20Tuple{N}}%20where%20{N,%20F<:GeneFinderMethod}), [`BioSequences.translate`](@extref)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/io.jl#L153-L198" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_fna-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_fna-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_fna</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_fna(input, output; kwargs...)
```


Write ORF nucleotide sequences in FASTA format.

Outputs the DNA sequences of all detected ORFs, with headers containing metadata (ID, coordinates, strand, frame, and features).

**Arguments**
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
  
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).
  

**Keywords**
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
  
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
  
- `minlen::Int64=6`: Minimum ORF length in nucleotides.
  

**Output Format**

```
>ORF01 id=01 start=35 stop=79 strand=+ frame=2 features=[]
ATGCATGCATGCATGCTAG
>ORF02 id=02 start=120 stop=180 strand=- frame=3 features=[]
ATGCTAGCTAGCTAGCTAA
```


**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file
write_orfs_fna(seq, "output.fna"; finder=NaiveFinder)

# Write to IO stream
open("output.fna", "w") do io
    write_orfs_fna(seq, io; finder=NaiveFinder, minlen=9)
end
```


See also: [`write_orfs_faa`](/api#GeneFinder.write_orfs_faa-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`write_orfs_gff`](/api#GeneFinder.write_orfs_gff-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`findorfs`](/api#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{F},%20Tuple{N}}%20where%20{N,%20F<:GeneFinderMethod})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/io.jl#L72-L113" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GeneFinder.write_orfs_gff-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}' href='#GeneFinder.write_orfs_gff-Union{Tuple{F}, Tuple{N}, Tuple{Union{BioSequences.LongDNA{N}, BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}, Union{IOStream, IOBuffer}}} where {N, F<:GeneFinderMethod}'><span class="jlbinding">GeneFinder.write_orfs_gff</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
write_orfs_gff(input, output; kwargs...)
```


Write ORF annotations in GFF3 (General Feature Format version 3).

GFF3 is a standard format for genomic annotations, compatible with genome browsers like IGV, JBrowse, and the UCSC Genome Browser.

**Arguments**
- `input::NucleicSeqOrView{DNAAlphabet{N}}`: Input DNA sequence to search for ORFs.
  
- `output::Union{IOStream, IOBuffer, String}`: Output destination (file path or IO stream).
  

**Keywords**
- `finder::Type{F}=NaiveFinder`: ORF finding algorithm to use.
  
- `seqname::String="Chr"`: Sequence/chromosome name for the first GFF column.
  
- `alternative_start::Bool=false`: Use alternative start codons (GTG, TTG) in addition to ATG.
  
- `minlen::Int64=6`: Minimum ORF length in nucleotides.
  

**Output Format**

```
##gff-version 3
##sequence-region Chr 1 1000
Chr	.	ORF	35	79	.	+	.	ID=ORF01;Name=ORF01;Frame=2;Features=[]
Chr	.	ORF	120	180	.	-	.	ID=ORF02;Name=ORF02;Frame=3;Features=[]
```


**Example**

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

# Write to file
write_orfs_gff(seq, "output.gff"; finder=NaiveFinder)

# Custom chromosome name for genome browser visualization
write_orfs_gff(seq, "output.gff"; finder=NaiveFinder, seqname="scaffold_1")

# Write to IO stream
open("output.gff", "w") do io
    write_orfs_gff(seq, io; finder=NaiveFinder)
end
```


See also: [`write_orfs_bed`](/api#GeneFinder.write_orfs_bed-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`write_orfs_fna`](/api#GeneFinder.write_orfs_fna-Union{Tuple{F},%20Tuple{N},%20Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}},%20Union{IOStream,%20IOBuffer}}}%20where%20{N,%20F<:GeneFinderMethod}), [`findorfs`](/api#GeneFinder.findorfs-Union{Tuple{Union{BioSequences.LongDNA{N},%20BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}}},%20Tuple{F},%20Tuple{N}}%20where%20{N,%20F<:GeneFinderMethod})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/camilogarciabotero/GeneFinder.jl/blob/9fe625d497391ff359a6096ea0669707e340ffc7/src/io.jl#L240-L285" target="_blank" rel="noreferrer">source</a></Badge>

</details>

