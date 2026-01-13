
## The ORF and ORFCollection Types {#The-ORF-and-ORFCollection-Types}

GeneFinder uses two complementary types for representing Open Reading Frames:
- **`ORF`**: A lightweight, coordinate-based structure storing location and metadata
  
- **`ORFCollection`**: A container bundling ORFs with a view of their source sequence
  

This design ensures ORFs are composable, serializable, and independent of global state.

### ORFCollection {#ORFCollection}

All `GeneFinderMethod` implementations return an `ORFCollection`:

```julia
using BioSequences, GeneFinder

seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
collection = findorfs(seq)

# The collection bundles ORFs with a view of the source
typeof(collection)  # ORFCollection{NaiveFinder, LongSubSeq{DNAAlphabet{4}}}

# Source is always stored as a view (memory efficient)
typeof(source(collection))  # LongSubSeq{DNAAlphabet{4}}

# Access ORFs
orfs(collection)      # Vector of ORFs
collection[1]         # First ORF
length(collection)    # Number of ORFs
```


::: tip Note

The source sequence is always stored as a `LongSubSeq` (view), which avoids copying the sequence data while maintaining a reference to the original. This is both memory-efficient and ensures the collection stays synchronized with the source.

:::

### Sequence Extraction {#Sequence-Extraction}

The `ORFCollection` provides clean sequence extraction using the `sequence` function:

```julia
# Extract by index
seq1 = sequence(collection, 1)

# Extract by ORF
orf = collection[1]
seq1 = sequence(collection, orf)
```


### Translation {#Translation}

Translate ORF sequences to amino acids using `BioSequences.translate`:

```julia
using BioSequences

# Translate a single ORF
protein = translate(sequence(collection, 1))

# Translate multiple ORFs
proteins = [translate(sequence(collection, i)) for i in eachindex(collection)]
```


### The ORF Type {#The-ORF-Type}

Individual ORFs store only coordinates and metadata:

```julia
orf = collection[1]

leftposition(orf)   # Start position
rightposition(orf)  # End position
strand(orf)         # PSTRAND or NSTRAND
frame(orf)          # Reading frame (1, 2, or 3)
features(orf)       # NamedTuple of features
finder(orf)         # The method that found it
```


### Validation Rules {#Validation-Rules}

The `ORF` constructor validates:
- Strand must be `PSTRAND` or `NSTRAND`
  
- Frame must be 1, 2, or 3
  
- Range length must be divisible by 3
  
- Range length must be ≥ 6 (start + stop codons)
  

```julia
# This will error - frame must be 1, 2, or 3
ORF{NaiveFinder}(1:33, PSTRAND, Int8(4))
# ERROR: ArgumentError: Invalid frame: 4, expected 1, 2, or 3
```


### Design Benefits {#Design-Benefits}

This architecture provides:
- **No global state**: Source sequence is explicitly bundled
  
- **Composability**: Works in any function or module context
  
- **Serialization**: Collections can be saved/loaded
  
- **Type safety**: Compile-time guarantees on finder method
  
- **Efficiency**: Views used where possible, validation at construction
  
