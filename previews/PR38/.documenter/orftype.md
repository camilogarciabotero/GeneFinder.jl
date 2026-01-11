
## The ORFI Type {#The-ORFI-Type}

The `ORFI` type in `GeneFinder` is designed to enforce stricter validation, preventing the creation of incompatible instances. This ensures greater stability and consistency when working with open reading frames (ORFs). For example, attempting to create an `ORFI` with invalid parameters will result in an error:

```julia
ORFI{4,NaiveFinder}(1:10, '+', 4) # Or any F <: GeneFinderMethod

ERROR: MethodError: no method matching OpenReadingFrameInterval{4, NaiveFinder}(::UnitRange{Int64}, ::Char, ::Int64)
The type `OpenReadingFrameInterval{4, NaiveFinder}` exists, but no method is defined for this combination of argument types when trying to construct it.

Closest candidates are:
  (::Type{OpenReadingFrameInterval{N, F}} where {N, F<:GeneFinder.GeneFinderMethod})(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:51
  OpenReadingFrameInterval{N, F}(::Type{F}, ::String, ::Int64, ::Int64, ::GenomicFeatures.Strand, ::Int8, ::BioSequences.LongSubSeq{BioSequences.DNAAlphabet{N}}, ::NamedTuple) where {N, F<:GeneFinder.GeneFinderMethod}
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:64

Stacktrace:
 [1] top-level scope
   @ REPL[3]:1
```


### Validation Rules {#Validation-Rules}

The `ORFI` type enforces the following rules:
- The strand must be either `+` or `-`. Any other value will trigger an error.
  
- Parameters must align with the expected structure of an open reading frame.
  

These safeguards help prevent the creation of invalid ORFs, making the system more robust and easier to extend.

### Example Usage {#Example-Usage}

Here is an example of creating a valid `ORFI` instance:

```julia
ORFI{4,NaiveFinder}(
    "seq", 
    1, 
    33, 
    STRAND_POS, 
    1, 
    convert(LongSubSeq, dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), 
    NamedTuple()
)
```


This instance can then be used to iterate over a sequence of interest and verify whether the `ORFI` is contained within the sequence.

::: warning Warning

While the `ORFI` type ensures stricter validation, it is still possible to create an `ORFI` that does not correspond to an actual open reading frame in the sequence. This limitation will be addressed in future versions of the package. However, the current implementation provides the benefit of retrieving the corresponding subsequence in a convenient 5&#39; to 3&#39; orientation, regardless of the strand.

:::
