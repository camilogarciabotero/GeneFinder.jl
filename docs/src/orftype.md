
## The ORFI type

For convenience, the `ORFI` type is more stringent in preventing the creation of incompatible instances. As a result, attempting to create an instance with incompatible parameters will result in an error. For instance, the following code snippet will trigger an error:

```julia
ORFI{4,NaiveFinder}(1:10, '+', 4) # Or any F <: GeneFinderMethod

ERROR: MethodError: no method matching OpenReadingFrameInterval{4, NaiveFinder}(::UnitRange{Int64}, ::Char, ::Int64)

Closest candidates are:
  (::Type{OpenReadingFrameInterval{N, F}} where {N, F})(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:49
  OpenReadingFrameInterval{N, F}(::Type{F}, ::String, ::Int64, ::Int64, ::Strand, ::Int8, ::LongSubSeq{DNAAlphabet{N}}, ::NamedTuple) where {N, F<:GeneFinder.GeneFinderMethod}
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:58

Stacktrace:
 [1] top-level scope
   @ REPL[21]:1
```
 
Similar behavior will be encountered when the strand is neither `+` nor `-`. This precautionary measure helps prevent the creation of invalid ORFs, ensuring greater stability and enabling the extension of its interface. For example, after creating a specific `ORFI`, users can seamlessly iterate over a sequence of interest and verify whether the ORFI is contained within the sequence.

```julia
ORFI{4,NaiveFinder}("seq", 1, 33, STRAND_POS, 1, convert(LongSubSeq, dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAG"), NamedTuple())
```

!!! warning
    It is still possible to create an `ORFI` and pass it to a sequence that does not necessarily contain an actual open reading frame. This will be addressed in future versions of the package. But the benefit of having it is that it will retrieve the corresponding subsequence of the sequence in a convinient way (5' to 3') regardless of the strand.

