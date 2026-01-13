## The ORF Type

The `ORF` type in `GeneFinder` is designed to enforce stricter validation, preventing the creation of incompatible instances. This ensures greater stability and consistency when working with open reading frames (ORFs). For example, attempting to create an `ORF` with invalid parameters will result in an error:

```julia
ORF{NaiveFinder}(:seq, 1:10, PSTRAND, Int8(4), (;)) # Frame 4 is expected to be invalid...

┌ Warning: The source sequence 'seq' is not defined. Make sure to define it or supply the correct source sequence identifier.
└ @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:46
ERROR: ArgumentError: Invalid frame: 4, expected 1, 2, or 3
Stacktrace:
 [1] _isvalidorf(seqid::Symbol, range::UnitRange{Int64}, strand::Strand, frame::Int8, features::@NamedTuple{})
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:49
 [2] OpenReadingFrame{NaiveFinder}(seqid::Symbol, range::UnitRange{Int64}, strand::Strand, frame::Int8, features::@NamedTuple{})
   @ GeneFinder ~/.julia/dev/GeneFinder/src/types.jl:100
 [3] top-level scope
   @ REPL[67]:1
```

Other sanity checks include verifying that the strand is either `PSTRAND` or `NSTRAND`, and ensuring that the features parameter is a `NamedTuple`. These checks help maintain the integrity of the `ORF` instances created within the package.

### Validation Rules

The `ORF` type enforces the following rules:
- Parameters must align with the expected structure of an open reading frame.
- The strand must be either `PSTRAND` or `NSTRAND`. Any other value will trigger an error.
- The frame must be an `Int8` with a value of 1, 2, or 3.
- The features parameter must be a `NamedTuple`. If not, an error will be raised.
- The length of the ORF must be at least 6 nucleotides to accommodate a start and stop codon.
- The sequence identifier must correspond to a defined sequence in the context where the `ORF` is used. If not defined, a warning will be issued.
- The range must be divisible by 3 to ensure it represents a valid reading frame.


These safeguards help prevent the creation of invalid ORFs, making the system more robust and easier to extend.

### Example Usage

Here is an example of creating a valid `ORF` instance:

```julia
seq = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"

ORF{NaiveFinder}(:seq, 1:33, PSTRAND, Int8(1), (;))
```

As long as the sequence identifier `:seq` corresponds to a defined sequence in your context, this instance can then be used to iterate over a sequence of interest and verify whether the `ORF` is contained within the sequence. If not defined, a warning will be issued to ensure clarity.

This instance can then be used to iterate over a sequence of interest and verify whether the `ORF` is contained within the sequence.

!!! warning
    While the `ORF` type ensures stricter validation, it is still possible to create an `ORF` that does not correspond to an actual open reading frame in the sequence. This limitation will be addressed in future versions of the package. However, the current implementation provides the benefit of retrieving the corresponding subsequence in a convenient 5' to 3' orientation, regardless of the strand.

