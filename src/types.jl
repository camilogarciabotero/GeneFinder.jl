# Structs associated with gene models 
# using GenomicFeatures
using BioSequences

abstract type Gene end

# abstract type exon end

# abstract type intron end

"""
    struct ORF
        location::UnitRange{Int64}
        strand::Char

        ORF(location, strand) = new(location, strand)
    end

The ORF struct represents an open reading frame in a DNA sequence. It has two fields: 

- `location`: which is a UnitRange{Int64} indicating the start and end locations of the ORF in the sequence
- `strand`:  is a Char type indicating whether the ORF is on the forward ('+') or reverse ('-') strand of the sequence.
"""
struct ORF <: Gene
    location::UnitRange{Int64} # Note that it is also called position for gene struct in GenomicAnotations
    strand::Char

    ORF(location, strand) = new(location, strand)
end

# # The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/
"""
    Codon <: BioSequence{DNAAlphabet{2}}

A `Struct` representing a codon, which is a subtype of `BioSequence` with
an `Alphabet` of type `DNAAlphabet{2}`. It has a single field `x` of type
`UInt8`. This was implemente in The following implementation is from https://biojulia.net/BioSequences.jl/stable/interfaces/
"""
struct Codon <: BioSequence{DNAAlphabet{2}}
    x::UInt8
end

"""
    Codon(iterable)

Constructor function for `Codon` instances. Takes an `iterable` as input and
creates a new `Codon` instance. It checks that the length of the `iterable` is
3, and if not, it throws an error. It then creates a new `UInt` variable `x`
and assigns it a value of 0. It then iterates over the `iterable` and encodes
each element using the `BioSequences.encode` function and the `Alphabet` of the
`Codon` type. It shifts the encoded value by 6 minus 2 times the index of the
element, and then ORs it with the current value of `x`. Finally, it creates a
new `Codon` instance with the value of `x` modulo `UInt8`.
"""
function Codon(iterable)
    length(iterable) == 3 || error("Must have length 3")
    x = zero(UInt)
    for (i, nt) in enumerate(iterable)
        x |= BioSequences.encode(Alphabet(Codon), convert(DNA, nt)) << (6-2i)
    end
    Codon(x % UInt8)
end

"""
    length(::Codon)

Returns the length of a `Codon` instance. It is defined as a constant with a
value of 3.
"""
Base.length(::Codon) = 3


"""
    extract_encoded_element(x::Codon, i::Int)

Part of the `BioSequences` module. Takes a `Codon` instance `x` and an integer
`i` as input. Returns the value of `x.x` shifted right by 6 minus 2 times `i`,
ANDed with 3, and then modulo `UInt`.
"""
function BioSequences.extract_encoded_element(x::Codon, i::Int)
    ((x.x >>> (6-2i)) & 3) % UInt
end

"""
    copy(seq::Codon)

A `Base` function that takes a `Codon` instance `seq` as input and returns a
new `Codon` instance with the same value as `seq.x`.
"""
Base.copy(seq::Codon) = Codon(seq.x)

BioSequences.has_interface(BioSequence, Codon, [DNA_C, DNA_T, DNA_G], false)

"""
stopcodons

A constant array of `Codon` instances representing the three stop codons in the
standard genetic code: TAG, TAA, and TGA.
"""
const stopcodons = [Codon("TAG"), Codon("TAA"), Codon("TGA")]

# const START_CODON_MATRIX = PWMSearchQuery([dna"ATG", dna"GTG", dna"TTG"], 1.0)

"""
    struct CDS
        orf::ORF
        sequence::LongDNA
    end

The `CDS` struct represents a coding sequence in a DNA sequence. It has three fields:

- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) displaying the location of the ORF and the associate strand: forward ('+') or reverse ('-')
- `sequence`: a `LongDNA` sequence representing the actual sequence of the CDS
"""
struct CDS <: Gene
    orf::ORF
    sequence::LongDNA
end

"""
    struct Protein
        orf::ORF
        sequence::LongDNA
    end


Similarly to the `CDS` struct, the `Protein` struct represents a encoded protein sequence in a DNA sequence. 
    It has three fields:

- `orf`: is the basic composible type (`location::UnitRange{Int}`, strand::Char) of the sequence
- `sequence`: a `LongAA` sequence representing the actual translated sequence of the CDS
"""
struct Protein <: Gene
    orf::ORF
    sequence::LongAA
end