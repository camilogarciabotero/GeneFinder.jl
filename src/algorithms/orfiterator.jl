seq = randdnaseq(1*10^7)

# regcds = biore"ATG(?:[N]{3})*?T(AG|AA|GA)"dna
# collect(matched(x) for x in eachmatch(regcds, seq))


function locationiterator(sequence::LongDNA; alternative_start::Bool=false)
    regcds = alternative_start ? EXTENDED_REGULAR_CDS : REGULAR_CDS
    finder(x) = findfirst(regcds, sequence, first(x)+3)
    itr = Iterators.takewhile(!isnothing, Iterators.drop(iterated(finder, -2:-2), 1))
    return itr
end

function orfiterator(sequence::LongDNA; alternative_start::Bool=false, min_len = 6)
    revseq = reverse_complement(sequence)
    @inbounds begin
        orfs = (ORF(location, strand) for strand in ['+', '-'] for location in locationiterator(strand == '+' ? sequence : revseq; alternative_start) if length(location) >= min_len) 
    end
    return orfs
end

function cdsgenerator(sequence::LongDNA; alternative_start::Bool=false, min_len=6)
    orfit = orfiterator(sequence; alternative_start, min_len)
    revseq = reverse_complement(sequence)
    @inbounds begin
        cds = (i.strand == '+' ? CDS(sequence[i.location], i) : CDS(revseq[i.location], i) for i in orfit)
    end
    return cds
end

function get_cds(sequence::LongDNA; alternative_start::Bool=false, min_len::Int64=6)
    cds = [i.sequence for i in cdsgenerator(sequence; alternative_start, min_len)]
    return cds
end


function protein_generator(sequence::LongDNA; alternative_start::Bool=false, code::GeneticCode=BioSequences.standard_genetic_code, min_len::Int64=6)
    orfit = orfiterator(sequence; alternative_start, min_len)
    reversedseq = reverse_complement(sequence)
    proteins = (i.strand == '+' ? Protein(translate(sequence[i.location]; alternative_start, code), i) : Protein(translate(reversedseq[i.location]; alternative_start, code), i) for i in orfit)
    return proteins
end