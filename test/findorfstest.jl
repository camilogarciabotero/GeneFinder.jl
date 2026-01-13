using Test
using BioSequences, GeneFinder, FASTX

# ════════════════════════════════════════════════════════════════════════════════
# Test Data
# ════════════════════════════════════════════════════════════════════════════════

const TEST_SEQ_SHORT = dna"ATGATGCATGCATGCATGCTAGTAACTAGCTAGCTAGCTAGTAA"
const TEST_SEQ_LONG = dna"AACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAACAGCACTGGCAATCTGACTGTGGGCGGTGTTACCAACGGCACTGCTACTACTGGCAACATCGCACTGACCGGTAACAATGCGCTGAGCGGTCCGGTCAATCTGAATGCGTCGAATGGCACGGTGACCTTGAACACGACCGGCAATACCACGCTCGGTAACGTGACGGCACAAGGCAATGTGACGACCAATGTGTCCAACGGCAGTCTGACGGTTACCGGCAATACGACAGGTGCCAACACCAACCTCAGTGCCAGCGGCAACCTGACCGTGGGTAACCAGGGCAATATCAGTACCGCAGGCAATGCAACCCTGACGGCCGGCGACAACCTGACGAGCACTGGCAATCTGACTGTGGGCGGCGTCACCAACGGCACGGCCACCACCGGCAACATCGCGCTGACCGGTAACAATGCACTGGCTGGTCCTGTCAATCTGAACGCGCCGAACGGCACCGTGACCCTGAACACAACCGGCAATACCACGCTGGGTAATGTCACCGCACAAGGCAATGTGACGACTAATGTGTCCAACGGCAGCCTGACAGTCGCTGGCAATACCACAGGTGCCAACACCAACCTGAGTGCCAGCGGCAATCTGACCGTGGGCAACCAGGGCAATATCAGTACCGCGGGCAATGCAACCCTGACTGCCGGCGGTAACCTGAGC"

# ════════════════════════════════════════════════════════════════════════════════
# ORF Type Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "ORF Construction" begin
    @testset "Valid ORF creation" begin
        orf = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        
        @test leftposition(orf) == 1
        @test rightposition(orf) == 33
        @test strand(orf) === PSTRAND
        @test frame(orf) == 1
        @test finder(orf) === NaiveFinder
        @test isempty(features(orf))
    end
    
    @testset "ORF with features" begin
        orf = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1), (score=0.95, gc=0.52))
        
        @test features(orf) == (score=0.95, gc=0.52)
        @test features(orf).score == 0.95
        @test features(orf).gc == 0.52
    end
    
    @testset "ORF on negative strand" begin
        orf = ORF{NaiveFinder}(10:42, NSTRAND, Int8(1))
        
        @test strand(orf) === NSTRAND
        @test leftposition(orf) == 10
        @test rightposition(orf) == 42
    end
    
    @testset "Invalid ORF - bad frame" begin
        @test_throws ArgumentError ORF{NaiveFinder}(1:33, PSTRAND, Int8(0))
        @test_throws ArgumentError ORF{NaiveFinder}(1:33, PSTRAND, Int8(4))
        @test_throws ArgumentError ORF{NaiveFinder}(1:33, PSTRAND, Int8(-1))
    end
    
    @testset "Invalid ORF - not divisible by 3" begin
        @test_throws ArgumentError ORF{NaiveFinder}(1:32, PSTRAND, Int8(1))
        @test_throws ArgumentError ORF{NaiveFinder}(1:34, PSTRAND, Int8(1))
    end
    
    @testset "Invalid ORF - too short" begin
        @test_throws ArgumentError ORF{NaiveFinder}(1:3, PSTRAND, Int8(1))
    end
    
    @testset "Valid minimum ORF (6 nt)" begin
        orf = ORF{NaiveFinder}(1:6, PSTRAND, Int8(1))
        @test rightposition(orf) - leftposition(orf) + 1 == 6
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# ORF Comparison Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "ORF Comparisons" begin
    @testset "Equality" begin
        orf1 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        orf2 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        orf3 = ORF{NaiveFinderLazy}(1:33, PSTRAND, Int8(1))  # Different finder
        
        @test orf1 == orf2
        @test orf1 == orf3  # Equality ignores finder type
    end
    
    @testset "Inequality" begin
        orf1 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        orf2 = ORF{NaiveFinder}(1:33, NSTRAND, Int8(1))  # Different strand
        orf3 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(2))  # Different frame
        orf4 = ORF{NaiveFinder}(4:36, PSTRAND, Int8(1))  # Different range
        
        @test orf1 != orf2
        @test orf1 != orf3
        @test orf1 != orf4
    end
    
    @testset "Sorting (isless)" begin
        orf1 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        orf2 = ORF{NaiveFinder}(10:42, PSTRAND, Int8(1))
        orf3 = ORF{NaiveFinder}(1:33, NSTRAND, Int8(1))
        
        @test orf1 < orf2  # Earlier start
        @test orf1 < orf3  # Same start, PSTRAND < NSTRAND
        
        orfs = [orf2, orf3, orf1]
        sorted = sort(orfs)
        @test sorted[1] == orf1
        @test sorted[2] == orf3
        @test sorted[3] == orf2
    end
    
    @testset "Hashing (Set/Dict support)" begin
        orf1 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        orf2 = ORF{NaiveFinder}(1:33, PSTRAND, Int8(1))
        orf3 = ORF{NaiveFinder}(4:36, PSTRAND, Int8(1))
        
        s = Set([orf1, orf2, orf3])
        @test length(s) == 2  # orf1 and orf2 are duplicates
        
        d = Dict(orf1 => "first", orf3 => "second")
        @test d[orf2] == "first"  # orf2 == orf1
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# ORFCollection Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "ORFCollection" begin
    collection = findorfs(TEST_SEQ_SHORT)
    
    @testset "Basic properties" begin
        @test length(collection) == 5
        @test !isempty(collection)
        @test finder(collection) === NaiveFinder
    end
    
    @testset "Source sequence" begin
        src = source(collection)
        @test src == TEST_SEQ_SHORT
        @test typeof(src) <: LongSubSeq  # Always a view
    end
    
    @testset "Iteration" begin
        count = 0
        for orf in collection
            count += 1
            @test orf isa ORF{NaiveFinder}
        end
        @test count == length(collection)
    end
    
    @testset "Indexing" begin
        @test collection[1] isa ORF{NaiveFinder}
        @test collection[end] isa ORF{NaiveFinder}
        @test length(collection[1:3]) == 3
        @test firstindex(collection) == 1
        @test lastindex(collection) == length(collection)
    end
    
    @testset "orfs accessor" begin
        orf_vec = orfvector(collection)
        @test orf_vec isa Vector{ORF{NaiveFinder}}
        @test length(orf_vec) == length(collection)
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# Sequence Extraction Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "Sequence Extraction" begin
    collection = findorfs(TEST_SEQ_SHORT)
    
    @testset "Extract by index" begin
        seq1 = GeneFinder.sequence(collection, 1)
        @test seq1 isa LongSubSeq
        @test length(seq1) == rightposition(collection[1]) - leftposition(collection[1]) + 1
    end
    
    @testset "Extract by ORF" begin
        orf = collection[1]
        seq1 = GeneFinder.sequence(collection, orf)
        seq2 = GeneFinder.sequence(collection, 1)
        @test seq1 == seq2
    end
    
    # @testset "Extract all sequences" begin
    #     seqs = sequences(collection)
    #     @test length(seqs) == length(collection)
    #     @test all(s -> s isa Union{LongSubSeq, LongDNA}, seqs)
    # end
    
    @testset "Sequence starts with ATG" begin
        for i in 1:length(collection)
            seq = GeneFinder.sequence(collection, i)
            @test seq[1:3] == convert(LongSubSeq{DNAAlphabet{4}}, dna"ATG")
        end
    end
    
    @testset "Sequence ends with stop codon" begin
        stops = (dna"TAA", dna"TAG", dna"TGA")
        for i in 1:length(collection)
            seq = GeneFinder.sequence(collection, i)
            @test seq[end-2:end] in stops
        end
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# Translation Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "Translation" begin
    collection = findorfs(TEST_SEQ_SHORT)
    
    @testset "Translate by index" begin
        protein = translate(collection, 1)
        @test protein isa LongAA
        @test protein[1] == AA_M  # Starts with Methionine
        @test protein[end] == AA_Term  # Ends with stop
    end
    
    @testset "Translate by ORF" begin
        orf = collection[1]
        protein1 = translate(collection, orf)
        protein2 = translate(collection, 1)
        @test protein1 == protein2
    end
    
    # @testset "Translate all" begin
    #     proteins = translations(collection)
    #     @test length(proteins) == length(collection)
    #     @test all(p -> p isa LongAA, proteins)
    #     @test all(p -> p[1] == AA_M, proteins)
    #     @test all(p -> p[end] == AA_Term, proteins)
    # end
    
    @testset "Translation consistency" begin
        # Manual translation should match
        for i in 1:length(collection)
            seq = GeneFinder.sequence(collection, i)
            manual = translate(seq)
            direct = translate(collection, i)
            @test manual == direct
        end
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# ORFCollection Comparison Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "ORFCollection Comparisons" begin
    c1 = findorfs(TEST_SEQ_SHORT, finder=NaiveFinder)
    c2 = findorfs(TEST_SEQ_SHORT, finder=NaiveFinderLazy)
    
    @testset "Equality between finders" begin
        @test c1 == c2  # Same ORFs found
    end
    
    @testset "isequal (strict)" begin
        @test isequal(c1, c2)  # Same source too
    end
    
    @testset "Different minlen produces different results" begin
        c_short = findorfs(TEST_SEQ_SHORT, minlen=6)
        c_long = findorfs(TEST_SEQ_SHORT, minlen=30)
        
        @test c_short != c_long
        @test length(c_short) >= length(c_long)
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# Set Operations Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "Set Operations" begin
    c_all = findorfs(TEST_SEQ_SHORT, minlen=6)
    c_long = findorfs(TEST_SEQ_SHORT, minlen=18)
    
    @testset "issubset" begin
        @test issubset(c_long, c_all)
        @test c_long ⊆ c_all
    end
    
    @testset "intersect" begin
        common = intersect(c_all, c_long)
        @test length(common) == length(c_long)
        @test all(orf -> orf in orfs(c_long), common)
    end
    
    @testset "setdiff" begin
        short_only = setdiff(c_all, c_long)
        @test all(orf -> orf ∉ orfs(c_long), short_only)
        @test length(short_only) == length(c_all) - length(c_long)
    end
    
    @testset "union" begin
        all_orfs = union(c_all, c_long)
        @test length(all_orfs) == length(c_all)  # c_long ⊆ c_all
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# NaiveFinder Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "NaiveFinder" begin
    @testset "Basic finding" begin
        collection = findorfs(TEST_SEQ_SHORT, finder=NaiveFinder)
        
        @test length(collection) == 5
        @test finder(collection) === NaiveFinder
        
        # First ORF
        @test leftposition(collection[1]) == 1
        @test rightposition(collection[1]) == 33
        @test strand(collection[1]) === PSTRAND
        @test frame(collection[1]) == 1
    end
    
    @testset "minlen filter" begin
        c_6 = findorfs(TEST_SEQ_SHORT, finder=NaiveFinder, minlen=6)
        c_30 = findorfs(TEST_SEQ_SHORT, finder=NaiveFinder, minlen=30)
        
        @test length(c_6) >= length(c_30)
        @test all(orf -> rightposition(orf) - leftposition(orf) + 1 >= 30, orfs(c_30))
    end
    
    @testset "No ORFs found" begin
        collection = findorfs(TEST_SEQ_LONG, minlen=1000)
        @test length(collection) == 0
        @test isempty(collection)
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# NaiveFinderLazy Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "NaiveFinderLazy" begin
    @testset "Same results as NaiveFinder" begin
        c1 = findorfs(TEST_SEQ_SHORT, finder=NaiveFinder)
        c2 = findorfs(TEST_SEQ_SHORT, finder=NaiveFinderLazy)
        
        @test length(c1) == length(c2)
        @test c1 == c2
    end
    
    @testset "Finder type preserved" begin
        collection = findorfs(TEST_SEQ_SHORT, finder=NaiveFinderLazy)
        @test finder(collection) === NaiveFinderLazy
        @test all(orf -> finder(orf) === NaiveFinderLazy, collection)
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# NaiveCollector Tests
# ════════════════════════════════════════════════════════════════════════════════

@testset "NaiveCollector" begin
    @testset "Basic finding" begin
        collection = findorfs(TEST_SEQ_SHORT, finder=NaiveCollector)
        
        @test !isempty(collection)
        @test finder(collection) === NaiveCollector
    end
    
    @testset "Non-overlapping by default" begin
        collection = findorfs(TEST_SEQ_SHORT, finder=NaiveCollector, overlap=false)
        @test !isempty(collection)
    end
    
    @testset "Overlapping option" begin
        c_no_overlap = findorfs(TEST_SEQ_SHORT, finder=NaiveCollector, overlap=false)
        c_overlap = findorfs(TEST_SEQ_SHORT, finder=NaiveCollector, overlap=true)
        
        @test length(c_overlap) >= length(c_no_overlap)
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# Lambda Phage Integration Test
# ════════════════════════════════════════════════════════════════════════════════

@testset "Lambda Phage (Integration)" begin
    lambdafile = "data/NC_001416.1.fasta"
    
    if isfile(lambdafile)
        lambdarecord = open(collect, FASTA.Reader, lambdafile)[1]
        lambdaseq = FASTX.sequence(LongDNA{4}, lambdarecord)
        
        @testset "NaiveFinder on lambda" begin
            collection = findorfs(lambdaseq, finder=NaiveFinder, minlen=75)
            
            @test length(collection) == 885
            @test all(orf -> strand(orf) in (PSTRAND, NSTRAND), collection)
            @test all(orf -> frame(orf) in (1, 2, 3), collection)
        end
        
        @testset "All sequences valid" begin
            collection = findorfs(lambdaseq, finder=NaiveFinder, minlen=75)
            
            # Check a sample of sequences
            sample_size = min(50, length(collection))
            for i in 1:sample_size
                seq = GeneFinder.sequence(collection, i)
                @test length(seq) >= 75
                @test seq[1:3] == convert(LongSubSeq{DNAAlphabet{4}}, dna"ATG")
            end
        end
        
        @testset "Translation works" begin
            collection = findorfs(lambdaseq, finder=NaiveFinder, minlen=75)
            
            # Translate a sample
            sample_size = min(50, length(collection))
            for i in 1:sample_size
                protein = translate(collection, i)
                @test protein[1] == AA_M
                @test protein[end] == AA_Term
            end
        end
    else
        @warn "Lambda phage test file not found: $lambdafile"
    end
end

# ════════════════════════════════════════════════════════════════════════════════
# Edge Cases
# ════════════════════════════════════════════════════════════════════════════════

@testset "Edge Cases" begin
    @testset "Minimum valid ORF (ATGTGA)" begin
        seq = dna"ATGTGA"
        collection = findorfs(seq, minlen=6)
        
        @test length(collection) == 1
        @test GeneFinder.sequence(collection, 1) == convert(LongSubSeq{DNAAlphabet{4}}, dna"ATGTGA")
        @test translate(collection, 1) == aa"M*"
    end
    
    @testset "No start codon" begin
        seq = dna"TGCTGCTGCTGA"
        collection = findorfs(seq)
        @test isempty(collection)
    end
    
    @testset "No stop codon" begin
        seq = dna"ATGATGATGATG"
        collection = findorfs(seq)
        @test isempty(collection)
    end
    
    @testset "Empty collection operations" begin
        seq = dna"TGCTGCTGCTGA"
        collection = findorfs(seq)
        
        @test isempty(collection)
        @test length(collection) == 0
        # @test isempty(GeneFinder.sequences(collection))
        # @test isempty(GeneFinder.translations(collection))
    end
end