## The log-odds ratio decision rule

The sequence probability given a transition probability model could be used as the source of a sequence classification based on a decision rule to classify whether a sequence correspond to a model or another. Now, imagine we got two DNA sequence transition models, a CDS model and a No-CDS model. The *log-odds ratio* decision rule could be establish as:

``` math
\begin{align}
S(X) = \log \frac{{P_C(X_1=i_1, \ldots, X_T=i_T)}}{{P_N(X_1=i_1, \ldots, X_T=i_T)}}  \begin{cases} > \eta & \Rightarrow \text{coding} \\ < \eta & \Rightarrow \text{noncoding} \end{cases}
\end{align}
```

Where the ``P_{C}`` is the probability of the sequence given a CDS model, ``P_{N}`` is the probability of the sequence given a No-CDS model, the decision rule is finally based on whether the ratio is greater or lesser than a given threshold *η* of significance level.

In this package we have implemented this rule and call some basic models of CDS and No-CDS of *E. coli* from Axelson-Fisk (2015) work (implemented in `BioMarkovChains.jl` package). To check whether a random sequence could be coding based on these decision we use the predicate `log_odds_ratio_decision_rule` with the `ECOLICDS` and `ECOLINOCDS` models:

```julia
julia> using GeneFinder, BioSequences

seq = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"

orfsdna = findorfs(seq, finder=NaiveFinder, minlen=75, alternative_start=true) .|> sequence

20-element Vector{NucSeq{4, DNAAlphabet{4}}}
 ATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCC…CGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGA
 GTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTC…TACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAG
 GTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAG…CGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGA
 GTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTC…TACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAG
 TTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTT…GCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGA
 ATGCTGCCGGTAATATGCATTATGAATTTCAGGTCGATAGTGAGGATGAAATCGGTTATCTGGCTTCCTGTCTGA
 ATGCCACAGAGCAGTATGCTGCCGGTAATATGCATTATG…ATAGTGAGGATGAAATCGGTTATCTGGCTTCCTGTCTGA
 ATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGT…TACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAG
 GTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAA…ATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAG
 ATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGT…ATCAAAATTCTCAATGACAACCACATTCTCCACATCTAA
 TTGGTGATTTCTTCGGTAATTTTGATTCGGATGTACTTACTGTATTTGCTCCCATTACCTCCAACTACAAGGTAA
 TTGTTGGTGATTTCTTCGGTAATTTTGATTCGGATGTACTTACTGTATTTGCTCCCATTACCTCCAACTACAAGGTAA
 ATGACAACCACATTCTCCACATCTAAGGGACTATTGGTA…CCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAG
 ATGCCGGCGGATTGTACACCAGTGAAACAACTTTGGAGACCGTTAAAATGCAGCTGGACTCCCTGGCGGTCTACCTGA
 TTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGG…TCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAA
 TTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAA
 ATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAAT…AAATACACCAAAAATAAAATATGCGAGTATAAATTTTAG
 TTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAA…ACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAA
 TTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAA
 ATGCCAACTGGAAGATCAGCACCATCTGGGGCATCGGCT…TACGTAATTAAGTCATACCAAAAGCAGGATTTGTTATGA

```

Now the question is which of those sequences can we consider as coding sequences. We can use the `iscoding` predicate to check whether a sequence is coding or not based on the *log-odds ratio* decision rule:

```julia
julia> iscoding.(orfsdna) # criteria = log_odds_ratio_decision_rule

20-element BitVector:
 0
 0
 0
 0
 0
 1
 1
 0
 0
 0
 0
 0
 0
 1
 0
 0
 0
 0
 0
 0
```

In this case, the sequence has 20 ORFIs and only 3 of them are classified as coding sequences. The classification is based on the *log-odds ratio* decision rule and the transition probability models of *E. coli* CDS and No-CDS. The `log_odds_ratio_decision_rule` method will return a boolean vector with the classification of each ORFI in the sequence. Now we can simply filter the ORFIs that are coding sequences:

```julia
orfs = filter(orf -> iscoding(orf), orfsdna)

3-element Vector{NucSeq{4, DNAAlphabet{4}}}
 ATGCTGCCGGTAATATGCATTATGAATTTCAGGTCGATAGTGAGGATGAAATCGGTTATCTGGCTTCCTGTCTGA
 ATGCCACAGAGCAGTATGCTGCCGGTAATATGCATTATG…ATAGTGAGGATGAAATCGGTTATCTGGCTTCCTGTCTGA
 ATGCCGGCGGATTGTACACCAGTGAAACAACTTTGGAGACCGTTAAAATGCAGCTGGACTCCCTGGCGGTCTACCTGA
```

Or in terms of the `ORFI` object:

```julia
julia> orfs = findorfs(seq, minlen=75, finder=NaiveFinder, alternative_start=true) # find ORFIs with alternative start as well
julia> orfs[iscoding.(orfsdna)]

3-element Vector{ORFI{4, NaiveFinder}}:
 ORFI{NaiveFinder}(194:268, '-', 2)
 ORFI{NaiveFinder}(194:283, '-', 2)
 ORFI{NaiveFinder}(650:727, '-', 2)
```

Or in a single line using another genome sequence:

```julia
julia> phi = dna"GTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTATATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTGAATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGCCGGGCAATAACGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTTCTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATGTTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTTATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATATGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCAGAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTTTCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAGTGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAAGGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAACCTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGA"

julia> filter(x -> iscoding(sequence(x), η=1e-10) && length(x) > 100, findorfs(phi))

34-element Vector{ORFI{4, NaiveFinder}}:
 ORFI{NaiveFinder}(636:1622, '+', 3)
 ORFI{NaiveFinder}(687:1622, '+', 3)
 ORFI{NaiveFinder}(774:1622, '+', 3)
 ORFI{NaiveFinder}(781:1389, '+', 1)
 ORFI{NaiveFinder}(814:1389, '+', 1)
 ORFI{NaiveFinder}(829:1389, '+', 1)
 ORFI{NaiveFinder}(861:1622, '+', 3)
 ORFI{NaiveFinder}(1021:1389, '+', 1)
 ORFI{NaiveFinder}(1386:1622, '+', 3)
 ORFI{NaiveFinder}(1447:1635, '+', 1)
 ORFI{NaiveFinder}(1489:1635, '+', 1)
 ORFI{NaiveFinder}(1501:1635, '+', 1)
 ORFI{NaiveFinder}(1531:1635, '+', 1)
 ORFI{NaiveFinder}(2697:3227, '+', 3)
 ORFI{NaiveFinder}(2745:3227, '+', 3)
 ⋮
 ORFI{NaiveFinder}(2874:3227, '+', 3)
 ORFI{NaiveFinder}(2973:3227, '+', 3)
 ORFI{NaiveFinder}(3108:3227, '+', 3)
 ORFI{NaiveFinder}(3142:3312, '+', 1)
 ORFI{NaiveFinder}(3481:3939, '+', 1)
 ORFI{NaiveFinder}(3659:3934, '+', 2)
 ORFI{NaiveFinder}(3734:3934, '+', 2)
 ORFI{NaiveFinder}(3772:3939, '+', 1)
 ORFI{NaiveFinder}(3806:3934, '+', 2)
 ORFI{NaiveFinder}(4129:4287, '+', 1)
 ORFI{NaiveFinder}(4160:4291, '-', 2)
 ORFI{NaiveFinder}(4540:4644, '+', 1)
 ORFI{NaiveFinder}(4690:4866, '+', 1)
 ORFI{NaiveFinder}(4741:4866, '+', 1)
 ORFI{NaiveFinder}(4744:4866, '+', 1)
```