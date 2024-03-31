## Scoring a sequence using a Markov model

A sequence of DNA could be scored using a Markov model of the transition probabilities of a known sequence. This could be done using a *log-odds ratio* score, which is the logarithm of the ratio of the transition probabilities of the sequence given a model and. The log-odds ratio score is defined as:

```math
\begin{align}
S(x) = \sum_{i=1}^{L} \beta_{x_{i}x} = \sum_{i=1} \log \frac{a^{\mathscr{m}_{1}}_{i-1} x_i}{a^{\mathscr{m}_{2}}_{i-1} x_i}
\end{align}
```

Where the ``a^{\mathscr{m}_{1}}_{i-1} x_i`` is the transition probability of the first model (in this case the calculated for the given sequence) from the state ``x_{i-1}`` to the state ``x_i`` and ``a^{\mathscr{m}_{2}}_{i-1} x_i`` is the transition probability of the second model from the state ``x_{i-1}`` to the state ``x_i``. The score is the sum of the log-odds ratio of the transition probabilities of the sequence given the two models.

In the current implementation the second model is a CDS transition probability model of *E. coli*. This classification score is implemented in the `naivescorefinder` method. This method will return ORFs with the associated score of the sequence given the CDS model of *E. coli*.

```julia
using GeneFinder, BioSequences

seq = dna"TTCGTCAGTCGTTCTGTTTCATTCAATACGATAGTAATGTATTTTTCGTGCATTTCCGGTGGAATCGTGCCGTCCAGCATAGCCTCCAGATATCCCCTTATAGAGGTCAGAGGGGAACGGAAATCGTGGGATACATTGGCTACAAACTTTTTCTGATCATCCTCGGAACGGGCAATTTCGCTTGCCATATAATTCAGACAGGAAGCCAGATAACCGATTTCATCCTCACTATCGACCTGAAATTCATAATGCATATTACCGGCAGCATACTGCTCTGTGGCATGAGTGATCTTCCTCAGAGGAATATATACGATCTCAGTGAAAAAGATCAGAATGATCAGGGATAGCAGGAACAGGATTGCCAGGGTGATATAGGAAATATTCAGCAGGTTGTTACAGGATTTCTGAATATCATTCATATCAGTATGGATGACTACATAGCCTTTTACCTTGTAGTTGGAGGTAATGGGAGCAAATACAGTAAGTACATCCGAATCAAAATTACCGAAGAAATCACCAACAATGTAATAGGAGCCGCTGGTTACGGTCGAATCAAAATTCTCAATGACAACCACATTCTCCACATCTAAGGGACTATTGGTATCCAGTACCAGTCGTCCGGAGGGATTGATGATGCGAATCTCGGAATTCAGGTAGACCGCCAGGGAGTCCAGCTGCATTTTAACGGTCTCCAAAGTTGTTTCACTGGTGTACAATCCGCCGGCATAGGTTCCGGCGATCAGGGTTGCTTCGGAATAGAGACTTTCTGCCTTTTCCCGGATCAGATGTTCTTTGGTCATATTGGGAACAAAAGTTGTAACAATGATGAAACCAAATACACCAAAAATAAAATATGCGAGTATAAATTTTAGATAAAGTGTTTTTTTCATAACAAATCCTGCTTTTGGTATGACTTAATTACGTACTTCGAATTTATAGCCGATGCCCCAGATGGTGCTGATCTTCCAGTTGGCATGATCCTTGATCTTCTC"

orfs = findorfs(seq, min_len=75, findermethod=naivefinderscored)
9-element Vector{ORF}:
 ORF(37:156, '+', 1, -0.0024384392084479912)
 ORF(194:268, '-', 2, -0.003991702119459298)
 ORF(194:283, '-', 2, -0.01431767026931985)
 ORF(249:347, '+', 3, -0.02024959025464718)
 ORF(426:590, '+', 3, -0.003289228147424537)
 ORF(565:657, '+', 1, -0.014806468147370438)
 ORF(650:727, '-', 2, -0.009087704913650461)
 ORF(786:872, '+', 3, -0.03486633273294755)
 ORF(887:976, '-', 2, -0.005778301450517392)
```

## The *log-odds ratio* decision rule

The sequence probability given a transition probability model could be used as the source of a sequence classification based on a decision rule to classify whether a sequence correspond to a model or another. Now, imagine we got two DNA sequence transition models, a CDS model and a No-CDS model. The *log-odds ratio* decision rule could be establish as:

``` math
\begin{align}
S(X) = \log \frac{{P_C(X_1=i_1, \ldots, X_T=i_T)}}{{P_N(X_1=i_1, \ldots, X_T=i_T)}}  \begin{cases} > \eta & \Rightarrow \text{coding} \\ < \eta & \Rightarrow \text{noncoding} \end{cases}
\end{align}
```

Where the ``P_{C}`` is the probability of the sequence given a CDS model, ``P_{N}`` is the probability of the sequence given a No-CDS model, the decision rule is finally based on whether the ratio is greater or lesser than a given threshold *η* of significance level.

In this package we have implemented this rule and call some basic models of CDS and No-CDS of *E. coli* from Axelson-Fisk (2015) work (implemented in `BioMarkovChains.jl` package). To check whether a random sequence could be coding based on these decision we use the predicate `isnaivecoding` with the `ECOLICDS` and `ECOLINOCDS` models:

```julia
orfsdna = get_orfs_dna(seq, findermethod=naivefinderscored, min_len=75, alternative_start=true);
isnaivecoding.(orfsdna)
```
```julia
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

In this case, the sequence has 20 ORFs and only 3 of them are classified as coding sequences. The classification is based on the *log-odds ratio* decision rule and the transition probability models of *E. coli* CDS and No-CDS. The `isnaivecoding` method will return a boolean vector with the classification of each ORF in the sequence.