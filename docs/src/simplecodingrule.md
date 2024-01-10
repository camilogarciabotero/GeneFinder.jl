## The *log-odds ratio* decision rule

The sequence probability given a transition probability model (eq. 2) could be used as the source of a sequence classification based on a decision rule to classify whether a sequence correspond to a model or another. Now, imagine we got two DNA sequence transition models, a CDS model and a No-CDS model. The *log-odds ratio* decision rule could be establish as:

``` math
\begin{align}
S(X) = \log \frac{{P_C(X_1=i_1, \ldots, X_T=i_T)}}{{P_N(X_1=i_1, \ldots, X_T=i_T)}}  \begin{cases} > \eta & \Rightarrow \text{coding} \\ < \eta & \Rightarrow \text{noncoding} \end{cases}
\end{align}
```

Where the ``P_{C}`` is the probability of the sequence given a
CDS model, ``P_{N}`` is the probability of the sequence given a
No-CDS model, the decision rule is finally based on whether the ratio is
greater or lesser than a given threshold *η* of significance level.

In the GeneFinder we have implemented this rule and a couple of basic
transition probability models of CDS and No-CDS of *E. coli* from
Axelson-Fisk (2015) work. To check whether a random sequence could be
coding based on these decision we use the predicate `iscoding` with the
`ECOLICDS` and `ECOLINOCDS` models:

``` julia
randseq = get_orfs_dna(randdnaseq(99))[1] # this will retrieved a random coding ORF

iscoding(randseq, ECOLICDS, ECOLINOCDS)
```

    true