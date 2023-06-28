# DNA as a Markov chain

Several packages (e.g. [MarkovChainsHammer](), [Markovians](), etc.) in
the Julia ecosystem have been implemented to work with Markov chains
with a *state space* of integers, those could be efficient in many ways,
but they are clumsy to work with a specialized biological types as in
the `BioJulia` ecosystem. Therefore, in the `GeneFinder` package we
dedicated some implementations to work with `BioSequence` types so that
we can expand the functionality in an efficient way (see complete
[API]()).

One important step towards many gene finding algorithms is to represent
a DNA sequence as a Markov chain. In this representation a DNA sequence
of a reduced alphabet 𝒜 = {*A*, *C*, *G*, *T*} is draw as a four-vertex
graph, where each letter of 𝒜 is a *state (vertex)* and the edges of the
graph represent *transitions* from one nucleotide to another in a
sequence (e.g. *A* → *T* represent a single nucleotide to nucleotide
transition). This is also considered more specifically as a Discrete
Markov chain (Axelson-Fisk 2015). The complete set of transitions and
states of a DNA sequence of alphabet 𝒜 could be seen as
<a href="#fig-dna-markov" class="quarto-xref">Figure 1</a>.

<img src="../assets/nucleotide-markov-chain.png" id="fig-dna-markov"
data-fig-align="center"
alt="Figure 1: DNA sequence as a Markov chain with a DNA alphabet \mathscr{A} = \{A, C, G, T\}" />

More formally a Markov chain is a random process where each state is a
random variable *X*<sub>*t*</sub> where *t* ∈ *T* is a discrete time in
a finite sequence *T* and the probability to jump from one state into
another is *only dependent of the current state.* Therefore a definition
of this *Markov property* is given by:

*P*(*X*<sub>*t*</sub>=*j*|*X*<sub>*t* − 1</sub>=*i*)

where *i*, *j* ∈ 𝒜 . This property led us to generalize a way to
calculate the probability of a sequence *T* from a process
*X*<sub>1</sub>...*X*<sub>*T*</sub> where each random variable is a
nucleotide from 𝒜 so that:

$$
P(X\_{1} = i\_{1},...,X\_{T} = i\_{T}) = P(X\_{1} = i\_{1}) \prod\_{t=2}^{T} P(X\_{t} = i\_{t} | X\_{t−1} = i\_{t−1})
$$

Note that previous equations has two terms, a initial probability
*P*(*X*<sub>1</sub>=*i*<sub>1</sub>) and the the product of all
transitions beginning at *t* = 2.

## Markov chain BioSequences

We can now calculate a transition matrix from a `LongDNA` sequence using
`transition_probability_matrix` method

``` julia
using BioSequences, GeneFinder

genome = randdnaseq(10^6)

cds = getcds(genome, min_len = 64)[1]

cds
```

    69nt DNA Sequence:
    ATGATTCCTGTACGTTCTAGCTGCTTCCTGGGTCTAAAACTCGTTTTCTGCCATTACTGCCACGACTGA

The `cds` object is a an ORF with the potential to encode a CDS from the
randomly generated `genome`. To see what is the transition probabilities
and the initial distribution, we can use build the Transition Model
using simply the constructor `TransitionModel`

``` julia
TransitionModel(cds)
```

    TransitionModel:
      - Transition Probability Matrix (Size: 4 × 4):
        0.25    0.417   0.083   0.25    
        0.105   0.211   0.158   0.526   
        0.231   0.308   0.154   0.308   
        0.167   0.25    0.292   0.292   
      - Initials (Size: 1 × 4):
        0.176   0.279   0.191   0.353   
      - order: 1

## References

Axelson-Fisk, Marina. 2015. *Comparative Gene Finding*. Vol. 20.
Computational Biology. London: Springer London.
<http://link.springer.com/10.1007/978-1-4471-6693-1>.
