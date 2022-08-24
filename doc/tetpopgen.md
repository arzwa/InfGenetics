---
geometry: a5paper, margin=1cm
fontfamily: stix
---

## Inbreeding coefficients in diploids

In diploids, one can speak of *the* inbreeding coefficient $F$ of an
individual, defined as the probability that two distinct genes at a locus are
identical by descent.
This probability is relevant, in that the pedigree of an individual with
respect to some *base generation* is sufficient to compute it.
Note that it is a probability, and hence a measure of uncertainty or degree of
belief.
Specifically, $F$ measures the degree of belief we have that two gene copies
are identical by descent *on the supposition of Mendelian segregation*.

The inbreeding coefficient in diploids is well-studied.
For various population models, e.g. the WF model, its expected value has across
the population has a well-defined law.
Note that the inbreeding coefficient of each *individual* in such a population
is still determined by the pedigree and not such a law.
However, to track the inbreeding coefficients of a diploid population, one does
not need to track the pedigree explicitly, as the latter can be adequately
represented by a matrix of coancestry coefficients (I think that is the right
term, but I am not completely sure).
This is for instance used in @barton2017 (where they term the matrix of these
coefficients the pedigree matrix).
Let $F_{i,j}$ be the coancestry coefficient for individual $i$ and $j$, i.e.
the probability that a random gene pair from $i$ and $j$ at a homologous locus
is IBD. 
$F_{i,i}$ is the usual inbreeding coefficient for individual $i$.
Let $P$ be the pedigree matrix, with entry $P_{i,k}$ be the probability that a
gene from individual $i$ came from individual $k$ in the *previous* generation.
$$F_{i,j}' = \sum_{k}\sum_{l} P_{i,k}P_{j,l} F_{k,l}$$
and
$$F_{i,i}' = \begin{cases}
    F_{k,l} & k \ne l\\
    \frac{1}{2}(1 + F_{k,k}) & k = l
    \end{cases}$$
where $k$ and $l$ are the parents of $i$.[^error]

[^error]: @barton2017 gave for diploids
$$F_{i,j}' = \sum_k \sum_l P_{i,k} P_{j,l} \begin{cases} 
    F_{k,l} & \text{if } k \ne l \\
    \frac{1}{2}(1 + F_{k,k}) & \text{if } k = l
    \end{cases}$$
However, this appears to be wrong.
Consider an individual $i$ which is an offspring from $k$ and $l$,
with $k \ne l$. There will be a term $P_{i,k}^2 \frac{1}{2}(1 + F_{k,k}) =
\frac{1}{8}(1+F_{k,k})$ as well as a term $P_{i,l}^2(1+F_{l,l})$ in the
sum for $F_{i,i}$, both of which are spurious since $F_{i,i}$ is the probability
that two distinct genes are IBD, and the probability that two *distinct* genes 
come from parent $k$ is not $P_{i,k}^2$ but 0. 

## Segregation variance in diploids

Consider a large number $n$ of diploid additive loci in linkage equilibrium,
each contributing $v$ to the segregation variance $V = \sum_{i=1}^n v$.
Segregation happens independently in each parent, with the gametic variance
for each $V/2$, and for each locus individually $v/2$.
When the inbreeding coefficient of an individual is $F$, the segregation
variance is reduced as $(1-F)v/2$ for each locus, and the total contribution of
each locus is (EVE law):
$$(1-F_i)v/2 + (1-F_j)v/2 = (1 - F_i + 1 - F_j)v/2 = (1 - (F_i + F_j)/2)v$$
and the sum over loci gives the expression of @barton2017.


## Inbreeding coefficients in tetraploids

In tetraploids there are five different degrees of homozygosity possible at any
locus. Symbolically: $abcd$, $aabc$, $aabb$, $aaab$, $aaaa$ (what's the law?).
There are, concomitantly, four inbreeding coefficients[^hexaploids].

[^hexaploids]: In hexaploids, we have 11 states, and 10 coefficients:
$abcdef$, $aabcde$, $aabbcd$, $aabbcc$, $aaabcd$, $aaabbc$, $aaabbb$, $aaaabc$,
$aaaabb$, $aaaaab$, $aaaaaa$.

## Segregation variance in autotetraploids

The segregation variance, here defined as the variance of *gametic values*, in
autotetraploids conditional on any of the five IBD states is given in the
following table

| state  | probability | segregation variance |
| -----  | ----------: | -------------------: |
| $abcd$ | $f_1$       | $V_0$                |
| $aabc$ | $f_2$       | $(1-\frac{1}{6})V_0$ |
| $aabb$ | $f_3$       | $(1-\frac{1}{3})V_0$ |
| $aaab$ | $f_4$       | $(1-\frac{1}{2})V_0$ |
| $aaaa$ | $f_5$       | $0$                  |

This is easily verified by Monte Carlo experiments.

When double reduction happens with probability $\alpha$, the segregation
variance is increased by a factor $(1+2\alpha)$ in my experiments.
However, it should be $(1+\alpha)$:
$(1-\alpha) V[X + X] + \alpha V[2X] = (1-\alpha) v_0 + 4 \alpha \frac{v_0}{2} =
v_0 (1 + \alpha)$? (this is for the total variance!)
Note that double reduction accelerates inbreeding, and hence the decline of the
segregation variance over time, but increases the segregation variance within
any given generation with respect to the segregation variance in the absence of
double reduction.

