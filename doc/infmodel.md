---
title: The infinitesimal model for polyploid and mixed-ploidy populations
author: Arthur Zwaenepoel
fontfamily: stix
fontsize: 12pt
header-includes: | 
 \usepackage{tikz}
 \usetikzlibrary{positioning}
 \usetikzlibrary{arrows}
---

\renewcommand{\Pr}{\mathbb{P}}
\newcommand{\Ex}{\mathbb{E}}
\newcommand{\var}{\mathrm{var}}
\newcommand{\No}{\mathcal{N}}


# The basic infinitesimal model

The infinitesimal model assumes that offspring trait values $z_{ij}$ of a
parental pair with trait values $z_i$ and $z_j$ are normally distributed:

$$ z_{ij} \sim \mathcal{N}\Big(\frac{z_i + z_j}{2}, V_{ij}\Big) $$

where $V_{ij}$ is called the *segregation variance* and is determined by the
hereditary process. The basic infinitesimal model can be derived as the limit
of a model where a quantitative trait is controlled by a large number $n$ of
Mendelian loci with additive gene action, each of small effect $\sim
O(\sqrt{n})$. 

Importantly, $V_{ij}$ is not a function of $z_i$ or $z_j$, but can evolve over
time. If populations are finite, inbreeding will cause the segregation
variance to decrease as a function of the inbreeding coefficients of the
parental individuals. 

A slightly different, and perhaps more insightful, way to specify the same
model is to write $z_{ij} = X_i + X_j$, where $X_i$ is the contribution of
parent $i$ to the genotypic value of the offspring and $X_j$ the same for $j$.
That is, $X_i$ is the genotypic value of a gamete from $i$.  For gametes
produced by a normal meiotic division, we assume $X_i \sim \No(z_i/2, V_i)$,
where $V_i$ is the contribution to the segregation variance from $i$. We
therefore have $V_{ij} = (V_i + V_j)/2$.  This way of formulating the model
stresses that segregation occurs independently in both parents, contributing
additively to the segregation variance (which is the variance among offspring
within a family).

While the basic phenotypic model holds for arbitrary ploidy levels, the
evolution of the segregation variance in finite populations differs for
different ploidy levels. In addition, when considering mixed-ploidy
populations, we need to consider how equilibrium variances in the model scale
with ploidy level, and how offspring of different ploidy levels are derived
from a parental pair.

## Haploids and diploids

For haploids the recursion for the IBD coefficients in terms of the pedigree
matrix is

$$F_{i,j}' = \begin{cases}
    \sum_k \sum_l P_{i,k} P_{j,l} F_{k,l} & \text{if } i \ne j \\
    0 & \text{if } i = j
    \end{cases} $$
    
For diploids, it should be:   
    
$$F_{i,j}' = \sum_k \sum_l P_{i,k} P_{j,l} \begin{cases} 
    F_{k,l} & \text{if } k \ne l \\
    \frac{1}{2}(1 + F_{k,k}) & \text{if } k = l
    \end{cases} $$
    
Although I wonder whether this is correct for $F_{i,i}'$, since this is defined
as the probability of IBD of distinct homologs in $i$. Consider an individual
$i$ offspring from $k$ and $l$, then there will be a term $P_{i,k}P_{i,k}
\frac{1}{2}(1 + F_{k,k}) = \frac{1}{8}(1+F_{k,k})$ in the sum for $F_{i,i}$,
which makes no sense?

## Tetraploids without double reduction

### The segregation variance

Under inbreeding, the expected segregation variance for a family with parents
$i$ and $j$ of identical ploidy levels can be decomposed into a contribution
from both parents

$$ \Ex[V_{ij}] = \Ex[V_i] + \Ex[V_j] $$

We will determine $\Ex[V_i]$. We define $V_0$ to be the segregation variance in
the base population consisting of unrelated individuals, so that $\Ex[V_i] =
V_0/2$ in the absence of inbreeding.

Assume the parents of $i$ in the previous generation were $k$ and $l$ and
consider the contribution of a single locus to $\Ex[V_i]$. We can consider
three mutually exclusive patterns of ancestry:

1. $i$ transmits both genes it inherited from $k$ (w.p. $1/6$)
2. $i$ transmits both genes it inherited from $l$ (w.p. $1/6$)
3. $i$ transmits two genes inherited from distinct parents (w.p. $2/3$)

When the two genes transmitted are IBD, they do not contribute to the 
segregation variance. Therefore, the expected contribution to the segregation
variance for a single locus $\Ex[v_i]$ conditional on, for instance, scenario 1
above would be 

$$\Ex[v_i] = 0 \times F_{k,k}(t-1) + (v_0/2) (1 - F_{k,k}(t-1))$$

Where $v_0 = V_0/n$ is the average per-locus segregation variance in the base
population. Combining all cases, we get

\begin{align*}
\Ex[V_i] &= \frac{1}{6} \frac{V_0}{2}\big(1 - F_{k,k}(t-1)\big) + 
        \frac{1}{6} \frac{V_0}{2} \big(1-F_{l,l}(t-1)\big) + 
        \frac{2}{3}\frac{V_0}{2}\big(1 - F_{k,l}(t-1)\big) \\
    &= \frac{V_0}{2} \bigg( 1 - \frac{1}{6}\big(F_{k,k}(t-1) + 
        F_{l,l}(t-1) + 4F_{k,l}(t-1)\big)\bigg) \\
    &\equiv \frac{V_0}{2} G_i
\end{align*}

Note that $G_i$ as defined here is the probability that the two genes
transmitted by $i$ to its offspring are not IBD. If $i$ had parents $k$
and $j$:

$$G_i = 1 - \frac{1}{6}\big(F_{k,k} + F_{l,l} + 4F_{k,l}\Big)$$

Where the inbreeding coefficients are from the generation of the parents of
$i$. The total segregation variance for the parental pair $(i,j)$ will be

$$\Ex[V_{ij}] = V_0\bigg(\frac{G_i + G_j}{2}\bigg)$$

That is, in autotetraploids that do not undergo double reduction, the
infinitesimal model (as a limit of Mendelian additive loci) under inbreeding is
defined as the model where offspring trait values of parents $i$ and $j$
follows a Gaussian distribution:

$$z_{ij} \sim \No\Bigg(\frac{z_i + z_j}{2}, 
    V_0\bigg(\frac{G_i + G_j}{2}\bigg)\Bigg)$$

Where $V_0$ is the segregation variance in the base population and where $G_i$
is the probability that the two genes at a locus transmitted by $i$ are not
IBD.  Note again that $G_i$ is a function of the inbreeding coefficients of the
*parents* of $i$, so that in the tetraploid model, the inbreeding coefficients
of *grandparents* are needed for computing the segregation variance.

### Recursions for the inbreeding coefficients

\begin{align*}
F_{i,i}(t) &= \sum_k \sum_l P_{i,k} P_{i,l} \frac{1}{6}\big(F_{k,k}(t-1) + 
    F_{l,l}(t-1) + 4F_{k,l}(t-1)\big) \\
F_{i,j}(t) &= \sum_k \sum_l P_{i,k} P_{j,l} F_{k,l}^\ast(t-1) 
\end{align*}

where 

$$ F_{k,l}^\ast(t) = \begin{cases} F_{k,l}(t) & \text{if } k \ne l \\
    \frac{1}{4}\big(1 + 3F_{k,k}(t)\big) & \text{if } k = l \end{cases} $$
    
To see the latter, note that $k = l$ means that $i$ and $j$ share a parent, in
which case for a given gene in $i$, there is a $1/2$ chance it came from te
shared parent, in which case there is a $1/2$ chance that $j$ inherited the
same gene independently and there is a $1/4$ chance that this gene is picked in
$j$, so that the probability a given gene is IBD with a random gene in $j$ is
$1/16 = P_{i,k}P_{j,k}(1/4)$. When the genes trace back to distinct lineages in
the same parent, with probability $P_{i,k}P_{j,k}(3/4)$, they are IBD with
probability $F_{k,k}(t-1)$. Note that the probability that both genes trace
back to the same parent is embodied in the pedigree matrix $P$, while
$F_{k,k}^\ast$ is the probability of being IBD *conditional* on tracing to the
same parent.

The recursions for the inbreeding coefficients and the law for the reduction of
segregation variance with inbreeding together enable simulation of the
infinitesimal model for tetraploid populations (without double reduction).
I verified the correctness of the derived recursions by comparing simulations
from the infinitesimal model to individual-based simulations of a random mating
tetraploid population with a large number ($n=1000$) of unlinked additive
Mendelian loci of small effect ($\sim O(1/\sqrt{4n}$).

![Comparison of the evolution of the segregation variance under inbreeding for
the infinitesimal model (black line) and the 1000 unlinked additive loci model
(grey). A population of 50 individuals is simulated over 1000 generations.  In
red the approximation $V_0 e^{-t/4N}$ is shown.](img/tet-segvar.pdf)

## Double reduction

When tetravalents are formed, a form of internal inbreeding may occur as a
result of the phenomenan called double reduction. Schematically, double 
reduction for a genotype $ABCD$ could look like:

    <a scheme>

**I don't think this is right**
Assume double reduction occurs at any locus with probability $d$ (what
*exactly* happens with probability $d$?). When a single double reduction takes
place at a locus, helf of the gametes will contain a pair of IBD genes at the
locus, so that the segregation variance becomes $\frac{v_0}{2}$. Following this
reasoning, the segregation variance can be found as

$$\Ex[V_i] = (1-d)^2 \Ex[V_i^\ast] + 2d(1-d) \frac{\Ex[V_i^\ast]}{2}$$

where $\Ex[V_i^\ast]$ is the contribution to the segregation variance of
individual $i$ in the absence of double reduction. This verifies that double
reduction leads to a kind of 'internal' inbreeding, as the segregation variance
is reduced with respect to $V_i^\ast$.

# The infinitesimal model for mixed-ploidy populations

## Scaling of genetic variance across ploidy levels

If we consider the infinitesimal model as the limit of a large number of
additive Mendelian loci, the genotypic value of an individual is defined as

$$ z = \sum_{k=1}^n \sum_{l=1}^m a_{k,l} $$

where $a_{k,m}$ is the allelic effect of homolog $l$ of locus $k$. Clearly, if
we assume identical allelic effects in diploids and congeneric polyploids, the
phenotypic range of tetraploids will be double that of diploids. While the
assumption of additive allelic effects is of course in itself problematic,
conditional on this assumption, the assumption of *equal* additive effects
across ploidy levels seems rather unbiological.

We introduce a scaler for allelic effects in an $m$-ploid $\sqrt{\beta_m}$, so that

$$ z = \sqrt{\beta_m} \sum_{k=1}^n \sum_{l=1}^m a_{k,l} $$

where we define $\beta_2 = 1$, so that $a_{k\cdot}$ is the effect the allele
would have in a diploid individual. With this model the HWLE additive genetic
variance $V_{A,m}$ and the base segregation variance $V_{0,m}$ in an $m$-ploid
population are $\frac{m}{2}\beta_m V_{A,2}$ and $\frac{m}{2}\beta_m V_{0,2}$
respectively.

## Unreduced gamete formation in diploids

Naively, one may think that an unreduced gamete contains the parental genome,
and that as a result the segregation variance for a $2n \times 2n \rightarrow
4n$ cross would be 0. However, the mechanisms of unreduced gamete formation do
not necessarily lead to a faithful transmission of the complete diploid genome.
Unreduced gametes are formed in two ways, depending on the meiotic abberation
that leads to their origin: (1) first division restitution (FDR) of (2) second
division restitution (SDR). Consider a locus in a diploid with two distinct
genes $A$ and $a$. Assume recombination happens with probability $c$ and
that conditional on unreduced gamete formation, formation is due to FDR with
probability $f$ while it is due to SDR with probability $1-f$. The different
unreduced gametes that are formed are represented in the following diagram:

\begin{center}
\begin{tikzpicture}
\node (rep) at (-6, 4.5) [] {Replication};
\node (rep) at (-6, 3.5) [] {Recombination};
\node (rep) at (-6, 2.5) [] {First division};
\node (rep) at (-6, 1.5) [] {Second division};
\node (rep) at (-6, 0) [] {Gametes};
\node (Aa)  at (0, 5) [] {$Aa$};
\node (AAaa)  at (0, 4) [] {$AA\ aa$};
\node (rAaAa) at (-2, 3) [] {$Aa\ Aa$};
\node (rAAaa) at ( 2, 3) [] {$AA\ aa$};
\path (Aa) edge (AAaa);
\path (AAaa) edge node[near end, above] {$c$} (rAaAa) ;
\path (AAaa) edge node[near end, above] {$1-c$} (rAAaa);

\node (fdr11) at (-3, 2) [] {$Aa\ Aa | \emptyset$};
\node (fdr12) at (-3, 1) [] {$A|a\  A|a$}; 
\node (AA) at (-3.8,0) [] {$AA$};
\node (p1) at (-3.8,-0.5) [] {$\frac{1}{4}$};
\node (Aa) at (-3,0) [] {$Aa$};
\node (p2) at (-3,-0.5) [] {$\frac{1}{2}$};
\node (aa) at (-2.2,0) [] {$aa$};
\node (p2) at (-2.2,-0.5) [] {$\frac{1}{4}$};
\path (fdr12) edge (AA);
\path (fdr12) edge (Aa);
\path (fdr12) edge (aa);
\path (fdr11) edge (fdr12);
\path (rAaAa) edge node[left] {$f$} (fdr11);

\node (sdr11) at (-1, 2) [] {$Aa | Aa$}; 
\node (sdr12) at (-1, 1) [] {$Aa | \emptyset$}; 
\node (Aa2) at (-1,0) [] {$Aa$};
\path (sdr12) edge (Aa2);
\path (sdr11) edge (sdr12);
\path (rAaAa) edge node[right] {$1-f$} (sdr11);

\node (fdr21) at (1,2) [] {$AA\ aa|\emptyset$};
\node (fdr22) at (1,1) [] {$A|A\ a|a$};
\node (Aa3) at (1,0) [] {$Aa$};
\path (rAAaa) edge node[left] {$f$} (fdr21);
\path (fdr21) edge (fdr22);
\path (fdr22) edge (Aa3);

\node (sdr21) at (3, 2) [] {$AA | aa$};
\node (sdr22) at (3, 1) [] {$AA| \emptyset, aa|\emptyset$};
\node (AA4) at (2.2,0) [] {$AA$};
\node (p1) at (2.2, -0.5) [] {$\frac{1}{2}$};
\node (aa4) at (3.8,0) [] {$aa$};
\node (p2) at (3.8, -0.5) [] {$\frac{1}{2}$};
\path (rAAaa) edge node[right] {$1-f$} (sdr21);
\path (sdr21) edge (sdr22);
\path (sdr22) edge (AA4);
\path (sdr22) edge (aa4);
\end{tikzpicture}
\end{center}

Clearly, the mechanism of unreduced gamete formation creates segregation
variance, as not all random tetraploid offspring from a single diploid parental
pair will receive the same pair of genes from each parent, depending on whether
or not recombination has occurred and FDR rather than SDR generates the
unreduced gamete. 

If we consider the diploid genotype $Aa$ depicted in the diagram above and consider
$X_A$ the number of $A$ alleles transmitted, we can easily find

$$\var X_A = \Big(1 - c - f + \frac{3}{2}cf \Big) := \xi$$

which is also the probability of transmitting two copies of the same allele.
We shall denote this derived quantity as $\xi$.

Of course, if the locus is IBD with probability $F_{i,i}$, we have that the
expected value of $\var X_A$ is $\xi(1-F_{i,i})$.  Under the
model delineated above with the scaling of allelic effects across ploidy
levels, we find that the contribution of a single parent $i$ to the segregation
variance among tetraploid offspring of a diploid parental pair is

$$ V_i = 2\beta_4 V_{0,2} \xi (1 - F_{i,i})$$ 

Notably, the mean offspring phenotype among tetraploid offspring from a diploid
cross can no longer be the average of the two parental phenotypes if we assume
the infinitesimal model as the limit of the Mendelian unlinked additive loci
model. The expected contribution of a single parent $i$ will be $X_i =
\sqrt{\beta_4}z_i$ so that for a cross of two diploids via unreduced gametes

$$z_{ij} \sim \No\Big(\sqrt{\beta_4}(z_i + z_j), V_i + V_j\Big)$$

## Inbreeding coefficients in the mixed-ploidy system

In a mixed-ploidy system tracking inbreeding coefficients becomes slightly more
complicated, as our recursions will differ whether some individual is derived
from parents of the same cytotype or not.

Recall that $\xi$ is the probability that an unreduced gamete transmits two
copies of the same allele at some locus. We still have that for a tetraploid
individual $i$

$$ F_{i,i}' = \sum_k \sum_l P_{i,k}P_{i,l} 
    \frac{1}{6}\big(F^\ast_{k,k} + F^\ast_{l,l} + 4F_{k,l}\big) $$

where 

$$ F_{k,k}^\ast = \begin{cases}
    F_{k,k} & \text{if } m_k = 4 \\ 
    F_{k,k} + \xi & \text{if } m_k = 2 
    \end{cases} $$

For $F_{i,j}'$ we still have 

$$ F_{i,j}' = \sum_k \sum_l P_{i,k}P_{i,l} F^\ast_{k,l} $$

where now

$$ F_{k,l}^\ast = \begin{cases}
    F_{k,l} & \text{if } k \ne l \\ 
    \frac{1}{2}(1 + F_{k,k}) & \text{if } k = l, m_k = 2 \\ 
    \frac{1}{4}(1 + 3F_{k,k}) & \text{if } k = l, m_k = 4 
    \end{cases}$$

(I think)

Lastly we need $G_i$, the probability that the two genes transmitted by a
tetraploid individual $i$ to an offspring are not IBD.

$$ G_i = 1 - \frac{1}{6}\Big(F_{k,k}^\ast + F_{l,l}^\ast + 4F_{k,l}\Big) $$

where 

$$ F_{k,k} = \begin{cases} 
        F_{k,k}(1-\xi) + \xi & \text{if } m_k = 2 \\
        F_{k,k} & \text{if } m_k = 4 
    \end{cases}$$

I think that's it (ignoring triploids and double reduction for now). Not sure
how to organize the calculations yet though.
