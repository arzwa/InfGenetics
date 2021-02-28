---
title: The infinitesimal model for polyploid and mixed-ploidy populations
author: Arthur Zwaenepoel
fontfamily: stix
fontsize: 12pt
---

\renewcommand{\Pr}{\mathbb{P}}
\newcommand{\Ex}{\mathbb{E}}
\newcommand{\var}{\mathrm{var}}
\newcommand{\No}{\mathcal{N}}

# The infinitesimal model

The infinitesimal model assumes that offspring trait values $z_{ij}$ of a
parental pair with trait values $z_i$ and $z_j$ is normally distributed:

$$ z_{ij} \sim \mathcal{N}\Big(\frac{z_i + z_j}{2}, V\Big) $$

where $V$ is called the *segregation variance* and is determined by the
hereditary process. The basic infinitesimal model can be derived as the limit
of a model where a quantitative trait is controlled by a large number $n$ of
Mendelian loci with additive gene action, each of small effect $\sim
O(\sqrt{n})$, as $n$ gets large. 

Importantly, $V$ is not a function of $z_i$ or $z_j$, but can evolve over
time. If populations are finite, inbreeding will cause the segregation
variance to decrease as a function of the inbreeding coefficients of the
parental individuals.

While the basic phenotypic model holds for arbitrary ploidy levels, the
evolution of the segregation variance in finite populations differs for
different ploidy levels.

# Haploids and diploids

# Tetraploids without double reduction

## The segregation variance

The expected segregation variance for a family with parents $i$ and $j$ of
identical ploidy levels can be decomposed into a contribution from both parents

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
transmitted by $i$ to its offspring are not IBD. The total segregation variance
for the parental pair $(i,j)$ will be

$$\Ex[V_{ij}] = V_0\bigg(\frac{G_i + G_j}{2}\bigg)$$

That is, in autotetraploids that do not undergo double reduction, the
infinitesimal model (as a limit of Mendelian additive loci) under inbreeding is
defined as the model where offspring trait values of parents $i$ and $j$
follows a Gaussian distribution:

$$z_{ij} \sim \No\Bigg(\frac{z_i + z_j}{2}, 
    V_0\bigg(\frac{G_i + G_j}{2}\bigg)\Bigg)$$

Where $V_0$ is the segregation variance in the base population and where $G_i$
is the probability that the two genes at a locus transmitted by $i$ are not IBD.
Note that $G_i$ is a function of the inbreeding coefficients of the *parents*
of $i$.

## Recursions for the inbreeding coefficients

\begin{align*}
F_{i,i}(t) &= \sum_k \sum_l P_{i,k} P_{i,l} \frac{1}{6}\big(F_{k,k}(t-1) + 
    F_{l,l}(t-1) + 4F_{k,l}(t-1)\big) \\
F_{i,j}(t) &= \sum_k \sum_l P_{i,k} P_{j,l} F_{k,l}^\ast(t-1) 
\end{align*}

where 

$$ F_{k,l}^\ast(t) = \begin{cases} F_{k,l}(t) & \text{if } k \ne l \\
    \frac{1}{4}\big(1 + 3F_{k,k}(t-1)\big) & \text{if } k = l \end{cases} $$
    
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
