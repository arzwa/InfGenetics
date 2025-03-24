
## Associate Editor


>Thank you for submitting your article to Evolution. It has now received two
reviews, both of which recognize its valuable contributions and essentially
recommend publication after relatively minor revision. I share these views.
While both reviews are positive, they make a number of worthwhile suggestions
for you to consider for revising the manuscript. One reviewer has suggested
combining the Results and the Discussion into a single section to avoid
unnecessary discussion of the results in the current Results section.
Personally, I would favor retaining both sections, **moving more discussive text
from the Results to the Discussion**. I leave the decision to you, but I do agree
with the reviewers that some resorting of the Results text is needed. I also
agree that the **Introduction could be better motivated by referring more to
empirical work**, and the **Discussion section could relate the results more to
other theoretical findings** in the literature. Overall, think that revision will
be relatively straightforward, but I would ask that you consider all the
suggestions made by the reviewers and, if you disagree with them, that you
justify and explain your reasons in a response letter.

I have now revised the Introduction to include more discussion of empirical
work, although I also, for the sake of brevity, refer the reader to the article
by Griswold (2021) to this end, given that he provides a nice overview relevant
studies. **TODO!**

Regarding the Results and the Discussion sections: I do not want to do away
with Discussion section altogether (as suggested by Reviewer 1), given that I
wish to discuss factors relevant to the problem which are not studied in my
article (dominance, inbreeding depression, *etc.*) and want to relate my work
to other recent theoretical work.
My aim is to describe and interpret the results in the Results section, and to
discuss the limitations and links with other work in the Discussion section.
I have sought to move 'more discussive text' to the Discussion section,
although I found there was fairly limited scope to do so when I seek to keep
the interpretation of the results more or less self-contained. **TODO!**

I have expanded the discussion slightly to include more of the previous
theoretical literature on the subject, in particular following the suggestions
of Reviewer 2. **TODO!**

## Reviewer 1

>This interesting manuscript investigated how the infinitesimal model can be
expanded to polyploidy populations (triploids and tetraploids here). It
proposes a mathematical formulation of a polyploid population's expected
genetic diversity for a given amount of inbreeding. The model is then used to
study how polyploidy can affect the adaptation to a new habitat (based on a
sink-source model).
I found the manuscript to be interesting and only have a few things to say
about it. I will for sure motivate further expansion of the model to understand
the evolution of autopolyploidy.
Having no line number does not help to comment precisely on the manuscript; I
did my best to be as precise as possible.

Thank you for these encouraging words, and my apologies for the oversight
when it comes to line numbers!

>Results are already well discussed when presented. I would either merge the
results and discussion sections (as in brief communications for example),
because it reads well as presented, or remove any explanations/discussion from
the results to put them in the discussion section. The first option sounds
better to me.

I appreciate the comment and have reworked the text with this in mind.
However, following the suggestion of the associate editor, I went with the
second option. **TODO**

>A small concluding section where the authors mentioned other questions to
study with their model will be useful in my opinion. The author proposed to
complexify the current model of adaptation with dominance or to relax some
hypotheses, but I think other biological scenarios could be studied too.

**TODO**

>Defining establishment when N=100: Is it a classic way to define establishment
in such models?  Or does it come from your experience with the model and you
see that Under N=100 extinction remains likely? More details would be needed
here.

$N=100$ is indeed an arbitrary choice. This threshold was also adopted in
@barton2018, so it makes comparisons of establishment probabilities *etc.* with
their paper somewhat more straightforward.
$N=100$ seems to be a reasonable threshold in that extinction becomes
sufficiently unlikely, while we do not have to simulate large populations (our
approach requires tracking the $N \times N$ matrix of identity coefficients). 
I did a check by putting the threshold at $N=150$, and the baseline predictions
(as in fig. 1) do not change significantly.

>Abstract: I am not sure the reference to Barton et al. (2017) is appropriate
here. Maybe a more general definition would fit better.

The reason for the reference is that their appears to be some confusion about
different models that bear the name 'infinitesimal model' (as outlined in for
instance @turelli2018 and @walsh2018). I wanted to make it clear that we use
the "Gaussian descendants" infinitesimal model, which is very carefully
outlined in @barton2017. However, I guess this is sufficiently clear from the
introduction and the methods sections, so the reference in the abstract can
indeed be omitted, as it is in the revised manuscript.

>Introduction: P2 (lower inbreeding depression (ronfort 1999, Otto & Whitton
2000): Some empirical data?  Husband et al. (2008), Clo & Kolar (2022).

I have now added the suggested references.

>Model and methods: P3-P4 (mixed ploidy population model): maybe write that the
rate of reduced gametes is fixed and genetically inherited from diploids to
tetraploids. It is clear from the equations but I suspect it could be clearer
for non-theoretical people who will not necessarily focus that much on
equations or matrices.

I added a sentence to emphasize this. I note that this is also emphasized in
the discussion "In addition, we have assumed a relatively high and constant
rate of unreduced gamete formation $u$ and triploid fertility $v$ in all our
simulations (5%), whereas these are known to be variable across the population,
and at least in part genetically determined (Kreiner et al., 2017a; Clo et al.,
2022)"

>P4 (tetraploids are not twice as big as diploids): Some references would be
needed to justify (Porturas et al. 2019 for example).

I have added the reference to Porturas et al. (2019).

>Discussion about inbreeding depression: again, some empirical data are
available to support some statements (Husband et al. (2008), Clo & Kolar
(2022)).

**TODO**

## Reviewer 2

>This paper is an important foundational work in the area of theoretical polyploid
evolutionary genetics. It presents an infinitesimal model of trait evolution for
triploids and tetraploids, and then applies the model to assess conditions under
which a tetraploid population will establish in a peripheral region connected to
a central core by dispersal. In terms of the infinitesimal model, several aspects
are very helpful
- Demonstrating how genetic effects are scaled across triploids and tetraploids
such that segregating variance is equivalent across ploidies, including with diploids
- How to combine the genetic effects in an individual with mixed-ploidy ancestry
- Calculating inbreeding and co-ancestry values for mixed-ploidy individuals
- Calculating the effect of unreduced gamete formation on segregation variance
across diploids, triploids and tetraploids
Other informative analyses include recursive or Markov models of null expecta-
tions for frequencies of ploidies, times to diploid ancestry and effective popula-
tion size across a population of mixed-ploidy.
Besides this set of results, the paper then examines conditions under which a
tetraploid population is established in a peripheral population connected to a
central core population that is predominantly composed of diploids. Given the
additive nature of gene-action in the infinitesimal model, this analysis will likely
be of interest to empiricists because additive effects are thought to be the least
constrained basis for something to evolve and there is empirical evidence that
the genetic basis of polyploid establishment in some species is additive. This is
in contrast to the theoretical results of Griswold (2021), which highlighted the
potential for a recessive-basis for autotetraploid establishment in a peripheral
region.
Main comments focus on mutation and adding background/context.

I thank the reviewer for this appraisal.

>Mutation:
For low migration rates, the possibility that autotetraploid relative to diploid es-
tablishment is driven by differences in how inbreeding accumulates is interesting
(figure 4), but the model does not include mutation. New mutation dissociates
identity by descent from identity in state and likely slows the rate of loss of
segregation variance. It would be helpful to address how scaled mutation may
affect dynamics of segregating variance between diploids and autotetraploids
and whether mutation may affect results. Is there a relationship between μVm
(Barton et al. 2017) and m and V on establishment?
Recognizing the point seems to be important, but its analysis can be for later.

I have not considered mutation at all so far because I assumed it would not be
relevant for the population sizes and timescales considered in my study.
Specifically, any individual at the time of establishment derives from a
completely outbred migrant individual a relatively short time in the past (10
to 100 generations, say), so the contribution of new mutation to differences in
establishment probability between diploids and tetraploids should be
negligible.
This is definitely the case when there is ongoing migration: as long as $m$ is
sufficiently greater than $\mu$ and $V_m$ is sufficiently smaller than the
genetic variance in the mainland population, we will have that for every
mutation contributing a very slight decrease in the loss of segregation
variance on the island, there will be many more migrant arrivals, which each
introduce a completely unrelated genotype.

Of course, when considering the mixed-ploidy infinitesimal model independently
of the establishment question, mutation does become a relevant topic. I mention
this in the revised discussion but I think it is out of scope to discuss this
in any detail in the paper. I think that my approach for scaling across ploidy
levels that is explicitly based on underlying allelic effects in a finite-locus
model extends readily to mutation. The main difference between ploidy levels is
of course that the mutational target size is proportional to the ploidy level,
so that, unless the mutational variance is exactly halved in tetraploids, the
contribution of mutation to genetic variance ($\mu V_m$) in tetraploids is
larger than in diploids.

>In the intro providing empirical examples of polyploid establishment in
peripheral habitats would broaden readership and connections to the
literature. Also, reviewing empirical work on the genetics of adaptation to
peripheral environments in polyploids would broaden connections to the
literature.

I agree with the relevance of empirical work on polyploid establishment and
adaptation to marginal habitats for the introduction of the paper. However, I
think that to present an actual review of empirical studies here is a bit out
of scope. Griswold presents a brief review in his 2021 paper, and I prefer to
refer to his paper instead of paraphrasing his paragraph. 

Similarly, I think that reviewing empirical work on adaptation to peripheral
habitats is somewhat out of scope (note that the paper is already at the word
limit for Evolution research articles). I refer to an important (if somewhat
old) review article of Kawecki (2008) and recent theoretical work by Sachdeva
et al. I think this gives a good sample of where to look for more references.

>The discussion could relate findings to other theoretical works more. For ex-
ample, it does not go back and place the results in the context of Oswald and
Nuismer (2011), which assumed a two-locus model and varied dosage effects, as
well as competition between diploids and tetraploids, including no competition.
In addition it allowed for assortative mating.

**TODO** (Discussion)

>It would be informative to expand on the the consequences of different scalings
of genetic effects. This paper made the sensible choice to scale such that segre-
gation variance is equal between ploidies as a baseline. Griswold (2021) scaled
genetic effects such that autotetraploid and diploid individuals have equal fit-
ness as a function of proportional allele count in an individual. A consequence
of Griswold’s scaling is autotetraploids have lower additive genetic variance at
HWE (Gallais 2003, p. 185), so are expected to respond more slowly in the
absence of another process. Gallais (2003, p. 186) mentions some studies that
compared additive genetic variance between diploids and autotetraploids.

To some extent, the consequences of different scalings is taken up in the
results and discussion throughout the paper where relevant. In the methods
section the relationship between allelic effects and additive genetic variance
is stated explicitly. In the results shown in Fig. 2, the effect of scaling
assumptions on establishment probability is explicitly shown and described. In
the section on establishment with recurrent migration, the phenomenon that
polyploid offspring is more extreme on average when allelic effects are not
exactly scaled by one half (as in Griswold) is taken up. 
**TODO** (Discussion) **I now added some further discussion, referring
explicitly to Griswold and Gallais in the Discussion section. I now also refer
to Porturas et al.** 

>Neither the intro nor the discussion revisit what seems to be an important point
in Barton & Etheridge (2018, p. 111) about local adaptation via infinitely small
effects.

I assume the reviewer is referring to the following issue: classical single
locus population genetics predicts that one needs $m < s$ for selection to
maintain a locally beneficial allele when migrants introduce the deleterious
allele, and that when adaptation is due to alleles of very small effect this
may seem to imply that there is little scope for local adaptation.

There are multiple reasons why this is not the case. The argument outlined in
@barton2018 is simply that the classical theory assumes complete divergence
between mainland and island, but that when adaptation is from standing
variation instead this doe snot hold, since allele frequency differences
between mainland and island are only slight. In the latter case, many alleles
of small effect can be divergently maintained by selection. Another reason why
this does not hold is linkage disequilibrium between locally beneficial
alleles, this is explored in detail in @sachdeva2022 and @zwaenepoel2024.

I do not think, however, that it is pertinent to discuss these issues in the
present paper -- they are somewhat technical and require quite some space to
explain and do not have much to do with the problem studied, i.e. polyploid
establishment. 


### Other comments:

>p.1, Abstract: Use the word “found”, but is “establish” more accurate? Overall,
the corresponding sentence is difficult to understand/parse. Suggest revision.

This sentence has been revised.

>p. 1: “negative frequency dependence” -> Isn’t MCE a positive-frequency de-
pendent process?

Indeed, it is! This has been corrected.

>p. 2: Suggest citing Bulmer (1971).

I am unsure whether this is a good idea, as it may lead to some confusion as to
which version of the infinitesimal model we are working with.
The model described in Bulmer (1971) is not exactly equivalent to the Gaussian
descendants approximation that I use (as outlined in Turelli (2017), which I
cite).
If the reference is suggested by the reviewer in order to cite a more
'original' reference, I don't think there is a reason to single out Bulmer
(1971), one could refer to Fisher (1918) or work by Alan Robertson in the
1960's.
I think the reference to Turelli (2017) which reviews the history of the
infinitesimal model suffices, and prevents confusion as to what we mean when we
refer to 'the infinitesimal model'.

>p.7: Suggest adding “, across gene copy arrangements and ploidies” at the end
of the last sentence before the “Establishment” section.

I do not fully understand this suggestion. I personally find this makes the
sentence more confusing, as I never talk of gene copy arrangements.

>p. 7, Establishment section: The organization of paragraphs two and three in
this section is a bit confusing. N is used before it is formally defined. In my
first reading I thought there was migration and then selection, whereas there is
selection and then migration. It would also help support a self-contained paper
to give the expression for E[eγ(...) ].

As far as I can tell, $N$ is defined when first mentioned. I slightly altered
the sentence to make this more clear.

I added the expression for the expected fitness (eq. 11).

The life cycle is indeed as follows: migrants arrive, they mate randomly and
offspring survives with a probability that is determined by its fitness. This
is indeed 'selection before migration'. I tried to make this more clear.

>p. 11, figure 3 (and later figure 5): The smallest value of m is 0.1 in these
figures, but based on figure 4, it looks like m = 0.01 is a transition point from
processes shaped by the accumulation of inbreeding versus migration becoming
more of a factor. Can m be extended to 0.01?

I extended the relevant figures as suggested, as $m=0.01$ is perhaps more
representative of weak migration than $m=0.1$
Weak migration means essentially that $m^2$ is negligible, so migrants and
early generation descendants thereof do not mate with migrants and the whole
process amounts more or less to independent establishment trials (there is no
swamping).

I do not think, however, that fig 4 suggests such a transition point though. As
noted in the main text (**line numbers**), the difference between the $m=0$ and
$m=0.01$ results in fig. 4 are mostly due to the fact that the latter assumes a
Gaussian distribution of migrant trait values, whereas the former assumes a
fixed initial migrant trait value.

>p. 12: Add qualifier “as m increases” at end of last sentence in last full para-
graph.

Done.

>p. 14, Figure 5 caption: “migration and selfing”

Fixed.

>p. 4, suppl., figure S4: “formation”

Fixed.

>p. 6, suppl., equation 6: Should the square only be in the denominator?

That's right, thank you!

>p. 10, suppl.: “coancestry coefficients in (10) is “ - missing the “(10)”

Fixed, thank you!

>p. 12, suppl.: add “between centromere and locus” after “. . . recombination
happens with probability c”

Done.

>p. 13-15, suppl.: Some notation seemed to not be defined, e.g. VS(2,2), V0,3, ...

The notation is introduced on p.11 in the beginning of S2.6. I have slightly
reorganized the first paragraphs of S2.6 to make the notation and general
approach more clear.


