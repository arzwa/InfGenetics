buildtex.sh paper.tex
latexdiff submission2.tex paper.tex > diff.tex
buildtex.sh diff.tex
