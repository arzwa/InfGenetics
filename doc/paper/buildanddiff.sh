buildtex.sh paper.tex
latexdiff submission1.tex paper.tex > diff.tex
buildtex.sh diff.tex
