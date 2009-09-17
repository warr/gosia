TEX = internals.tex
DVI = internals.dvi
AUX = internals.aux
TOC = internals.toc
LOG = internals.log
PS  = internals.ps
PDF = internals.pdf

DEP += adiabacity.tex
DEP += collision.tex
DEP += common.tex
DEP += coupling.tex
DEP += ideas.tex
DEP += kepler.tex
DEP += overview.tex
DEP += spherical.tex
DEP += functions.tex
DEP += load.tex

all: $(PDF)

$(PDF) : $(PS)
	ps2pdf -sPAPERSIZE=a4 $(PS)

$(PS) : $(DVI)
	dvips $(DVI)
	
$(DVI) : $(TEX) $(DEP)
	latex $(TEX)
	latex $(TEX)
	latex $(TEX)

clean:
	rm -f *~ $(DVI) $(AUX) $(TOC) $(LOG) $(PS)
	
proper:
	rm -f *~ $(DVI) $(AUX) $(TOC) $(LOG) $(PS) $(PDF)
