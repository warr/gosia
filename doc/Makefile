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
DEP += gamma_ray_attenuation.tex

IMG += geometry_for_attenuation.eps

all: $(PDF)
	
$(PDF) : $(TEX) $(DEP) $(IMG)
	pdflatex $(TEX) && pdflatex $(TEX) && pdflatex $(TEX)

clean:
	rm -f *~ $(DVI) $(AUX) $(TOC) $(LOG) $(PS) $(PDF) *-eps-converted-to.pdf

tidy:
	rm -f *~ $(DVI) $(AUX) $(TOC) $(LOG) $(PS) *-eps-converted-to.pdf
