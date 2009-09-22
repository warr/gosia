TEX = internals.tex
DVI = internals.dvi
AUX = internals.aux
TOC = internals.toc
LOG = internals.log
PS  = internals.ps
PDF = internals.pdf

DEP += adhoc.tex
DEP += alloc.tex
DEP += ampder.tex
DEP += angula.tex
DEP += apram.tex
DEP += arccos.tex
DEP += arctg.tex
DEP += branr.tex
DEP += bricc.tex
DEP += cclkup.tex
DEP += cegry.tex
DEP += chmem.tex
DEP += cmlab.tex
DEP += code7.tex
DEP += conv.tex
DEP += coord.tex
DEP += decay.tex
DEP += djmm.tex
DEP += double.tex
DEP += effix.tex
DEP += elmt.tex
DEP += expon.tex
DEP += fakp.tex
DEP += faza1.tex
DEP += faza.tex
DEP += f.tex
DEP += fhip.tex
DEP += fiint1.tex
DEP += fiint.tex
DEP += ftbm.tex
DEP += func1.tex
DEP += func.tex
DEP += fxis1.tex
DEP += fxis2.tex
DEP += gamatt.tex
DEP += gcf.tex
DEP += gf.tex
DEP += gkk.tex
#DEP += gkvac.tex
#DEP += gosia.tex
DEP += half.tex
#DEP += intg.tex
#DEP += invkin.tex
#DEP += klopot.tex
#DEP += kontur.tex
#DEP += lagran.tex
#DEP += laiamp.tex
#DEP += laisum.tex
#DEP += leadf.tex
#DEP += limits.tex
DEP += load.tex
DEP += lsloop.tex
#DEP += mem.tex
#DEP += mini.tex
#DEP += mixr.tex
#DEP += mixup.tex
#DEP += newcat.tex
DEP += newcnv.tex
#DEP += newlv.tex
#DEP += openf.tex
#DEP += path.tex
#DEP += podziel.tex
#DEP += pol4.tex
#DEP += pomnoz.tex
#DEP += prelm.tex
#DEP += prim.tex
#DEP += pticc.tex
DEP += qe.tex
#DEP += qfit.tex
DEP += qm.tex
#DEP += qrange.tex
DEP += rangel.tex
DEP += ready.tex
#DEP += recoil.tex
#DEP += reset.tex
#DEP += rk4.tex
#DEP += rndm.tex
#DEP += rotate.tex
#DEP += select.tex
DEP += seq.tex
#DEP += setin.tex
#DEP += simin.tex
#DEP += sixel.tex
DEP += snake.tex
DEP += spline.tex
DEP += splint.tex
DEP += splner.tex
#DEP += stamp.tex
#DEP += sting.tex
#DEP += szereg.tex
#DEP += tacos.tex
#DEP += tapma.tex
#DEP += tasin.tex
#DEP += tcabs.tex
#DEP += tcexp.tex
#DEP += tenb.tex
#DEP += tens.tex
#DEP += trint.tex
DEP += wsixj.tex
DEP += wthrej.tex
DEP += xstatic.tex
#DEP += ylm1.tex
#DEP += ylm.tex

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
