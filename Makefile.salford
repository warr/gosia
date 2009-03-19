BINDIR=$(ROOT)/usr/bin
MANDIR=$(ROOT)/usr/share/man/man1

BASE=gosia
EXE=$(BASE).exe
MAN=$(BASE).1
SRCS=$(BASE).f

FC=wine ftn77.exe
LINK=wine slink.exe

# Turn on debugging - note that -Wall and -O2 together gives warnings
# about variables being possibly used without being initialised. These
# cases seem to be harmless. i.e. if (condition) then a=1 else a=2 endif
# and then using a gives this warning, but in fact a is always set to
# something
FFLAGS = 

DEPS=Makefile

ALL: $(EXE)

OBJS += $(BASE).obj
OBJS += adhoc.obj
OBJS += alloc.obj
OBJS += ampder.obj
OBJS += angula.obj
OBJS += apram.obj
OBJS += arccos.obj
OBJS += arctg.obj
OBJS += ats.obj
OBJS += branr.obj
OBJS += bricc.obj
OBJS += cclkup.obj
OBJS += cegry.obj
OBJS += chmem.obj
OBJS += cmlab.obj
OBJS += code7.obj
OBJS += conv.obj
OBJS += coord.obj
OBJS += decay.obj
OBJS += djmm.obj
OBJS += double.obj
OBJS += effix.obj
OBJS += elmt.obj
OBJS += expon.obj
OBJS += f.obj
OBJS += fakp.obj
OBJS += faza.obj
OBJS += faza1.obj
OBJS += fhip.obj
OBJS += fiint.obj
OBJS += fiint1.obj
OBJS += ftbm.obj
OBJS += func.obj
OBJS += func1.obj
OBJS += fxis1.obj
OBJS += fxis2.obj
OBJS += gamatt.obj
OBJS += gcf.obj
OBJS += gf.obj
OBJS += gkk.obj
OBJS += gkvac.obj
OBJS += half.obj
OBJS += intg.obj
OBJS += klopot.obj
OBJS += kontur.obj
OBJS += lagran.obj
OBJS += laiamp.obj
OBJS += laisum.obj
OBJS += leadf.obj
OBJS += limits.obj
OBJS += load.obj
OBJS += lsloop.obj
OBJS += mem.obj
OBJS += mini.obj
OBJS += mixr.obj
OBJS += mixup.obj
OBJS += newcat.obj
OBJS += newcnv.obj
OBJS += newlv.obj
OBJS += openf.obj
OBJS += path.obj
OBJS += podziel.obj
OBJS += pol4.obj
OBJS += pomnoz.obj
OBJS += prelm.obj
OBJS += prim.obj
OBJS += pticc.obj
OBJS += qe.obj
OBJS += qfit.obj
OBJS += qm.obj
OBJS += qrange.obj
OBJS += rangel.obj
OBJS += ready.obj
OBJS += recoil.obj
OBJS += reset.obj
OBJS += rk4.obj
OBJS += rndm.obj
OBJS += rotate.obj
OBJS += select.obj
OBJS += seq.obj
OBJS += setin.obj
OBJS += simin.obj
OBJS += sixel.obj
OBJS += snake.obj
OBJS += spline.obj
OBJS += splint.obj
OBJS += splner.obj
OBJS += stamp.obj
OBJS += sting.obj
OBJS += szereg.obj
OBJS += tacos.obj
OBJS += tapma.obj
OBJS += tasin.obj
OBJS += tcabs.obj
OBJS += tcexp.obj
OBJS += tenb.obj
OBJS += tens.obj
OBJS += trint.obj
OBJS += wsixj.obj
OBJS += wthrej.obj
OBJS += xstatic.obj
OBJS += ylm.obj
OBJS += ylm1.obj

SRCS += arccos.f arctg.f load.f lsloop.f leadf.f mem.f cmlab.f qe.f qm.f \
snake.f fhip.f alloc.f rangel.f qrange.f ampder.f laisum.f expon.f faza.f \
setin.f sting.f laiamp.f faza1.f trint.f pol4.f stamp.f reset.f half.f \
double.f path.f intg.f newlv.f code7.f apram.f newcat.f pomnoz.f tenb.f \
tens.f djmm.f ftbm.f mini.f cegry.f fakp.f prim.f seq.f gf.f f.f conv.f \
wthrej.f wsixj.f lagran.f func.f func1.f gkvac.f gkk.f xstatic.f ats.f ylm.f \
decay.f angula.f ready.f branr.f limits.f szereg.f sixel.f prelm.f recoil.f \
rotate.f ylm1.f fiint.f fiint1.f tapma.f simin.f mixup.f fxis1.f fxis2.f \
podziel.f klopot.f mixr.f coord.f chmem.f pticc.f rndm.f kontur.f rk4.f \
qfit.f gamatt.f gcf.f tcexp.f tcabs.f tasin.f tacos.f openf.f effix.f \
adhoc.f elmt.f select.f bricc.f newcnv.f splner.f spline.f splint.f cclkup.f
	
include: include.c
	gcc -o $@ $<

DATE=$(shell date +%04Y%02m%02d)
SINGLE_FILE = $(BASE)_$(DATE).f

%.obj: %.f
	$(FC) $(FFLAGS) $*.f

$(EXE): $(OBJS) $(DEPS)
	$(LINK) $(LDFLAGS) -file:$@ $(OBJS)

clean:
	rm -f *~ *.obj $(EXE) $(BASE)_20*.f include

install: $(EXE) $(MAN)
	install -m 755 -d $(BINDIR)
	install -m 755 -s $(EXE) $(BINDIR)
	install -m 755 -d $(MANDIR)
	install -m 644 $(MAN) $(MANDIR)
	gzip -f $(MANDIR)/$(MAN)

single_file: include
	./include $(SRCS) > $(SINGLE_FILE)

