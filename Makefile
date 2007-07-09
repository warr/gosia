BINDIR=$(ROOT)/usr/bin
MANDIR=$(ROOT)/usr/share/man/man1

UID=0
GID=0

EXE=gosia
MAN=gosia.1

FFLAGS = -g -Wall

DEPS=Makefile

ALL: $(EXE)

OBJS += adhoc.o
OBJS += alloc.o
OBJS += ampder.o
OBJS += angula.o
OBJS += apram.o
OBJS += arccos.o
OBJS += arctg.o
OBJS += ats.o
OBJS += branr.o
OBJS += cegry.o
OBJS += chmem.o
OBJS += cmlab.o
OBJS += code7.o
OBJS += conv.o
OBJS += coord.o
OBJS += decay.o
OBJS += djmm.o
OBJS += double.o
OBJS += effix.o
OBJS += elmt.o
OBJS += expon.o
OBJS += f.o
OBJS += fakp.o
OBJS += faza.o
OBJS += faza1.o
OBJS += fhip.o
OBJS += fiint.o
OBJS += fiint1.o
OBJS += ftbm.o
OBJS += func.o
OBJS += func1.o
OBJS += fxis1.o
OBJS += fxis2.o
OBJS += gamatt.o
OBJS += gcf.o
OBJS += gf.o
OBJS += gkk.o
OBJS += gkvac.o
OBJS += gosia.o
OBJS += half.o
OBJS += intg.o
OBJS += klopot.o
OBJS += kontur.o
OBJS += lagran.o
OBJS += laiamp.o
OBJS += laisum.o
OBJS += leadf.o
OBJS += limits.o
OBJS += load.o
OBJS += lsloop.o
OBJS += mem.o
OBJS += mini.o
OBJS += mixr.o
OBJS += mixup.o
OBJS += newcat.o
OBJS += newlv.o
OBJS += openf.o
OBJS += path.o
OBJS += podziel.o
OBJS += pol4.o
OBJS += pomnoz.o
OBJS += prelm.o
OBJS += prim.o
OBJS += pticc.o
OBJS += qe.o
OBJS += qfit.o
OBJS += qm.o
OBJS += qrange.o
OBJS += rangel.o
OBJS += ready.o
OBJS += recoil.o
OBJS += reset.o
OBJS += rk4.o
OBJS += rndm.o
OBJS += rotate.o
OBJS += seq.o
OBJS += setin.o
OBJS += simin.o
OBJS += sixel.o
OBJS += snake.o
OBJS += stamp.o
OBJS += sting.o
OBJS += szereg.o
OBJS += tacos.o
OBJS += tapma.o
OBJS += tasin.o
OBJS += tcabs.o
OBJS += tcexp.o
OBJS += tenb.o
OBJS += tens.o
OBJS += trint.o
OBJS += wsixj.o
OBJS += wthrej.o
OBJS += xstatic.o
OBJS += ylm.o
OBJS += ylm1.o

adhoc.o: adhoc.f $(DEPS)
	g77 $(FFLAGS) -c adhoc.f

alloc.o: alloc.f $(DEPS)
	g77 $(FFLAGS) -c alloc.f

ampder.o: ampder.f $(DEPS)
	g77 $(FFLAGS) -c ampder.f

angula.o: angula.f $(DEPS)
	g77 $(FFLAGS) -c angula.f

apram.o: apram.f $(DEPS)
	g77 $(FFLAGS) -c apram.f

arccos.o: arccos.f $(DEPS)
	g77 $(FFLAGS) -c arccos.f

arctg.o: arctg.f $(DEPS)
	g77 $(FFLAGS) -c arctg.f

ats.o: ats.f $(DEPS)
	g77 $(FFLAGS) -c ats.f

branr.o: branr.f $(DEPS)
	g77 $(FFLAGS) -c branr.f

cegry.o: cegry.f $(DEPS)
	g77 $(FFLAGS) -c cegry.f

chmem.o: chmem.f $(DEPS)
	g77 $(FFLAGS) -c chmem.f

cmlab.o: cmlab.f $(DEPS)
	g77 $(FFLAGS) -c cmlab.f

code7.o: code7.f $(DEPS)
	g77 $(FFLAGS) -c code7.f

conv.o: conv.f $(DEPS)
	g77 $(FFLAGS) -c conv.f

coord.o: coord.f $(DEPS)
	g77 $(FFLAGS) -c coord.f

decay.o: decay.f $(DEPS)
	g77 $(FFLAGS) -c decay.f

djmm.o: djmm.f $(DEPS)
	g77 $(FFLAGS) -c djmm.f

double.o: double.f $(DEPS)
	g77 $(FFLAGS) -c double.f

effix.o: effix.f $(DEPS)
	g77 $(FFLAGS) -c effix.f

elmt.o: elmt.f $(DEPS)
	g77 $(FFLAGS) -c elmt.f

expon.o: expon.f $(DEPS)
	g77 $(FFLAGS) -c expon.f

f.o: f.f $(DEPS)
	g77 $(FFLAGS) -c f.f

fakp.o: fakp.f $(DEPS)
	g77 $(FFLAGS) -c fakp.f

faza.o: faza.f $(DEPS)
	g77 $(FFLAGS) -c faza.f

faza1.o: faza1.f $(DEPS)
	g77 $(FFLAGS) -c faza1.f

fhip.o: fhip.f $(DEPS)
	g77 $(FFLAGS) -c fhip.f

fiint.o: fiint.f $(DEPS)
	g77 $(FFLAGS) -c fiint.f

fiint1.o: fiint1.f $(DEPS)
	g77 $(FFLAGS) -c fiint1.f

ftbm.o: ftbm.f $(DEPS)
	g77 $(FFLAGS) -c ftbm.f

func.o: func.f $(DEPS)
	g77 $(FFLAGS) -c func.f

func1.o: func1.f $(DEPS)
	g77 $(FFLAGS) -c func1.f

fxis1.o: fxis1.f $(DEPS)
	g77 $(FFLAGS) -c fxis1.f

fxis2.o: fxis2.f $(DEPS)
	g77 $(FFLAGS) -c fxis2.f

gamatt.o: gamatt.f $(DEPS)
	g77 $(FFLAGS) -c gamatt.f

gcf.o: gcf.f $(DEPS)
	g77 $(FFLAGS) -c gcf.f

gf.o: gf.f $(DEPS)
	g77 $(FFLAGS) -c gf.f

gkk.o: gkk.f $(DEPS)
	g77 $(FFLAGS) -c gkk.f

gkvac.o: gkvac.f $(DEPS)
	g77 $(FFLAGS) -c gkvac.f

gosia.o: gosia.f $(DEPS)
	g77 $(FFLAGS) -c gosia.f

half.o: half.f $(DEPS)
	g77 $(FFLAGS) -c half.f

intg.o: intg.f $(DEPS)
	g77 $(FFLAGS) -c intg.f

klopot.o: klopot.f $(DEPS)
	g77 $(FFLAGS) -c klopot.f

kontur.o: kontur.f $(DEPS)
	g77 $(FFLAGS) -c kontur.f

lagran.o: lagran.f $(DEPS)
	g77 $(FFLAGS) -c lagran.f

laiamp.o: laiamp.f $(DEPS)
	g77 $(FFLAGS) -c laiamp.f

laisum.o: laisum.f $(DEPS)
	g77 $(FFLAGS) -c laisum.f

leadf.o: leadf.f $(DEPS)
	g77 $(FFLAGS) -c leadf.f

limits.o: limits.f $(DEPS)
	g77 $(FFLAGS) -c limits.f

load.o: load.f $(DEPS)
	g77 $(FFLAGS) -c load.f

lsloop.o: lsloop.f $(DEPS)
	g77 $(FFLAGS) -c lsloop.f

mem.o: mem.f $(DEPS)
	g77 $(FFLAGS) -c mem.f

mini.o: mini.f $(DEPS)
	g77 $(FFLAGS) -c mini.f

mixr.o: mixr.f $(DEPS)
	g77 $(FFLAGS) -c mixr.f

mixup.o: mixup.f $(DEPS)
	g77 $(FFLAGS) -c mixup.f

newcat.o: newcat.f $(DEPS)
	g77 $(FFLAGS) -c newcat.f

newlv.o: newlv.f $(DEPS)
	g77 $(FFLAGS) -c newlv.f

openf.o: openf.f $(DEPS)
	g77 $(FFLAGS) -c openf.f

path.o: path.f $(DEPS)
	g77 $(FFLAGS) -c path.f

podziel.o: podziel.f $(DEPS)
	g77 $(FFLAGS) -c podziel.f

pol4.o: pol4.f $(DEPS)
	g77 $(FFLAGS) -c pol4.f

pomnoz.o: pomnoz.f $(DEPS)
	g77 $(FFLAGS) -c pomnoz.f

prelm.o: prelm.f $(DEPS)
	g77 $(FFLAGS) -c prelm.f

prim.o: prim.f $(DEPS)
	g77 $(FFLAGS) -c prim.f

pticc.o: pticc.f $(DEPS)
	g77 $(FFLAGS) -c pticc.f

qe.o: qe.f $(DEPS)
	g77 $(FFLAGS) -c qe.f

qfit.o: qfit.f $(DEPS)
	g77 $(FFLAGS) -c qfit.f

qm.o: qm.f $(DEPS)
	g77 $(FFLAGS) -c qm.f

qrange.o: qrange.f $(DEPS)
	g77 $(FFLAGS) -c qrange.f

rangel.o: rangel.f $(DEPS)
	g77 $(FFLAGS) -c rangel.f

ready.o: ready.f $(DEPS)
	g77 $(FFLAGS) -c ready.f

recoil.o: recoil.f $(DEPS)
	g77 $(FFLAGS) -c recoil.f

reset.o: reset.f $(DEPS)
	g77 $(FFLAGS) -c reset.f

rk4.o: rk4.f $(DEPS)
	g77 $(FFLAGS) -c rk4.f

rndm.o: rndm.f $(DEPS)
	g77 $(FFLAGS) -c rndm.f

rotate.o: rotate.f $(DEPS)
	g77 $(FFLAGS) -c rotate.f

seq.o: seq.f $(DEPS)
	g77 $(FFLAGS) -c seq.f

setin.o: setin.f $(DEPS)
	g77 $(FFLAGS) -c setin.f

simin.o: simin.f $(DEPS)
	g77 $(FFLAGS) -c simin.f

sixel.o: sixel.f $(DEPS)
	g77 $(FFLAGS) -c sixel.f

snake.o: snake.f $(DEPS)
	g77 $(FFLAGS) -c snake.f

stamp.o: stamp.f $(DEPS)
	g77 $(FFLAGS) -c stamp.f

sting.o: sting.f $(DEPS)
	g77 $(FFLAGS) -c sting.f

szereg.o: szereg.f $(DEPS)
	g77 $(FFLAGS) -c szereg.f

tacos.o: tacos.f $(DEPS)
	g77 $(FFLAGS) -c tacos.f

tapma.o: tapma.f $(DEPS)
	g77 $(FFLAGS) -c tapma.f

tasin.o: tasin.f $(DEPS)
	g77 $(FFLAGS) -c tasin.f

tcabs.o: tcabs.f $(DEPS)
	g77 $(FFLAGS) -c tcabs.f

tcexp.o: tcexp.f $(DEPS)
	g77 $(FFLAGS) -c tcexp.f

tenb.o: tenb.f $(DEPS)
	g77 $(FFLAGS) -c tenb.f

tens.o: tens.f $(DEPS)
	g77 $(FFLAGS) -c tens.f

trint.o: trint.f $(DEPS)
	g77 $(FFLAGS) -c trint.f

wsixj.o: wsixj.f $(DEPS)
	g77 $(FFLAGS) -c wsixj.f

wthrej.o: wthrej.f $(DEPS)
	g77 $(FFLAGS) -c wthrej.f

xstatic.o: xstatic.f $(DEPS)
	g77 $(FFLAGS) -c xstatic.f

ylm.o: ylm.f $(DEPS)
	g77 $(FFLAGS) -c ylm.f

ylm1.o: ylm1.f $(DEPS)
	g77 $(FFLAGS) -c ylm1.f


gosia: $(OBJS) $(DEPS)
	g77 $(LDFLAGS) -o gosia $(OBJS)

clean:
	rm -f *~ *.o $(EXE)

install: $(EXE) $(MAN)
	install -m 755 -o $(UID) -g $(GID) -d $(BINDIR)
	install -m 755 -o $(UID) -g $(GID) -s $(EXE) $(BINDIR)
	install -m 755 -o $(UID) -g $(GID) -d $(MANDIR)
	install -m 644 -o $(UID) -g $(GID) $(MAN) $(MANDIR)
	gzip -f $(MANDIR)/$(MAN)

