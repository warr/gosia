GOSIADIR=/home/warr/devel/gosia

INP += simple-test.inp
INP += fort.9
INP += fort.22.good

OUT += fort.14
OUT += fort.15
OUT += fort.17
OUT += fort.22

all: compare

compare: fort.22
	diff fort.22 fort.22.good

$(GOSIADIR)/gosia:
	make -C $(GOSIADIR)

fort.22: $(GOSIADIR)/gosia
	$(GOSIADIR)/gosia < simple-test.inp

clean:
	rm -f *~ $(OUT)

