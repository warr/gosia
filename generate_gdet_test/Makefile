GDET = det.gdt
INP  = generate_gdet.inp
OUT  = generate_gdet.out
GOOD = det.gdt.good # This is the one provided by Adam

all: compare

$(GDET): $(INP)
	gosia < $(INP)

compare: $(GDET) $(GOOD)
	diff $(GDET) $(GOOD)


clean:
	rm -f *~ $(OUT) $(GDET)
