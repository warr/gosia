BINDIR=$(ROOT)/usr/bin
MANDIR=$(ROOT)/usr/share/man/man1

UID=0
GID=0

EXE=gosia
MAN=gosia.1

FFLAGS = -g -Wall

DEPS=Makefile

ALL: $(EXE)

OBJS += gosia.o

gosia.o: gosia.f $(DEPS)
	g77 $(FFLAGS) -c gosia.f

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

