LDLIBS = -lgsl -lm

CC=g++

DESTDIR    ?=
prefix     ?= $(DESTDIR)/usr
bindir     ?= $(prefix)/bin

all: fit_res

fit_res: fit_res.o fit.o
fit_res.o: fit.h

install:
	mkdir -p ${bindir}
	install -m755 fit_res ${bindir}

clean:
	rm -f fit_res *.o

