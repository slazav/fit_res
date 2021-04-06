LDLIBS = -lgsl -lm

CC=g++

all: fit_res

fit_res: fit_res.o fit.o
fit_res.o: fit.h

bindir ?= /usr/bin
install:
	mkdir -p ${bindir}
	install -m755 fit_res ${bindir}

clean:
	rm -f fit_res *.o

