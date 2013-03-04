SHELL := /bin/bash

default: dSFMT.o

.c.o:
	${CXX} -DDSFMT_MEXP=19937 -march=native -c $<

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il

help:
	@echo "Usage: make"

