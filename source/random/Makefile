SHELL := /bin/bash

default: CFLAGS += -O3 -march=native -Wno-unused
default: dSFMT.o

debug: CFLAGS += -g3
debug: default

profile: CFLAGS += -pg
profile: default

.c.o:
	${CXX} ${CFLAGS} -DDSFMT_MEXP=19937 -march=native -c $<

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il

help:
	@echo "Usage: make"

