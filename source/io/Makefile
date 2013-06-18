SHELL := /bin/bash

FORCE_FILES = $(shell ls *.cpp)

OBJECTS := $(subst .cpp,.o,${FORCE_FILES})

default: CFLAGS += -O3 -march=native -Wno-unused
default: ${OBJECTS}

debug: CFLAGS += -g3
debug: default

profile: CFLAGS += -pg
profile: default

# generic compilation rule
.cpp.o:
	${CXX} ${CFLAGS} -c $<

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il

help:
	@echo "Usage: make"

