SHELL := /bin/bash

FORCE_FILES = $(shell ls *.cpp)

OBJECTS := $(subst .cpp,.o,${FORCE_FILES})

default: headers ${OBJECTS}

headers:
	@for files in force_*.h; do echo "#include \"$${files}\""; done > list_forces.h

# generic compilation rule
.cpp.o:
	${CXX} ${CFLAGS} -c $<

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il list_forces.h

help:
	@echo "Usage: make"

