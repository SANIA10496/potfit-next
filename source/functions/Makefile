SHELL := /bin/bash -O extglob

FORCE_FILES = $(shell ls *.cpp)

OBJECTS := $(subst .cpp,.o,${FORCE_FILES})

default: headers ${OBJECTS}

debug: CFLAGS += -g3
debug: default

headers:
	@for files in func_*.h; do echo "#include \"$${files}\""; done > list_functions.h

# generic compilation rule
.cpp.o:
	${CXX} ${CFLAGS} -c $<

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il list_functions.h

help:
	@echo "Usage: make"

