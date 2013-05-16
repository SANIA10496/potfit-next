SHELL := /bin/bash

TABLE_FILES = $(shell ls *.cpp)

OBJECTS := $(subst .cpp,.o,${TABLE_FILES})

default: headers ${OBJECTS}
debug: CFLAGS += -g3
debug: default

headers:
	@for files in opt_*.h; do echo "#include \"$${files}\""; done > list_opt.h

# generic compilation rule
.cpp.o:
	${CXX} ${CFLAGS} -c $<

clean:
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il list_opt.h

help:
	@echo "Usage: make"

