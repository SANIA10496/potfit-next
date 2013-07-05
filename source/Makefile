############################################################################
#
# potfit -- The ITAP Force Matching Program
# 	Copyright 2002-2013
#
# 	Institute for Theoretical and Applied Physics,
# 	University of Stuttgart, D-70550 Stuttgart, Germany
# 	http://potfit.itap.physik.uni-stuttgart.de/
#
############################################################################
#
#   This file is part of potfit.
#
#   potfit is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   potfit is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with potfit; if not, see <http://www.gnu.org/licenses/>.
#
############################################################################

include Makefile.in

# optimization flags
# OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS} -DNDEBUG
# DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}

# optimization or debug
CFLAGS := ${FLAGS}

default: git potfit strip

debug: CFLAGS += ${DEBUG_FLAGS} -DDEBUG
debug: TARGET = debug
debug: git potfit
#potfit: CFLAGS += ${OPT_FLAGS}

profile: CFLAGS += ${PROF_FLAGS}
profile: LIBS += -pg
profile: TARGET = profile
profile: git potfit

SUBDIRS = config forces functions io opt random tables
CLEANDIRS = $(SUBDIRS:%=clean-%)

.PHONY: subdirs ${SUBDIRS}

subdirs: ${SUBDIRS}

${SUBDIRS}:
	${MAKE} -C $@ ${TARGET}

SOURCES = $(shell ls *.cpp)

CONFIG_SOURCES = $(shell ls config/*.cpp)
CONFIG_OBJECTS = $(subst .cpp,.o,${CONFIG_SOURCES})
FORCE_SOURCES = $(shell ls forces/*.cpp)
FORCE_OBJECTS = $(subst .cpp,.o,${FORCE_SOURCES})
FUNCTION_SOURCES = $(shell ls functions/*.cpp)
FUNCTION_OBJECTS = $(subst .cpp,.o,${FUNCTION_SOURCES})
IO_SOURCES = $(shell ls io/*.cpp)
IO_OBJECTS = $(subst .cpp,.o,${IO_SOURCES})
OPT_SOURCES = $(shell ls opt/*.cpp)
OPT_OBJECTS = $(subst .cpp,.o,${OPT_SOURCES})
RANDOM_SOURCES = $(shell ls random/*.c)
RANDOM_OBJECTS = $(subst .c,.o,${RANDOM_SOURCES})
TABLE_SOURCES = $(shell ls tables/*.cpp)
TABLE_OBJECTS = $(subst .cpp,.o,${TABLE_SOURCES})

OBJECTS := $(subst .cpp,.o,${SOURCES})

git:
	@echo -e "#define POTFIT_VERSION \"git (r`git rev-list HEAD | wc -l`)\"" > version.h
	@echo -e "#define POTFIT_DATE \"`date +%Y-%m-%d\ %H:%M:%S\ %z`\"" >> version.h

# potfit is the default target
potfit: subdirs ${OBJECTS}
	${CXX} ${LIBS} ${LFLAGS_${PARALLEL}} -o $@ ${OBJECTS} ${FORCE_OBJECTS} ${RANDOM_OBJECTS} ${TABLE_OBJECTS} ${FUNCTION_OBJECTS} ${CONFIG_OBJECTS} ${OPT_OBJECTS} ${IO_OBJECTS} -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -Wl,--as-needed

strip:
	${STRIP} --strip-unneeded -R .comment potfit

RANDOM_HEADERS:
	${MAKE} -C random $1

FORCE_HEADERS:
	${MAKE} -C forces

TABLE_HEADERS:
	${MAKE} -C tables

# all objects depend on headers
#${OBJECTS}: ${HEADERS}

utils.o: utils.cpp
	${CXX} ${OPT_FLAGS} ${CFLAGS} -I/common/linux/paket/intel/compiler-11.0/cc/mkl/include -c utils.cpp

# generic compilation rule
.cpp.o:
	${CXX} ${OPT_FLAGS} ${CFLAGS} -c $<

clean: ${CLEANDIRS}
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il potfit version.h

${CLEANDIRS}:
	$(MAKE) -C $(@:clean-%=%) clean

help:
	@echo "Usage: make"

