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
#
# Beware: This Makefile works only with GNU make (gmake)!
#
# Usage:  make
#
###########################################################################

###########################################################################
#
#  Adjust these variables to your system
#
###########################################################################

# Currently the following systems are available:
# x86_64-icc  	64bit Intel Compiler
# x86_64-gcc    64bit GNU Compiler
# i686-icc 	32bit Intel Compiler
# i686-gcc  	32bit GNU Compiler
#
#SYSTEM 		= x86_64-icc 	# Use this as fallback
SYSTEM 		= $(shell uname -m)-gcc

# Base directory of your installation of the MKL or ACML
MKLDIR          = /common/linux/paket/intel/compiler-11.0/cc/mkl
ACML4DIR  	= /common/linux/paket/acml4.4.0/ifort64
ACML5DIR  	= /opt/acml5.2.0/ifort64
LIBMDIR 	= /opt/amdlibm

###########################################################################
#
#  Defaults for some variables
#
###########################################################################

STRIP 		= $(shell which strip 2> /dev/null)

LIBS		+= -lm -lstdc++
DEBUG_FLAGS	+= -DDEBUG
MKLPATH         = ${MKLDIR}/lib
ACML4PATH 	= ${ACML4DIR}/lib
ACML5PATH 	= ${ACML5DIR}/lib

###########################################################################
#
#  flags for 64bit
#
###########################################################################

ifeq (x86_64-icc,${SYSTEM})
# compiler
  OMPI_CXX      = icc
  OMPI_CLINKER  = icc

# general optimization flags
  OPT_FLAGS     += -fast -xHost

# profiling and debug flags
  PROF_FLAGS    += -prof-gen
  PROF_LIBS     += -prof-gen
  DEBUG_FLAGS   += -g -Wall

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/em64t
  CINCLUDE 	+= -I${MKLDIR}/include
  LIBS 		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_lp64_sequential.a \
		   -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential \
		   -lmkl_core -Wl,--end-group -lpthread
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE 	+= -I${ACML4DIR}/include
  LIBS		= -L${ACML4PATH} -lpthread -lacml -lacml_mv
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
   LIBMPATH 	= ${LIBMDIR}/lib/static
   CINCLUDE     += -I${ACML5DIR}/include -I${LIBMDIR}/include
   LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm
endif

 export        OMPI_CXX OMPI_CLINKER
endif

ifeq (x86_64-gcc,${SYSTEM})
# compiler
  OMPI_CXX      = gcc
  OMPI_CXXLINKER  = gcc

# general optimization flags
  OPT_FLAGS     += -O3 -march=native -Wno-unused

# profiling and debug flags
  PROF_FLAGS    += -g3 -pg
  PROF_LIBS     += -g3 -pg
  DEBUG_FLAGS   += -g3 -Wall

# Intel Math Kernel Library
# ifeq (,$(strip $(findstring acml,${MAKETARGET})))
#   MKLPATH       = ${MKLDIR}/lib/em64t/
#   CINCLUDE      += -I${MKLDIR}/include
#   LIBS 		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_lp64_sequential.a \
#                    -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
#                    -Wl,--end-group -lpthread -Wl,--as-needed
# endif

# AMD Core Math Library
# ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
#   CINCLUDE     		+= -I${ACML4DIR}/include
#   LIBS			+= -L${ACML4PATH} -lpthread -lacml -lacml_mv -Wl,--as-needed
# endif
# ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
#   LIBMPATH 	= ${LIBMDIR}/lib/static
#   CINCLUDE     	+= -I${ACML5DIR}/include -I${LIBMDIR}/include
#   LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm -Wl,--as-needed
# endif

 export        OMPI_CXX OMPI_CXXLINKER
endif


###########################################################################
#
#  flags for 32bit
#
###########################################################################

ifeq (i686-icc,${SYSTEM})
# compiler
  OMPI_CXX      = icc
  OMPI_CLINKER  = icc

# general optimization flags
  OPT_FLAGS	+= -fast -xHost

# profiling and debug flags
  PROF_FLAGS	+= -prof-gen
  PROF_LIBS 	+= -prof-gen
  DEBUG_FLAGS	+= -g -Wall

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/32
  CINCLUDE      += -I${MKLDIR}/include
  LIBS 		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_sequential.a \
		   -Wl,--start-group -lmkl_intel -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I$(ACML4DIR)/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml -lacml_mv
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/static
  CINCLUDE     	+= -I$(ACML5DIR)/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm
endif

  export        OMPI_CXX OMPI_CLINKER
endif

ifeq (i686-gcc,${SYSTEM})
# compiler
  OMPI_CXX     	= gcc
  OMPI_CLINKER 	= gcc

# general optimization flags
  OPT_FLAGS	+= -O3 -march=native -Wno-unused

# profiling and debug flags
  PROF_FLAGS	+= -g3 -pg
  PROF_LIBS	+= -g3 -pg
  DEBUG_FLAGS	+= -g3 -Wall

# Intel Math Kernel Library
ifeq (,$(strip $(findstring acml,${MAKETARGET})))
  MKLPATH       = ${MKLDIR}/lib/32
  CINCLUDE      += -I${MKLDIR}/include
  LIBS		+= -L${MKLPATH} ${MKLPATH}/libmkl_solver_sequential.a \
		   -Wl,--start-group -lmkl_intel -lmkl_sequential -lmkl_core \
		   -Wl,--end-group -lpthread -Wl,--as-needed
endif

# AMD Core Math Library
ifneq (,$(strip $(findstring acml4,${MAKETARGET})))
  CINCLUDE     	+= -I$(ACML4DIR)/include
  LIBS		+= -L${ACML4PATH} -lpthread -lacml -lacml_mv -Wl,--as-needed
endif
ifneq (,$(strip $(findstring acml5,${MAKETARGET})))
  LIBMPATH 	= ${LIBMDIR}/lib/static
  CINCLUDE     	+= -I$(ACML5DIR)/include -I${LIBMDIR}/include
  LIBS		+= -L${ACML5PATH} -L${LIBMPATH} -lpthread -lacml -lamdlibm -Wl,--as-needed
endif

  export        OMPI_CXX OMPI_CLINKER
endif

###########################################################################
#
#  Compiler, flags, libraries
#
###########################################################################

# compiler; if empty, we issue an error later
CXX = mpicxx

# optimization flags
OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS} -DNDEBUG
DEBUG_FLAGS += ${${PARALLEL}_FLAGS} ${DEBUG_${PARALLEL}_FLAGS}

# optimization or debug
CFLAGS := ${FLAGS}

default: git potfit strip

debug: CFLAGS += ${DEBUG_FLAGS}
debug: TARGET = debug
debug: potfit
potfit: CFLAGS += ${OPT_FLAGS}

SUBDIRS = forces functions random tables
CLEANDIRS = $(SUBDIRS:%=clean-%)

.PHONY: subdirs ${SUBDIRS}

subdirs: ${SUBDIRS}

${SUBDIRS}:
	${MAKE} -C $@ ${TARGET}

SOURCES = $(shell ls *.cpp)

FORCE_SOURCES = $(shell ls forces/*.cpp)
FORCE_OBJECTS = $(subst .cpp,.o,${FORCE_SOURCES})
FUNCTION_SOURCES = $(shell ls functions/*.cpp)
FUNCTION_OBJECTS = $(subst .cpp,.o,${FUNCTION_SOURCES})
RANDOM_SOURCES = $(shell ls random/*.c)
RANDOM_OBJECTS = $(subst .c,.o,${RANDOM_SOURCES})
TABLE_SOURCES = $(shell ls tables/*.cpp)
TABLE_OBJECTS = $(subst .cpp,.o,${TABLE_SOURCES})

# Substitute .o for .c to get the names of the object files
OBJECTS := $(subst .cpp,.o,${SOURCES})

git:
	@echo -e "#define POTFIT_VERSION \"git (r`git rev-list HEAD | wc -l`)\"" > version.h
	@echo -e "#define POTFIT_DATE \"`date +%Y-%m-%d\ %H:%M:%S\ %z`\"" >> version.h

# potfit is the default target
potfit: subdirs ${OBJECTS}
	${CXX} ${LIBS} ${LFLAGS_${PARALLEL}} -o $@ ${OBJECTS} ${FORCE_OBJECTS} ${RANDOM_OBJECTS} ${TABLE_OBJECTS} ${FUNCTION_OBJECTS}

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

# generic compilation rule
.cpp.o:
	${CXX} ${CFLAGS} -c $<

clean: ${CLEANDIRS}
	rm -f *.o *.u *~ \#* *.V *.T *.O *.il potfit

${CLEANDIRS}:
	$(MAKE) -C $(@:clean-%=%) clean

help:
	@echo "Usage: make"

