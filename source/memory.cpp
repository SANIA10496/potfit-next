/****************************************************************
 *
 * memory.cpp:
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include <mpi.h>
#include <cstdlib>

#include "memory.h"
#include "io.h"

using namespace POTFIT_NS;

Memory::Memory(POTFIT *ptf) : Pointers(ptf) {}

// safe memory (de-)allocation

void *Memory::smalloc(int nbytes, const char *name)
{
  if (nbytes == 0) return NULL;

  void *ptr = malloc(nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate %d bytes for array %s", nbytes,name);
    io->error(str);
  }
  return ptr;
}

void *Memory::srealloc(void *ptr, int nbytes, const char *name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return NULL;
  }

  ptr = realloc(ptr,nbytes);
  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate %d bytes for array %s",nbytes,name);
    io->error(str);
  }
  return ptr;
}

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}
