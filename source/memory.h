/****************************************************************
 *
 * memory.h: misc io routines
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

#ifndef PTF_MEMORY_H
#define PTF_MEMORY_H

#include <cstdio>

#include "pointers.h"

namespace POTFIT_NS {

  class Memory : protected Pointers {
  public:
    Memory(class POTFIT *);

    void *smalloc(int , const char *);
    void *srealloc(void *, int n, const char *);
    void sfree(void *);

    // create 1D array
    template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name)
    {
      int nbytes = sizeof(TYPE) * n;
      array = (TYPE *) smalloc(nbytes,name);
      return array;
    }

    // grow/shrink 1D array
    template <typename TYPE>
    TYPE *grow(TYPE *&array, int n, const char *name)
    {
      if (array == NULL) return create(array,n,name);

      int nbytes = sizeof(TYPE) * n;
      array = (TYPE *) srealloc(array,nbytes,name);
      return array;
    }

    // destroy 1D array
    template <typename TYPE>
    void destroy(TYPE *array)
    {
      sfree(array);
    }

    // 2+3D arrays are still missing

  private:
  };

}

#endif /* PTF_MEMORY_H */
