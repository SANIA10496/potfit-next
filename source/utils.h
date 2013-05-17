/****************************************************************
 *
 * utils.h:
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

#ifndef PTF_UTILS_H
#define PTF_UTILS_H

#include "pointers.h"
#include "types.h"

namespace POTFIT_NS {

  class Utils : protected Pointers {
  public:
    Utils(class POTFIT *);
    ~Utils();

    char *tolowercase(char *);
    double vect_dist(vector, vector);
    void quicksort(double *, int, int, double **);

    char flagfile[255];
  private:
    int partition(double *, int, int, int, double **);
    void swap_population(double *, double *);
  };

}

template<class T>
T square(T a) {return a*a;}


#endif /* PTF_UTILS_H */
