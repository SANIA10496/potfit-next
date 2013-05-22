/****************************************************************
 *
 * splines.h:
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

#ifndef PTF_SPLINES_H
#define PTF_SPLINES_H

#include <vector>

#include "../pointers.h"

namespace POTFIT_NS {
  class Splines : protected Pointers {
  public:
    Splines(class POTFIT *);
    ~Splines(void);

    // init splines with equidistand or non-equidistant sampling points
    void  spline_ed(double, double *, int, double, double, double *);
    void  spline_ne(double *, double *, int, double, double, double *);

    // calculate values with spline interpolation
//    double splint_ed(pot_table_t *, double *, int, double);
//    double splint_ne(pot_table_t *, double *, int, double);

    // calculate gradients with spline interpolation
//    double splint_grad_ed(pot_table_t *, double *, int, double);
//    double splint_grad_ne(pot_table_t *, double *, int, double);

    // calculate values and gradients with spline interpolation
//    double splint_comb_ed(pot_table_t *, double *, int, double, double *);
//    double splint_comb_ne(pot_table_t *, double *, int, double, double *);

    // calculate ???
    double splint_dir(double *, double *, int, double, double);
    double splint_comb_dir(double *, double *, int, double, double, double *);
//    double splint_grad_dir(pot_table_t *, double *, int, double, double);

    std::vector<double> u;
  };
}

#endif // PTF_SPLINES_H
