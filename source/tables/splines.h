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

    double splint_dir(double *, double *, int, double, double);
    double splint_comb_dir(double *, double *, int, const double &, const double &, double *);
  private:
    double a, p1, p2, d21, d22;
    std::vector<double> u;
  };

  inline double Splines::splint_dir(double *xi, double *d2tab, int k, double b, double step) {
    a = 1.0 - b;
    p1 = xi[k];
    d21 = d2tab[k++];
    p2 = xi[k];
    d22 = d2tab[k];

    return a * p1 + b * p2 + ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
  }

  inline double Splines::splint_comb_dir(double *xi, double *d2tab, int k, const double &b, const double &step, double *grad) {
    a = 1.0 - b;
    p1 = xi[k];
    d21 = d2tab[k++];
    p2 = xi[k];
    d22 = d2tab[k];
    *grad = (p2 - p1) / step + ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) * step / 6.0;

    return a * p1 + b * p2 + ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
  }
}

#endif // PTF_SPLINES_H
