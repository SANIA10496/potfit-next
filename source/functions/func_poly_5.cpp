/****************************************************************
 *
 * func_poly_5.cpp:
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

#include "func_poly_5.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncPOLY_5::FuncPOLY_5() : dr(0.0) {}

FuncPOLY_5::~FuncPOLY_5() {}

int FuncPOLY_5::num_params(void) {
  return 5;
}

void FuncPOLY_5::calc(const double &r, const std::vector<double> &p, double *f) {
  dr = (r - 1.) * (r - 1.);

  *f = p[0] + .5 * p[1] * dr + p[2] * (r - 1.) * dr + p[3] * (dr * dr) + p[4] * (dr * dr) * (r - 1.);

  return;
}
