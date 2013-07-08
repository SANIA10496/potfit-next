/****************************************************************
 *
 * func_morse.cpp:
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

#include "func_morse.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncMORSE::FuncMORSE() {}

FuncMORSE::~FuncMORSE() {}

int FuncMORSE::num_params(void) {
  return 3;
}

void FuncMORSE::calc(const double &r, const std::vector<double> &p, double *f) {
  *f = p[0] * (exp(-2 * p[1] * (r - p[2])) - 2. * exp(-p[1] * (r - p[2])));

  return;
}
