/****************************************************************
 *
 * func_lj.cpp:
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

#include "func_lj.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncLJ::FuncLJ() : r6(0.0), r12(0.0) {}

FuncLJ::~FuncLJ() {}

int FuncLJ::num_params(void) {
  return 2;
}

void FuncLJ::calc(const double &r, const std::vector<double> &p, double *f) {
  r6 = (p[1] * p[1]) / (r * r);
  r6 = r6 * r6 * r6;
  r12 = square(r6);

  *f = 4.0 * p[0] * (r12 - r6);

  return;
}
