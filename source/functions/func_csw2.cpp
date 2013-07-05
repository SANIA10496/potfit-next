/****************************************************************
 *
 * func_csw2.cpp:
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

#include "func_csw2.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncCSW2::FuncCSW2() : power(0) {}

FuncCSW2::~FuncCSW2() {}

int FuncCSW2::num_params(void) {
  return 4;
}

void FuncCSW2::calc(const double &r, double *p, double *f) {
  power_1(power, r, p[3]);

  *f = (1. + p[0] * cos(p[1] * r + p[2])) / power;

  return;
}
