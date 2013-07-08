/****************************************************************
 *
 * func_mishin.cpp:
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

#include "func_mishin.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncMISHIN::FuncMISHIN() : power(0.0), temp(0.0), z(0.0) {}

FuncMISHIN::~FuncMISHIN() {}

int FuncMISHIN::num_params(void) {
  return 6;
}

void FuncMISHIN::calc(const double &r, const std::vector<double> &p, double *f) {
  z = r - p[3];
  temp = exp(-p[5] * r);

  power_1(power, z, p[4]);

  *f = p[0] * power * temp * (1. + p[1] * temp) + p[2];

  return;
}
