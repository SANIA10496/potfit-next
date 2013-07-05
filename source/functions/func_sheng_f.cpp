/****************************************************************
 *
 * func_sheng_f.cpp:
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

#include "func_sheng_f.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncSHENG_F::FuncSHENG_F() : power(0.0) {}

FuncSHENG_F::~FuncSHENG_F() {}

int FuncSHENG_F::num_params(void) {
  return 4;
}

void FuncSHENG_F::calc(const double &r, double *p, double *f) {
  power_1(power, r, p[1]);

  *f = p[0] * power + p[2] * r + p[3];

  return;
}
