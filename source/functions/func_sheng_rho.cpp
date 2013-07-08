/****************************************************************
 *
 * func_sheng_rho.cpp:
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

#include "func_sheng_rho.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncSHENG_RHO::FuncSHENG_RHO() : power(0.0), r6(0.0), r12(0.0) {}

FuncSHENG_RHO::~FuncSHENG_RHO() {}

int FuncSHENG_RHO::num_params(void) {
  return 5;
}

void FuncSHENG_RHO::calc(const double &r, const std::vector<double> &p, double *f) {
  if (r > 1.45) {
    power_1(power, r, p[1]);

    *f = p[0] * power + p[2];
  } else {
    r6 = (p[4] * p[4]) / (r * r);
    r6 = r6 * r6 * r6;
    r12 = r6 * r6;

    *f = 4.0 * p[3] * (r12 - r6);
  }

  return;
}
