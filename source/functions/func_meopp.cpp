/****************************************************************
 *
 * func_meopp.cpp:
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

#include "func_meopp.h"

#include "../utils.h"

using namespace POTFIT_NS;

FuncMEOPP::FuncMEOPP() {
  x[0] = 0.0;
  x[1] = 0.0;
  y[0] = 0.0;
  y[1] = 0.0;
  power[0] = 0.0;
  power[1] = 0.0;

  return;
}

FuncMEOPP::~FuncMEOPP() {}

int FuncMEOPP::num_params(void) {
  return 7;
}

void FuncMEOPP::calc(const double &r, const std::vector<double> &p, double *f) {
  x[0] = r - p[6];
  x[1] = r;
  y[0] = p[1];
  y[1] = p[3];

  power_m(2, power, x, y);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);

  return;
}
