/****************************************************************
 *
 * func_eopp.cpp:
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

#include "func_eopp.h"

using namespace POTFIT_NS;

FuncEOPP::FuncEOPP() {
  power[0] = 0.0;
  power[1] = 0.0;

  return;
}

FuncEOPP::~FuncEOPP() {
  return;
}

int FuncEOPP::num_params(void) {
  return 6;
}

void FuncEOPP::calc(double r, double *p, double *f) {
  power[0] = pow(r, p[1]);
  power[1] = pow(r, p[3]);

  *f = p[0] / power[0] + (p[2] / power[1]) * cos(p[4] * r + p[5]);

  return;
}
