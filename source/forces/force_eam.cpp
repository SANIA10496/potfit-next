/****************************************************************
 *
 * force_eam.cpp:
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

#include "force_eam.h"
#include "../structures.h"

using namespace POTFIT_NS;

ForceEAM::ForceEAM(POTFIT *ptf): Force(ptf) {
  return;
}

ForceEAM::~ForceEAM() {
  return;
}

int ForceEAM::num_slots(void) {
  return 2;
}

int ForceEAM::neigh_type(void) {
  return 2;
}

int ForceEAM::get_col(int col, int a, int b) {
  return 0;
}

int ForceEAM::cols() {
}

void ForceEAM::read_additional_data(FILE *infile) {
}

double ForceEAM::calc_forces(double *forces) {
  return 0.0;
}
