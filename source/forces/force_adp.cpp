/****************************************************************
 *
 * force_adp.cpp:
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

#include "force_adp.h"
#include "../structures.h"

using namespace POTFIT_NS;

ForceADP::ForceADP(POTFIT *ptf): Force(ptf) {
  return;
}

ForceADP::~ForceADP() {
  return;
}

void ForceADP::read_additional_data(FILE *infile) {
  return;
};

int ForceADP::num_slots(void) {
  return 4;
}

int ForceADP::neigh_type(void) {
  return 2;
}

int ForceADP::get_col(const int &slot, const int &a, const int &b) {
  return 0;
}

void ForceADP::update_min_dist(double *min_dist) {
  return;
}

int ForceADP::cols(void) {
  int n = structures->get_ntypes();
  return (int)n*(n+1)/2.;
}

double ForceADP::calc_forces(void) {
  return 0.0;
}

void ForceADP::write_imd_pot(void) {

  return;
}

void ForceADP::write_lammps_pot(void)  {

  return;
}
