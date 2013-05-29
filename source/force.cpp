/****************************************************************
 *
 * force.cpp:
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

#include "force.h"
#include "structures.h"

using namespace POTFIT_NS;

Force::Force(POTFIT *ptf) : Pointers(ptf) {
  force_vect = NULL;
  energy_p = 0;
  stress_p = 0;
  dummy_p = 0;
  limit_p = 0;
  punish_par_p = 0;
  punish_pot_p = 0;

  fcalls = 0;

  return;
}

Force::~Force() {
  delete [] force_vect;

  return;
}

void Force::calc_pointers(void) {
  energy_p = 3 * structures->get_num_total_atoms();
  stress_p = energy_p + structures->get_num_total_configs();
  limit_p = stress_p + 6 * structures->get_num_total_configs();
  dummy_p = limit_p + structures->get_num_total_configs();
  //  TODO
//  punish_par_p = dummy_p + 2 * ntypes;
//  punish_pot_p = punish_par_p + apot_table.total_par - apot_table.invar_pots;
//  punish_par_p = stress_p + 6 * nconf;
//  punish_pot_p = punish_par_p + apot_table.total_par;

  force_vect = new double[dummy_p];

  return;
}
