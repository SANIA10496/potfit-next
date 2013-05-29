/****************************************************************
 *
 * error.cpp:
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

#include <cmath>

#include "error.h"
#include "interaction.h"
#include "io.h"
#include "settings.h"
#include "structures.h"
#include "utils.h"

using namespace POTFIT_NS;

Error::Error(POTFIT *ptf) : Pointers(ptf) {
  total_sum = 1.0;
  force_sum = 0.2;
  energy_sum = 0.0;
  stress_sum = 0.0;
  punish_sum = 0.0;
  rms_force = 0.0;
  rms_energy = 0.0;
  rms_stress = 0.0;
  total_contrib = 0;
  num_forces = 0;
  num_energies = 0;
  num_stresses = 0;
  fcalls = 1;

  force_vect = NULL;

  return;
}

Error::~Error() {
  return;
}

void Error::write_report(void) {
  calc_errors();

  io->write("\n###### error report ######\n");
  io->write("total error sum %f\n",total_sum);
  io->write("\t%d contributions ",total_contrib);
  io->write("(%d forces, %d energies, %d stresses)\n", num_forces, num_energies, num_stresses);

  io->write("sum of force-errors  = %e\t(%7.3f%% - avg: %f)\n", force_sum, force_sum / total_sum * 100.,
	force_sum / structures->get_num_total_atoms());
  io->write("sum of energy-errors = %e\t(%7.3f%%)\n", energy_sum, energy_sum / total_sum * 100.);
  io->write("sum of stress-errors = %e\t(%7.3f%%)\n", stress_sum, stress_sum / total_sum * 100.);
  io->write("sum of punishments   = %e\t(%7.3f%%)\n", punish_sum, punish_sum / total_sum * 100.);

  if (0 != punish_sum && 1 == settings->opt ) {
    io->warning("This sum contains punishments! Check your results.\n");
  } else {
    io->write("\n");
  }

  io->write("rms-errors:\n");
  io->write("force\t%e\t(%10.3f meV/A)\n", rms_force, rms_force * 1000.);
  io->write("energy\t%e\t(%10.3f meV)\n", rms_energy, rms_energy * 1000.);
  io->write("stress\t%e\t(%10.3f MPa)\n", rms_stress, rms_stress / 160.2 * 1000.);

  if (1 == settings->opt) {
    io->write("\nRuntime: %d hours, %d minutes and %d seconds.\n",
	utils->timediff() / 3600, (utils->timediff() % 3600) / 60, utils->timediff() % 60);
    io->write("%d force calculations, each took %f seconds.\n", fcalls, utils->timediff() / fcalls);
  }

  return;
}

void Error::calc_errors(void) {
  double temp = 0.0;
  int i, j, count = 0;

  // calculate forces with current potential
  total_sum = interaction->calc_forces();

  num_forces = 3 * structures->get_num_contrib_atoms();
  num_energies = structures->get_num_contrib_energies();
  num_stresses = structures->get_num_contrib_stresses();
  total_contrib = num_forces + num_energies + num_stresses;

  // calculate errors
  for (i=0;i<structures->get_num_total_configs();i++) {
    for (j=0;j<3 * structures->config[i]->num_atoms;j++) {
      force_sum += structures->config[i]->conf_weight * square(interaction->force->force_vect[count++]);
    }
    energy_sum += structures->config[i]->conf_weight * settings->eweight *
	    square(interaction->force->force_vect[interaction->force->energy_p + i]);
    if (1 == structures->config[i]->use_stresses) {
      for (j=0;j<6;j++) {
        stress_sum += structures->config[i]->conf_weight * settings->sweight *
		square(interaction->force->force_vect[interaction->force->stress_p + 6 * i + j]);
      }
    }
  }

  // calculate rms errors
  for (i=0;i<structures->get_num_total_configs();i++) {
    for (j=0;j<structures->config[i]->num_atoms;j++) {
      rms_force += square(interaction->force->force_vect[count++]);
    }
    rms_energy += square(interaction->force->force_vect[interaction->force->energy_p + i]);
    if (1 == structures->config[i]->use_stresses) {
      for (j=0;j<6;j++) {
        rms_stress += square(interaction->force->force_vect[interaction->force->stress_p + 6 * i + j]);
      }
    }
  }
  rms_force = sqrt(rms_force / (3. * structures->get_num_contrib_atoms()));
  rms_energy = sqrt(rms_energy / structures->get_num_total_configs());
  rms_stress = sqrt(rms_stress / ( 6. * structures->get_num_total_configs()));

  return;
}
