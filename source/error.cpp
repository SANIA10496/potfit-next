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

#include "error.h"
#include "io.h"
#include "settings.h"
#include "utils.h"

using namespace POTFIT_NS;

Error::Error(POTFIT *ptf) : Pointers(ptf) {
  total_sum = 0.0;
  force_sum = 0.0;
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

  io->write("sum of force-errors  = %f\t\t (%f%% - av: %f)\n", force_sum, force_sum / total_sum * 100., force_sum);
  io->write("sum of energy-errors = %f\t\t (%f%%)\n", energy_sum, energy_sum / total_sum * 100.);
  io->write("sum of stress-errors = %f\t\t (%f%%)\n", stress_sum, stress_sum / total_sum * 100.);
  io->write("sum of punishments   = %f\t\t (%f%%)\n", punish_sum, punish_sum / total_sum * 100.);

  if (0 == punish_sum) {
    io->warning("This sum contains punishments! Check your results.\n");
  }

  io->write("rms-errors:\n");
  io->write("force\t%f\t(%f meV/A)\n", rms_force, rms_force * 1000.);
  io->write("energy\t%f\t(%f meV/A)\n", rms_energy, rms_energy * 1000.);
  io->write("stress\t%f\t(%f meV/A)\n", rms_stress, rms_stress / 160.2 * 1000.);

  if (1 == settings->opt) {
    io->write("\nRuntime: %d hours, %d minutes and %d seconds.\n",
	utils->timediff() / 3600, (utils->timediff() % 3600) / 60, utils->timediff() % 60);
    io->write("%d force calculations, each took %f seconds\n", fcalls, utils->timediff() / fcalls);
  }

  return;
}

void Error::calc_errors(void) {

  return;
}
