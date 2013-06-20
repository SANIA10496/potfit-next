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
#include "interaction.h"
#include "io.h"
#include "settings.h"
#include "structures.h"
#include "utils.h"

using namespace POTFIT_NS;

Error::Error(POTFIT *ptf) :
  Pointers(ptf),
  force_vect(NULL),
  total_sum(1.0),
  force_sum(0.0),
  energy_sum(0.0),
  stress_sum(0.0),
  punish_sum(0.0),
  rms_force(0.0),
  rms_energy(0.0),
  rms_stress(0.0),
  total_contrib(0),
  num_forces(0),
  num_energies(0),
  num_stresses(0),
  fcalls(0)
{}

Error::~Error() {}

void Error::write_report(void) {
  calc_errors();

  io->write << std::endl << "###### error report ######" << std::endl;
  io->write << "total error sum " << std::fixed << std::setprecision(6) << total_sum << std::endl;
  io->write << "  " << total_contrib << " contributions ";
  io->write << "(" << num_forces << " forces, " << num_energies << " energies, " << num_stresses << " stresses)" << std::endl;

  io->write << "sum of force-errors    = " << std::setw(12) << std::setprecision(6);
  io->write << std::scientific << force_sum << "\t(";
  io->write << std::fixed << std::setw(6) << std::setprecision(2) << force_sum / total_sum * 100.;
  io->write << " % )" << std::endl;
  io->write << "sum of energy-errors   = " << std::setw(12) << std::setprecision(6);
  io->write << std::scientific << energy_sum << "\t(";
  io->write << std::fixed << std::setw(6) << std::setprecision(2) << energy_sum / total_sum * 100.;
  io->write << " % )" << std::endl;
  io->write << "sum of stress-errors   = " << std::setw(12) << std::setprecision(6);
  io->write << std::scientific << stress_sum << "\t(";
  io->write << std::fixed << std::setw(6) << std::setprecision(2) << stress_sum / total_sum * 100.;
  io->write << " % )" << std::endl;
  io->write << "sum of punishments     = " << std::setw(12) << std::setprecision(6);
  io->write << std::scientific << punish_sum << "\t(";
  io->write << std::fixed << std::setw(6) << std::setprecision(2) << punish_sum / total_sum * 100.;
  io->write << " % )" << std::endl;

  if (0 != punish_sum && 1 == settings->get_opt() )
    io->warning << "This sum contains punishments! Check your results." << std::endl;

  io->write << std::endl << "rms-errors:" << std::endl;
  io->write << "force\t" << std::setw(12) << std::setprecision(6) << std::scientific << rms_force << "\t(";
  io->write << std::fixed << std::setw(8) << std::setprecision(2) << rms_force * 1000. << " meV/\u212B)" << std::endl;
  io->write << "energy\t" << std::setw(12) << std::setprecision(6) <<  std::scientific << rms_energy << "\t(";
  io->write << std::fixed << std::setw(8) << std::setprecision(2) << rms_energy * 1000. << " meV/atom)" << std::endl;
  io->write << "stress\t" << std::setw(12) << std::setprecision(6) <<  std::scientific << rms_stress << "\t(";
  io->write << std::fixed << std::setw(8) << std::setprecision(2) << rms_stress / 160.2 * 1000. << " MPa)" << std::endl;

  if (1 == settings->get_opt()) {
    int time = utils->timediff();
    io->write << std::endl << "Runtime: " << time / 3600 << " hours, ";
    io->write << (time % 3600) / 60 << " minutes and ";
    io->write << time % 60 << " seconds." << std::endl;
    io->write << fcalls << " force calculations, each took ";
    io->write << std::scientific << static_cast<double>(time) / fcalls << " seconds." << std::endl;
  }

  return;
}

void Error::calc_errors(void) {
  int i, j, count = 0;

  total_sum = interaction->force->get_error_sum();

  num_forces = 3 * structures->get_num_contrib_atoms();
  num_energies = structures->get_num_contrib_energies();
  num_stresses = structures->get_num_contrib_stresses();
  total_contrib = num_forces + num_energies + num_stresses;

  // calculate errors
  for (i=0; i<structures->get_num_total_configs(); i++) {
    for (j=0; j<3 * structures->config[i]->num_atoms; j++) {
      force_sum += structures->config[i]->conf_weight * square(interaction->force->force_vect[count++]);
    }
    energy_sum += structures->config[i]->conf_weight * settings->get_eweight() *
                  square(interaction->force->force_vect[interaction->force->energy_p + i]);
    if (1 == structures->config[i]->use_stresses) {
      for (j=0; j<6; j++) {
        stress_sum += structures->config[i]->conf_weight * settings->get_sweight() *
                      square(interaction->force->force_vect[interaction->force->stress_p + 6 * i + j]);
      }
    }
  }

  // calculate rms errors
  count = 0;
  for (i=0; i<structures->get_num_total_configs(); i++) {
    for (j=0; j<3*structures->config[i]->num_atoms; j++) {
      rms_force += square(interaction->force->force_vect[count++]);
    }
    rms_energy += square(interaction->force->force_vect[interaction->force->energy_p + i]);
    if (1 == structures->config[i]->use_stresses) {
      for (j=0; j<6; j++) {
        rms_stress += square(interaction->force->force_vect[interaction->force->stress_p + 6 * i + j]);
      }
    }
  }
  rms_force = sqrt(rms_force / num_forces);
  rms_energy = sqrt(rms_energy / num_energies);
  rms_stress = sqrt(rms_stress / num_stresses);

  fcalls = interaction->force->get_fcalls();

  return;
}
