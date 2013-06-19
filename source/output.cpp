/****************************************************************
 *
 * output.cpp:
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

#include <fstream>
#include <string>

#include "interaction.h"
#include "io.h"
#include "output.h"
#include "potential.h"
#include "settings.h"
#include "structures.h"

using namespace POTFIT_NS;

Output::Output(POTFIT *ptf) :
  Pointers(ptf),
  plotmin(0.1),
  imdpotsteps(1000),
  enable_distfile(0),
  enable_output_files(0),
  enable_imd_pot(0),
  enable_lammps_pot(0),
  enable_plotfile(0),
  enable_log(0)
{}


Output::~Output() {}

void Output::write_output(void) {
  write_endpot();

  if (1 == enable_imd_pot)
    write_imdpot();

  if (1 == enable_plotfile)
    write_plotfile();

  if (1 == enable_lammps_pot)
    write_lammpspot();

  if (1 == enable_output_files)
    write_output_files();

  return;
}

void Output::write_tempfile(void) {
  return;
}

void Output::write_endpot(void) {
  std::ofstream outfile;

  outfile.exceptions(std::ofstream::failbit);

  try {
    outfile.open(endpot.c_str(), std::ofstream::out);
  } catch (std::ofstream::failure e) {
    io->error << "Could not open " << endpot << " file for writing." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  potential->write_potentials(outfile);

  outfile.close();

  io->write << "Final potential written to file\t\t\t" << endpot << std::endl;

  return;
}

void Output::write_imdpot(void) {
  io->write << "IMD potential written to file\t\t\t" << imdpot << std::endl;
  return;
}

void Output::write_plotfile(void) {
  io->write << "Potential plotting data written to file\t\t" << plotfile << std::endl;
  return;
}

void Output::write_lammpspot(void) {
  io->write << "Potential in LAMMPS format written to\t\t" << plotfile << std::endl;
  return;
}

void Output::write_output_files(void) {
  // always write force output file
  write_output_file_forces();

  // only write energy output file if eweight != 0
  if (0 != settings->get_eweight())
    write_output_file_energies();

  // only write energy output file if sweight != 0
  if (0 != settings->get_sweight())
    write_output_file_stresses();

  // always write punishments file
  write_output_file_punishments();
  return;
}

void Output::write_output_file_forces(void) {
  int i, j, count = 0;
  double sqr = 0.0;
  double *force = interaction->force->force_vect;
  std::string filename(output_prefix);
  Config *conf;

  filename.append(".force");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

  output << "#conf:atom\ttype\tdf^2\t\tf\t\tf0\t\tdf/f0" << std::endl;

  for (i=0; i<structures->get_num_total_configs(); i++) {
    conf = structures->config[i];
    for (j=0; j<conf->num_atoms; j++) {
      sqr = conf->conf_weight * square(force[count + j]);
      output << std::setw(4) << i << ":";
      output << std::setw(6) << j << ":x";
      output << "\t" << potential->elements[conf->atoms[j]->type];
      output << "\t" << std::scientific << sqr;
      output << "\t" << force[count + 3 * j] + conf->atoms[j]->force.x;
      output << "\t" << conf->atoms[j]->force.x;
      output << "\t" << force[count + 3 * j] / conf->atoms[j]->force.x;
      output << std::endl;
      sqr = conf->conf_weight * square(force[count + j + 1]);
      output << std::setw(4) << i << ":";
      output << std::setw(6) << j << ":y";
      output << "\t" << potential->elements[conf->atoms[j]->type];
      output << "\t" << std::scientific << sqr;
      output << "\t" << force[count + 3 * j + 1] + conf->atoms[j]->force.y;
      output << "\t" << conf->atoms[j]->force.y;
      output << "\t" << force[count + 3 * j + 1] / conf->atoms[j]->force.y;
      output << std::endl;
      sqr = conf->conf_weight * square(force[count + j + 2]);
      output << std::setw(4) << i << ":";
      output << std::setw(6) << j << ":z";
      output << "\t" << potential->elements[conf->atoms[j]->type];
      output << "\t" << std::scientific << sqr;
      output << "\t" << force[count + 3 * j + 2] + conf->atoms[j]->force.z;
      output << "\t" << conf->atoms[j]->force.z;
      output << "\t" << force[count + 3 * j + 2] / conf->atoms[j]->force.z;
      output << std::endl;
    }
    if (i<structures->get_num_total_configs()-1)
      output << std::endl;
    count += structures->config[i]->num_atoms;
  }

  output.close();
  io->write << "Force data written to\t\t\t\t" << filename << std::endl;

  return;
}

void Output::write_output_file_energies(void) {
  double sqr = 0.0;
  double *force = interaction->force->force_vect + interaction->force->energy_p;
  std::string filename(output_prefix);

  filename.append(".energy");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

  output << "# global energy weight w is " << settings->get_eweight() << std::endl;
  output << "# nr.\tconf_w\tw*de^2\t\te\t\te0\t\te-e0\t\tde/e0 [%]" << std::endl;

  for (int i=0; i<structures->get_num_total_configs(); i++) {
    sqr = settings->get_eweight() * structures->config[i]->conf_weight *
          square(force[i]);
    output << std::setw(4) << i << "\t";
    output << std::fixed << std::setprecision(2) << structures->config[i]->conf_weight << "\t";
    output << std::setw(8) << std::setprecision(6) << sqr << "\t";
    output << force[i] + structures->config[i]->coh_energy << "\t";
    output << structures->config[i]->coh_energy << "\t";
    output << force[i] << "\t";
    output << std::fabs(force[i] / structures->config[i]->coh_energy) * 100. << std::endl;
  }

  output.close();

  io->write << "Energy data written to\t\t\t\t" << filename << std::endl;
  return;
}

void Output::write_output_file_stresses(void) {
  double sqr = 0.0;
  double *force = interaction->force->force_vect + interaction->force->stress_p;
  std::string filename(output_prefix);
  const std::string component[] = {"xx", "yy", "zz", "xy", "yz", "zx"};

  filename.append(".stress");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

  output << "# global stress weight w is " << settings->get_sweight() << std::endl;
  output << "# nr.\tconf_w\tw*ds^2\t\ts\t\ts0\t\ts-s0\t\tds/s0 [%]" << std::endl;

  for (int i=0; i<structures->get_num_total_configs(); i++) {
    for (int j=0; j<6; j++) {
      sqr = settings->get_sweight() * structures->config[i]->conf_weight *
            square(force[6*i + j]);
      output << std::setw(4) << i << ":" << component[j] << "\t";
      output << std::fixed << std::setprecision(2) << structures->config[i]->conf_weight << "\t";
      output << std::setw(8) << std::setprecision(6) << sqr << "\t";
      output << force[6 * i + j] + *(structures->config[i]->dstress[j]) << "\t";
      output << *(structures->config[i]->dstress[j]) << "\t";
      output << force[6 * i + j] << "\t";
      output << std::fabs(force[6 * i + j] / *(structures->config[i]->dstress[j])) * 100. << std::endl;
    }
  }

  output.close();
  io->write << "Stress data written to\t\t\t\t" << filename << std::endl;

  return;
}

void Output::write_output_file_punishments(void) {
  std::string filename(output_prefix);
  filename.append(".punish");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

  output.close();

  io->write << "Punishments written to\t\t\t\t" << filename << std::endl;

  return;
}

void Output::set_endpot(std::string str) {
  endpot = str;

  return;
}

std::string Output::get_endpot(void) {
  return endpot;
}

void Output::set_output_prefix(std::string str) {
  output_prefix = str;
  if (!output_prefix.empty())
    enable_output_files = 1;

  return;
}

std::string Output::get_output_prefix(void) {
  return output_prefix;
}

void Output::set_imdpot(std::string str) {
  imdpot = str;
  enable_imd_pot = 1;

  return;
}

std::string Output::get_imdpot(void) {
  return imdpot;
}

void Output::set_plotfile(std::string str) {
  plotfile = str;
  enable_plotfile = 1;

  return;
}

std::string Output::get_plotfile(void) {
  return plotfile;
}

void Output::set_distfile(std::string str) {
  distfile = str;
  enable_distfile = 1;

  return;
}

std::string Output::get_distfile(void) {
  return distfile;
}

void Output::set_imdpotsteps(int n) {
  imdpotsteps = n;

  return;
}

int Output::get_imdpotsteps(void) {
  return imdpotsteps;
}

void Output::set_plotmin(double x) {
  plotmin = x;

  return;
}

double Output::get_plotmin(void) {
  return plotmin;
}

void Output::set_tempfile(std::string str) {
  tempfile = str;

  return;
}

std::string Output::get_tempfile(void) {
  return tempfile;
}

void Output::set_enable_pairdist(int i) {
  enable_pairdist = i;

  return;
}

int Output::get_enable_pairdist(void) {
  return enable_pairdist;
}

void Output::set_plotpointfile(std::string str) {
  plotpointfile = str;

  return;
}

std::string Output::get_plotpointfile(void) {
  return plotpointfile;
}

void Output::set_lammps_pot(int i) {
  enable_lammps_pot = i;

  return;
}
