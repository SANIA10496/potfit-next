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

#include <string>
#include <fstream>

#include "io.h"
#include "output.h"
#include "potential.h"
#include "settings.h"

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
  std::string component[3];
  std::string filename(output_prefix);
  component[0] = "x";
  component[1] = "y";
  component[2] = "z";
  filename.append(".force");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

  output.close();
  io->write << "Force data written to\t\t\t\t" << filename << std::endl;

  return;
}

void Output::write_output_file_energies(void) {
  std::string filename(output_prefix);
  filename.append(".energy");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

  output.close();

  io->write << "Energy data written to\t\t\t\t" << filename << std::endl;
  return;
}

void Output::write_output_file_stresses(void) {
  std::string component[6];
  std::string filename(output_prefix);
  component[0] = "x";
  component[1] = "y";
  component[2] = "z";
  filename.append(".stress");

  std::ofstream output;
  output.open (filename.c_str(), std::ofstream::out);

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
