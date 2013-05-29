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
#include "settings.h"

using namespace POTFIT_NS;

Output::Output(POTFIT *ptf) : Pointers(ptf) {
  imdpotsteps = 1000;
  enable_distfile = 0;
  enable_output_files = 0;
  enable_imd_pot = 0;
  enable_lammps_pot = 0;
  enable_plot_file = 0;
  enable_log = 0;

  plotmin = 0.1;

  return;
}

Output::~Output() {
  return;
}

void Output::write_output(void) {
  write_endpot();

  if (1 == enable_imd_pot)
    write_imdpot();

  if (1 == enable_plot_file)
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
  io->write("Final potential written to file\t\t\t%s\n",endpot.c_str());
  return;
}

void Output::write_imdpot(void) {
  io->write("IMD potential written to file\t\t\t%s\n",imdpot.c_str());
  return;
}

void Output::write_plotfile(void) {
  io->write("Potential plotting data written to file\t\t%s\n",plotfile.c_str());
  return;
}

void Output::write_lammpspot(void) {
  io->write("Potential in LAMMPS format written to\t\t%s\n",plotfile.c_str());
  return;
}

void Output::write_output_files(void) {
  // always write force output file
  write_output_file_forces();

  // only write energy output file if eweight != 0
  if (0 != settings->eweight)
    write_output_file_energies();

  // only write energy output file if sweight != 0
  if (0 != settings->sweight)
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
  output.open ("test.txt", std::ofstream::out | std::ofstream::app);


  output.close();
  io->write("Force data written to\t\t\t\t%s.force\n",output_prefix.c_str());

  return;
}

void Output::write_output_file_energies(void) {
  io->write("Energy data written to\t\t\t\t%s.energy\n",output_prefix.c_str());
  return;
}

void Output::write_output_file_stresses(void) {
  io->write("Stress data written to\t\t\t\t%s.stress\n",output_prefix.c_str());
  return;
}

void Output::write_output_file_punishments(void) {
  io->write("Punishments written to\t\t\t\t%s.punish\n",output_prefix.c_str());
  return;
}
