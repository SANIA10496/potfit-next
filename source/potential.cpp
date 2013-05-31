/****************************************************************
 *
 * potential.cpp:
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
#include "interaction.h"
#include "io.h"
#include "memory.h"
#include "potential.h"
#include "structures.h"
#include "templates.h"

#include "tables/table.h"
#include "tables/list_tables.h"

using namespace POTFIT_NS;

Potential::Potential(POTFIT *ptf) : Pointers(ptf) {
  enable_cp = 0;
  num_globals = 0;

  num_pots = 0;
  num_free_pots = 0;

  num_params = 0;
  num_free_params = 0;

  invar_pot = NULL;

  idxpot = NULL;
  idxparam = NULL;

  rcut = NULL;
  rmin = NULL;
  rcut_max = 0.0;

  pots = NULL;
  opt = NULL;
  global_params = NULL;
  chem_pot = NULL;

  return;
}

Potential::~Potential() {
  delete [] invar_pot;
  delete [] rcut;
  delete [] rmin;

  delete [] idxpot;
  delete [] idxparam;

  for (unsigned i = 0; i < elements.size(); ++i)
    delete [] elements[i];
  elements.clear();

  for (int i=0; i<num_pots; i++)
    delete pots[i];
  delete [] pots;
  delete opt;

  if (global_params)
    delete global_params;

  if (chem_pot)
    delete chem_pot;

  return;
}

void Potential::init(int size) {
  num_pots = size;
  num_free_pots = size;
  invar_pot = new int[interaction->force->cols()];
  for (int i=0; i<interaction->force->cols(); i++)
    invar_pot[i] = 0;
  for (int i = 0; i < structures->get_ntypes(); ++i) {
    elements.push_back(new char[11]);
    sprintf(elements[i],"%d",i);
  }

  return;
}

void Potential::read_globals(FILE *infile) {
  char  buffer[1024], name[255];
  double val, min, max;
  int ret_val;
  fpos_t filepos;

  do {
    fgetpos(infile, &filepos);
    ret_val = fscanf(infile, "%s", buffer);
  } while (strcmp(buffer, "globals") != 0 && !feof(infile));
  fsetpos(infile, &filepos);

  // check for global keyword
  if (strcmp(buffer, "globals") == 0) {
    if (2 > fscanf(infile, "%s %d", buffer, &num_globals)) {
      io->error << "Global parameters are missing in the potential file." << std::endl;
      io->exit(EXIT_FAILURE);
    }

    global_params = new GlobalsTable(ptf, num_globals);

    // read the global parameters
    for (int j = 0; j < num_globals; j++) {
      ret_val = fscanf(infile, "%s %lf %lf %lf", name, &val, &min, &max);
      if (4 > ret_val)
        if (strcmp(name, "type") == 0) {
          io->error << "Not enough global parameters!" << std::endl;
          io->error << "You specified " << j << " parameter(s), but needed are " << num_globals << "." << std::endl;
	  io->exit(EXIT_FAILURE);
        }
      global_params->add_param(j, name, val, min, max);
    }
  } else {
    global_params = new GlobalsTable(ptf, 0);
  }

  io->write << "- Read " << num_globals << " global parameters" << std::endl;

  num_params += global_params->get_number_params();

  return;
}

void Potential::read_potentials(FILE *infile) {
  char  buffer[255], name[255];
  int   ret_val;
  fpos_t filepos;

  pots = new Table*[num_pots];

  for (int i=0; i<num_pots; i++) {
    do {
      fgetpos(infile, &filepos);
      ret_val = fscanf(infile, "%s", buffer);
    } while (strcmp(buffer, "type") != 0 && !feof(infile));
    fsetpos(infile, &filepos);
    // read type
    if (2 > fscanf(infile, "%s %s", buffer, name)) {
      io->error << "Premature end of potential file!" << std::endl;
      io->exit(EXIT_FAILURE);
    }
    if (strcmp(buffer, "type") != 0) {
      io->error << "Unknown keyword in potential file, expected \"type\" but found \"" << buffer << "\"." << std::endl;
      io->exit(EXIT_FAILURE);
    }
    if (strcmp(buffer, "table3") == 0) {
      pots[i] = new TableTab3(ptf);
    } else if (strcmp(buffer, "table4") == 0) {
      pots[i] = new TableTab4(ptf);
    } else
      pots[i] = new TableAnalytic(ptf);
    pots[i]->init(name, i);
    pots[i]->read_potential(infile);

    num_params += pots[i]->get_number_params();
    if (invar_pot[i] == 0)
      num_free_params += pots[i]->get_number_free_params();
    rcut_max = MAX(rcut_max, pots[i]->get_cutoff());
  }
  rcut = new double[square(structures->get_ntypes())];
  rmin = new double[square(structures->get_ntypes())];

  // TODO: this only works for pair interactions
  for (int i=0; i<structures->get_ntypes(); i++) {
    for (int j=0; j<structures->get_ntypes(); j++) {
      int k = (i <= j) ? i * structures->get_ntypes() + j - ((i * (i + 1)) / 2)
              : j * structures->get_ntypes() + i - ((j * (j + 1)) / 2);
      rcut[i*structures->get_ntypes()+j] = pots[k]->get_cutoff();
      rmin[i*structures->get_ntypes()+j] = pots[k]->get_rmin();
    }
  }

  io->write << "- Sucessfully read " << num_pots << " potentials" << std::endl;

  opt = new OptTable(ptf, num_free_params);
  init_opt_table();

  return;
}

void Potential::init_opt_table(void) {
  int count = 0;

  idxpot = new int[num_free_params];
  idxparam = new int[num_free_params];

  for (int i=0; i<num_pots;i++) {
    if (0 == invar_pot[i]) {
      for (int j=0; j<pots[i]->get_number_params(); j++) {
        if (pots[i]->invar_par[j] == 0) {
	  idxpot[count] = i;
	  idxparam[count] = j;
	  opt->values[count++] = pots[i]->values[j];
        }
      }
    }
  }

  if (count != num_free_params) {
    io->error << "Number of free parameters is inconsistent!" << std::endl;
    io->exit(EXIT_FAILURE);
  }

  return;
}

void Potential::update_potentials(int update) {
  for (int i=0;i<num_pots;i++) {
    if (0 == invar_pot[i] || update)
      pots[i]->update_calc_table(update);
  }

  return;
}

void Potential::set_enable_cp(int i) {
  enable_cp = i;

  return;
}

int Potential::get_enable_cp(void) {
  return enable_cp;
}
