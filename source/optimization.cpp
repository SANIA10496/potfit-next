/****************************************************************
 *
 * optimization.cpp:
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

#include <cstdlib>

#include "interaction.h"
#include "io.h"
#include "optimization.h"
#include "potential.h"
#include "settings.h"

#include "opt/list_opt.h"

using namespace POTFIT_NS;

Optimization::Optimization(POTFIT *ptf) : Pointers(ptf) {
  num_algs = 0;
  opt = NULL;

  return;
}

Optimization::~Optimization() {
  return;
}

void Optimization::run(void) {

  if (1 == settings->get_opt()) {
  io->write << "Starting optimization with the following parameters:" << std::endl;
  if (settings->get_eweight() != 0.0)
    io->write << " - Global energy weight " << settings->get_eweight() << std::endl;
  if (settings->get_sweight() != 0.0)
    io->write << " - Global stress weight "<< settings->get_sweight() << std::endl;
  io->write << " - " << potential->num_free_params << " free parameters" << std::endl << std::endl;

  if (0 == settings->get_myid()) {
    for (int i=0; i<num_algs; i++) {
      io->write << "Starting " << i + 1 << ". optimization" << std::endl;
      opt = init_algorithm(algorithms[i].c_str());
      if (NULL == opt) {
        io->error << "Could not find algorithm \"" << algorithms[i] << "\"." << std::endl;
	io->exit(EXIT_FAILURE);
      }
      opt->init(params[i]);
      io->write << " type: " << algorithms[i] << "\t parameters: " << params[i][0];
      io->write << ", " << params[i][1] << ", " << params[i][2] << std::endl;
      opt->run();
      io->write << "Finished " << i + 1 << ". optimization" << std::endl << std::endl;
      delete opt;
    }
  } else {
    interaction->calc_forces();
  }
  } else {
    io->write << "Optimization disabled." << std::endl << std::endl;
  }

  return;
}

void Optimization::add_algorithm(std::ifstream &infile) {
  std::string temp;

  // read algorithm name
  infile >> temp;
  algorithms.push_back(temp);
  params.push_back(new double[3]);

  // read params
  for (int i=0; i<3; i++) {
    infile >> temp;
    params[num_algs][i] = atof(temp.c_str());
  }
  num_algs++;

  return;
}

BaseOpt *Optimization::init_algorithm(const char *algo_type) {
  if (strcmp(algo_type,"none") == 0)
    return NULL;
#define OPTALG_TYPE
#define OptAlgType(key,Class) \
  else if (strcmp(algo_type,#key) == 0) \
    return new Class(ptf);
#include "opt/list_opt.h"
#undef OPTALG_TYPE
  return NULL;
}
