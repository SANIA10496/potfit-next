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

  if (1 == settings->opt) {
  io->write("Starting optimization with the following parameters:\n");
  if (settings->eweight != 0.0)
    io->write(" - Global energy weight %f\n",settings->eweight);
  if (settings->sweight != 0.0)
    io->write(" - Global stress weight %f\n",settings->sweight);
  io->write(" - %d free parameters\n\n",potential->num_free_params);

  if (0 == settings->myid) {
    for (int i=0; i<num_algs; i++) {
      io->write("Starting %d. optimization\n",i+1);
      opt = init_algorithm(algorithms[i].c_str());
      if (NULL == opt)
        io->error("Could not find algorithm \"%s\"",algorithms[i].c_str());
      opt->init(params[i]);
      io->write(" type: %s\t parameters: %f %f %f\n\n",algorithms[i].c_str(),params[i][0],params[i][1],params[i][2]);
      opt->run();
      io->write("Finished %d. optimization\n\n",i+1);
      delete opt;
    }
  } else {
    interaction->calc_forces();
  }
  } else {
    io->write("Optimization disabled.\n\n");
  }

  return;
}

void Optimization::add_algorithm(void) {
  char *str;

  // read algorithm name
  str = strtok(NULL, " \t\r\n");
  if (str == NULL)
    io->error("Algorithm name is missing!");
  algorithms.push_back(str);
  params.push_back(new double[3]);

  // read params
  for (int i=0; i<3; i++) {
    str = strtok(NULL, " \t\r\n");
    if (str != NULL)
      params[num_algs][i] = atof(str);
    else
      params[num_algs][i] = 0;
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
