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

#include<cstdlib>

#include "io.h"
#include "optimization.h"
#include "potential.h"
#include "settings.h"

using namespace POTFIT_NS;

Optimization::Optimization(POTFIT *ptf) : Pointers(ptf) {}

Optimization::~Optimization() {}

void Optimization::run(void) {

  io->write("Starting optimization with the following parameters:\n");
  if (settings->eweight != 0.0)
    io->write(" - Global energy weight %f\n",settings->eweight);
  if (settings->sweight != 0.0)
    io->write(" - Global stress weight %f\n",settings->sweight);
  io->write(" - %d free parameters\n\n",potential->num_free_params);

  io->write("num_algs = %d\n",num_algs);
  for (int i=0;i<num_algs;i++)
	  io->write("alg_%d: %s\n",i,algorithms[i]);
  return;
}

void Optimization::add_algorithm(void) {
  char *str;

  // read algorithm name
  str = strtok(NULL, " \t\r\n");
  if (str == NULL)
    io->error("Algorithm name is missing!");
  algorithms.push_back(new char[20]);
  sprintf(algorithms[num_algs++], str);

  // read maxsteps
  str = strtok(NULL, " \t\r\n");
  if (str == NULL)
    io->error("Algorithm name is missing!");
  maxsteps.push_back(atoi(str));

  return;
}
