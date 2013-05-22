/****************************************************************
 *
 * communication.cpp:
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

#include <mpi.h>

#include "communication.h"
#include "io.h"
#include "settings.h"
#include "structures.h"

using namespace POTFIT_NS;

Communication::Communication(POTFIT *ptf) : Pointers(ptf) {
  return;
}

Communication::~Communication() {
  return;
}

void Communication::init(void) {
  if (1 < settings->num_cpus)
    io->write("Starting up MPI with %d processes.\n",settings->num_cpus);

  return;
}

void Communication::broadcast_params(void) {
  return;
}

void Communication::set_config_per_cpu(void) {
  if (1 == settings->num_cpus) {
    structures->firstconf = 0;
    structures->nconf = structures->total_num_conf;
  } else {
    int each = (structures->total_num_conf / settings->num_cpus);
    int odd = (structures->total_num_conf % settings->num_cpus) - settings->num_cpus;

    structures->firstconf = settings->myid * each;
    structures->nconf = each;
    if (settings->myid < odd)
      structures->nconf++;
  }

  return;
}
