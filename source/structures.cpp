/****************************************************************
 *
 * structures.cpp:
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

#include <cstring>
#include <cstdlib>

#include "communication.h"
#include "interaction.h"
#include "io.h"
#include "potential.h"
#include "structures.h"

using namespace POTFIT_NS;

Structures::Structures(POTFIT *ptf) : Pointers(ptf) {
  ntypes = 0;
  total_num_atoms = 0;
  total_num_conf = 0;

  using_forces = 0;
  using_stresses = 0;

  min_dist = NULL;

  line = 0;

  return;
}

Structures::~Structures() {
  delete [] min_dist;

  for (unsigned i = 0; i < config.size(); ++i)
    delete config[i];
  config.clear();
  num_per_type.clear();

  return;
}

void Structures::init(void) {
  if (ntypes == 0) {
    io->error("ntypes is 0!\n");
  }

  min_dist = new double[ntypes*ntypes];

  for (int i=0; i<ntypes; i++) {
    num_per_type.push_back(0);
    for (int j=0; j<ntypes; j++)
      min_dist[i*ntypes+j] = 999.9;
  }

  return;
}

void Structures::read_config(FILE *infile) {
  do {
    config.push_back(new Config(ptf, total_num_conf));
    config[total_num_conf]->read(infile, &line);
    total_num_atoms += config[total_num_conf]->num_atoms;
    for (int i=0;i<ntypes;i++)
      num_per_type[i] += config[total_num_conf]->num_per_type[i];
    total_num_conf++;
  } while (!feof(infile));

  communication->set_config_per_cpu();
  interaction->force->calc_pointers();

  return;
}

void Structures::print_mindist(void) {
  int i, j, k;

  io->write("Minimal Distances Matrix:\n");
  io->write("Atom\t");
  for (i = 0; i < ntypes; i++)
    io->write("%8s\t", potential->elements[i]);
  io->write("with\n");
  for (i = 0; i < ntypes; i++) {
    io->write("%s\t", potential->elements[i]);
    for (j = 0; j < ntypes; j++) {
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2) : j * ntypes + i - ((j * (j + 1)) / 2);
      io->write("%f\t", min_dist[k]);

    }
    io->write("\n");
  }
  io->write("\n");

  return;
}
