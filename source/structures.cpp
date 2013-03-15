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

#include "io.h"
#include "structures.h"

using namespace POTFIT_NS;

Structures::Structures(POTFIT *ptf) : Pointers(ptf) {
  ntypes = 0;
  num_atoms = 0;
  num_conf = 0;

  using_forces = 0;
  using_stresses = 0;

  elements = NULL;
  min_dist = NULL;

  line = 0;

  return;
}

Structures::~Structures() {
  for (int i=0; i<ntypes; i++)
    delete [] elements[i];
  delete [] elements;
  delete [] min_dist;

  for (unsigned i = 0; i < config.size(); ++i)
    delete config[i];
  config.clear();
  num_per_type.clear();

  return;
}

void Structures::init(void) {
  if (ntypes == 0) {
    io->error("Ntypes is 0!\n");
  }

  elements = new char*[ntypes];
  for (int i=0; i<ntypes; i++) {
    num_per_type.push_back(0);
    elements[i] = new char[3];
    strcpy(elements[i],"\0");
  }
  min_dist = new double[ntypes*ntypes];

  return;
}

void Structures::read_config(FILE *infile) {
  do {
    config.push_back(new Config(ptf));
    config[num_conf]->read(infile, &line);
    num_atoms += config[num_conf]->num_atoms;
    for (int i=0;i<ntypes;i++)
      num_per_type[i] += config[num_conf]->num_per_type[i];
    num_conf++;
  } while (!feof(infile));

  return;
}
