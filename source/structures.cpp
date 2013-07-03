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
#include <iomanip>

#include "communication.h"
#include "interaction.h"
#include "io.h"
#include "potential.h"
#include "structures.h"

using namespace POTFIT_NS;

Structures::Structures(POTFIT *ptf) :
  Pointers(ptf),
  using_forces(0),
  using_stresses(0),
  min_dist(NULL),
  ntypes(0),
  line(0),
  total_num_conf(0),
  total_num_atoms(0)
{}

Structures::~Structures() {
  if (NULL != min_dist)
    delete [] min_dist;

  for (unsigned i = 0; i < config.size(); ++i)
    delete config[i];
  config.clear();
  num_per_type.clear();

  return;
}

void Structures::init(void) {
  if (ntypes == 0) {
    io->error << "ntypes is 0!" << std::endl;
    io->pexit(EXIT_FAILURE);
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
    for (int i=0; i<ntypes; i++)
      num_per_type[i] += config[total_num_conf]->num_per_type[i];
    total_num_conf++;
  } while (!feof(infile));

  update_pointers(total_num_conf-1);

  communication->set_config_per_cpu();
  interaction->force->calc_pointers();

  return;
}

void Structures::update_pointers(const int &index) {
  int count = 0;
  int count_neigh = 0;
  int neigh_type = interaction->force->neigh_type();

  for (int i=0;i<=index;i++) {
    config[i]->atoms = &atoms[count];
    count += config[i]->num_atoms;
    if (2 == neigh_type) {
      for (int j=0;j<config[i]->num_atoms;j++) {
        config[i]->atoms[j].neighs = &neigh_2[count_neigh];
	count_neigh += config[i]->atoms[j].num_neighbors;
      }
    } else {
      for (int j=0;j<config[i]->num_atoms;j++) {
        config[i]->atoms[j].neighs = &neigh_3[count_neigh];
	count_neigh += config[i]->atoms[j].num_neighbors;
      }
    }
  }

  return;
}

void Structures::print_mindist(void) {
  int i, j, k;

  io->write << "Minimal Distances Matrix:" << std::endl;
  io->write << "Atom\t";
  for (i = 0; i < ntypes; i++)
    io->write << std::setw(8) << potential->elements[i] << "\t";
  io->write << std::endl;
  for (i = 0; i < ntypes; i++) {
    io->write << potential->elements[i] << "\t";
    for (j = 0; j < ntypes; j++) {
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2) : j * ntypes + i - ((j * (j + 1)) / 2);
      io->write << std::fixed << std::setw(8) << std::setprecision(4) << min_dist[k] << "\t";

    }
    io->write << std::endl;
  }
  io->write << std::endl;

  interaction->force->update_min_dist(min_dist);

  return;
}

const int Structures::get_num_contrib_atoms(void) {
  int count = 0;

  for (int i=0; i<total_num_conf; i++) {
    if (1 == config[i]->use_forces) {
      for (int j=0; j<config[i]->num_atoms; j++) {
        if (1 == config[i]->atoms[j].contrib)
          count++;
      }
    }
  }

  return count;
}

const int Structures::get_num_contrib_energies(void) {
  return total_num_conf;
}

const int Structures::get_num_contrib_stresses(void) {
  int count = 0;

  for (int i=0; i<total_num_conf; i++) {
    if (1 == config[i]->use_stresses)
      count++;
  }

  return 6 * count;
}

const int Structures::get_num_total_atoms(void) {
  return total_num_atoms;
}

const int Structures::get_num_total_configs(void) {
  return total_num_conf;
}

void Structures::set_ntypes(const int &n) {
  if (n < 1) {
    io->error << "ntypes cannot be smaller than 1!" << std::endl;
    io->pexit(EXIT_FAILURE);
  } else {
    ntypes = n;
  }

  return;
}

int Structures::get_ntypes(void) {
  return ntypes;
}
