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

#include "communication.h"
#include "interaction.h"
#include "io.h"
#include "potential.h"
#include "settings.h"
#include "structures.h"

using namespace POTFIT_NS;

Communication::Communication(POTFIT *ptf) :
  Pointers(ptf)
{
  MPI_VECTOR = MPI::DOUBLE.Create_contiguous(3);
  MPI_VECTOR.Commit();
  MPI_STENS = MPI::DOUBLE.Create_contiguous(6);
  MPI_STENS.Commit();

  return;
}

Communication::~Communication() {
  MPI_VECTOR.Free();
  MPI_STENS.Free();

  return;
}

void Communication::init(void) {
  if (1 < settings->get_num_cpus())
    io->write << "Starting up MPI with " << settings->get_num_cpus() << " processes." << std::endl;

  return;
}

void Communication::broadcast_params(void) {
  MPI::COMM_WORLD.Bcast(potential->opt->val_p,opt_pot_len,MPI_DOUBLE,0);

  return;
}

double Communication::gather_forces(const double &tmpsum) {
  double sum(0.0);

  if (settings->get_num_cpus() > 1) {
    MPI::COMM_WORLD.Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0);
    // gather forces, energies, stresses
    if (settings->get_myid() == 0) {		// root node already has data in place
      // forces
      MPI::COMM_WORLD.Gatherv(MPI_IN_PLACE, structures->myatoms, MPI_VECTOR,
                              interaction->force->force_vect, structures->atom_len, structures->atom_dist, MPI_VECTOR, 0);
      // energies
      MPI::COMM_WORLD.Gatherv(MPI_IN_PLACE, structures->nconf, MPI_DOUBLE,
                              interaction->force->force_vect + interaction->force->energy_p,
                              structures->conf_len, structures->conf_dist, MPI_DOUBLE, 0);
      // stresses
      MPI::COMM_WORLD.Gatherv(MPI_IN_PLACE, structures->nconf, MPI_STENS,
                              interaction->force->force_vect + interaction->force->stress_p,
                              structures->conf_len, structures->conf_dist, MPI_STENS, 0);
    } else {
      // forces
      MPI::COMM_WORLD.Gatherv(interaction->force->force_vect + structures->firstatom, structures->myatoms, MPI_VECTOR,
                              interaction->force->force_vect, structures->atom_len, structures->atom_dist, MPI_VECTOR, 0);
      // energies
      MPI::COMM_WORLD.Gatherv(interaction->force->force_vect + interaction->force->energy_p + structures->firstconf, structures->nconf,
                              MPI_DOUBLE, interaction->force->force_vect + interaction->force->energy_p,
                              structures->conf_len, structures->conf_dist, MPI_DOUBLE, 0);
      // stresses
      MPI::COMM_WORLD.Gatherv(interaction->force->force_vect + interaction->force->stress_p + structures->firstconf,
                              structures->nconf, MPI_STENS, interaction->force->force_vect + interaction->force->stress_p,
                              structures->conf_len, structures->conf_dist, MPI_STENS, 0);
    }
  } else {
    sum = tmpsum;
  }

  return sum;
}

void Communication::set_config_per_cpu(void) {
  const int id = settings->get_myid();
  const int cpus = settings->get_num_cpus();
  const int nconf = structures->get_num_total_configs();
  const int each = (nconf / cpus);
  const int odd = (nconf % cpus) - cpus;
  int i(0);

  if (1 == cpus) {
    structures->firstconf = 0;
    structures->nconf = nconf;
  } else {
    structures->firstconf = id * each + ((id + odd) > 0 ? (id + odd) : 0);
    structures->nconf = each + ((id + odd) >= 0 ? 1 : 0);
  }

  opt_pot_len = potential->get_num_free_params();

  for (i=0; i<structures->nconf; i++) {
    structures->myatoms += structures->config[structures->firstconf + i]->num_atoms;
  }

  structures->atom_len = new int[cpus];
  structures->atom_dist = new int[cpus];
  structures->conf_len = new int[cpus];
  structures->conf_dist = new int[cpus];

  for (i = 0; i < cpus; i++)
    structures->conf_dist[i] = i * each + ((i + odd) > 0 ? (i + odd) : 0);
  for (i = 0; i < cpus - 1; i++)
    structures->conf_len[i] = structures->conf_dist[i + 1] - structures->conf_dist[i];
  structures->conf_len[cpus - 1] = structures->get_num_total_configs() - structures->conf_dist[cpus - 1];
  for (i = 0; i < cpus; i++)
    structures->atom_dist[i] = structures->config[structures->conf_dist[i]]->cnfstart;
  for (i = 0; i < cpus - 1; i++)
    structures->atom_len[i] = structures->atom_dist[i + 1] - structures->atom_dist[i];
  structures->atom_len[cpus - 1] = structures->get_num_total_atoms() - structures->atom_dist[cpus - 1];

  structures->firstatom = structures->atom_dist[id];

  return;
}
