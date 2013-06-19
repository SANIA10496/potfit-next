/****************************************************************
 *
 * force_pair.cpp:
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

#include "force_pair.h"
#include "../communication.h"
#include "../io.h"
#include "../potential.h"
#include "../settings.h"
#include "../tables/table_analytic.h"
#include "../templates.h"

using namespace POTFIT_NS;

ForcePair::ForcePair(POTFIT *ptf):
  Force(ptf),
  h(0),i(0),j(0),k(0),l(0),
  self(0), uf(0),
  us(0), stresses(0),
  tmpsum(0.0),sum(0.0),
  phi_val(0.0),phi_grad(0.0),
  atom(NULL),
  conf(NULL),
  neigh(NULL),
  pot(NULL)
{}

ForcePair::~ForcePair() {}

int ForcePair::num_slots(void) {
  return 1;
}

int ForcePair::neigh_type(void) {
  return 2;
}

int ForcePair::get_col(int slot, int a, int b) {
  int col;

  if (slot != 0) {
    io->error << "Pair potentials only use 1 slot." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  col = (a <= b) ? a * structures->get_ntypes() + b - ((a * (a + 1)) / 2)
        : b * structures->get_ntypes() + a - ((b * (b + 1)) / 2);

  return col;
}

int ForcePair::cols() {
  int n = structures->get_ntypes();

  return (int)n*(n+1)/2.;
}

void ForcePair::update_min_dist(double *min_dist) {
  int n = structures->get_ntypes();
  k = 0;

  for (int i=0;i<n;i++)
    for (int j=0;j<n;j++) {
      k = (i <= j) ? i * n + j - ((i * (i + 1)) / 2) : j * n + i - ((j * (j + 1)) / 2);
      potential->pots[i*n+j]->begin = 0.95 * min_dist[k];
    }

  potential->update_potentials(1);
  potential->update_slots();

  return;
}

/****************************************************************
 *
 * ForcePair::read_additional_data(FILE *infile)
 * 	read the chemical potentials (if present)
 *
 ****************************************************************/

void ForcePair::read_additional_data(FILE *infile) {
  char  buffer[255], name[255];
  char *token;
  double val, min, max;
  int ret_val;

  if (potential->get_enable_cp() == 1) {

    // search for chempot keyword
    do {
      fscanf(infile, "%s", buffer);
    } while (strcmp(buffer, "chempot") != 0 && !feof(infile));

    potential->chem_pot = new ChempotTable(ptf, structures->get_ntypes());

    // loop over all atom types
    for (int j = 0; j < structures->get_ntypes(); j++) {

      // read one line
      if (4 > fscanf(infile, "%s %lf %lf %lf", buffer, &val, &min, &max)) {
        io->error << "Could not read chemical potential for atomtype " << j + 1 << "." << std::endl;
        io->pexit(EXIT_FAILURE);
      }

      // split cp and _#
      token = strchr(buffer, '_');
      if (token != NULL) {
        strncpy(name, buffer, strlen(buffer) - strlen(token));
        name[strlen(buffer) - strlen(token)] = '\0';
      }
      if (strcmp("cp", name) != 0) {
        io->error << "Found \"" << name << "\" instead of \"cp\"" << std::endl;
        io->error << "No chemical potentials found in potential file." << std::endl;
	io->pexit(EXIT_FAILURE);
      }
      potential->chem_pot->add_value(j, buffer, val, min, max);
    }
    io->write << "- Enabled chemical potentials for " << structures->get_ntypes() << " elements." << std::endl;
  }

  return;
}

/****************************************************************
 *
 * ForcePair::calc_forces(double *forces)
 * 	read the chemical potentials (if present)
 *
 ****************************************************************/

double ForcePair::calc_forces(void) {

  // This is the start of an infinite loop
  while (1) {
    communication->broadcast_params();
    potential->update_potentials(0);

    tmpsum = 0.0;		// sum of squares of local process

    atom_count = 0;
    phi_val = 0.0;
    phi_grad = 0.0;
    conf = NULL;
    atom = NULL;
    neigh = NULL;
    pot = NULL;

    // loop over configurations
    for (h = structures->firstconf; h < structures->firstconf + structures->nconf; h++) {
      conf = structures->config[h];
      uf = conf->use_forces;
      us = conf->use_stresses;

      // reset energies and stresses
      force_vect[energy_p + h] = 0.;
      for (i = 0; i < 6; i++)
	force_vect[stress_p + 6 * h + i] = 0.;

      //      TODO: chemical potentials
//      force_vect[energy_p + h] += chemical_potential(ntypes, na_type[h], xi_opt + cp_start);

      // first loop over atoms: reset forces, densities
      for (i = 0; i < conf->num_atoms; i++) {
        atom = conf->atoms[i];
	val = force_vect + 3 * (conf->cnfstart + i);
        if (uf && 1 == atom->contrib) {
	  *(val++) = -atom->force.x;
	  *(val++) = -atom->force.y;
	  *val = -atom->force.z;
        } else {
	  *(val++) = 0.0;
	  *(val++) = 0.0;
	  *val = 0.0;
        }
      }
      // end first loop

      // 2nd loop: calculate pair forces and energies
      for (i = 0; i < conf->num_atoms; i++) {
        atom = conf->atoms[i];
        k = 3 * (conf->cnfstart + i);
        // loop over neighbors
        for (j = 0; j < atom->num_neighbors; j++) {
          neigh = atom->neighs[j];
	  l = neigh->col[0];
	  pot = potential->pots[l];

          // pair potential part
          if (neigh->r < pot->end) {
            // fn value and grad are calculated in the same step
            if (uf)
              phi_val = pot->splines.splint_comb_dir(pot->table, pot->d2tab, neigh->slot[0],
			      neigh->shift[0], neigh->step[0], &phi_grad);
            else
              phi_val = pot->splines.splint_dir(pot->table, pot->d2tab, neigh->slot[0],
			      neigh->shift[0], neigh->step[0]);
            // avoid double counting if atom is interacting with a copy of itself
	    if (neigh->self) {
              phi_val *= 0.5;
              phi_grad *= 0.5;
            }
            // not double force: cohesive energy
            force_vect[energy_p + h] += phi_val;

            if (uf) {
              tmp_force.x = neigh->dist.x * phi_grad;
              tmp_force.y = neigh->dist.y * phi_grad;
              tmp_force.z = neigh->dist.z * phi_grad;
	      if (1 == atom->contrib) {
		val = force_vect + k;
	        *(val++) += tmp_force.x;
  	        *(val++) += tmp_force.y;
	        *(val) += tmp_force.z;
                // actio = reactio
	        val = force_vect + 3 * (conf->cnfstart + neigh->nr);
		*(val++) -= tmp_force.x;
		*(val++) -= tmp_force.y;
		*(val) -= tmp_force.z;
	      }
              // also calculate pair stresses
              if (us) {
		val = force_vect + stress_p + 6 * h;
		*(val++) -= neigh->rdist.x * tmp_force.x;
		*(val++) -= neigh->rdist.y * tmp_force.y;
		*(val++) -= neigh->rdist.z * tmp_force.z;
		*(val++) -= neigh->rdist.x * tmp_force.y;
		*(val++) -= neigh->rdist.y * tmp_force.z;
		*(val) -= neigh->rdist.z * tmp_force.x;
              }
            }
          }
        }			// loop over neighbours

        // then we can calculate contribution of forces right away
        if (uf) {
		// TODO
          /* Weigh by absolute value of force */
//          forces[k] /= FORCE_EPS + atom->absforce;
//          forces[k + 1] /= FORCE_EPS + atom->absforce;
//          forces[k + 2] /= FORCE_EPS + atom->absforce;
          /* sum up forces */
          if (atom->contrib)
            tmpsum += conf->conf_weight * (
			    square(force_vect[k]) + square(force_vect[k + 1]) + square(force_vect[k + 2]));
        }			/* second loop over atoms */
      }

      /* energy contributions */
      force_vect[energy_p + h] *= conf->inv_num_atoms;
      force_vect[energy_p + h] -= conf->coh_energy;

      tmpsum += conf->conf_weight * settings->get_eweight() * square(force_vect[energy_p + h]);

      /* stress contributions */
      if (uf && us) {
        for (i = 0; i < 6; i++) {
          force_vect[stress_p + 6 * h + i] *= conf->inv_volume;
          force_vect[stress_p + 6 * h + i] -= *(conf->dstress[i]);
          tmpsum +=
            conf->conf_weight * settings->get_sweight() * square(force_vect[stress_p + 6 * h + i]);
        }
      }
      // limiting constraints per configuration
    }				// loop over configurations

//    /* dummy constraints (global) */
//#ifdef APOT
//    /* add punishment for out of bounds (mostly for powell_lsq) */
//    if (myid == 0) {
//      tmpsum += apot_punish(xi_opt, forces);
//    }
//#endif /* APOT */

    sum = tmpsum;		/* global sum = local sum  */
//#ifdef MPI
//    /* reduce global sum */
//    sum = 0.;
//    MPI_Reduce(&tmpsum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    /* gather forces, energies, stresses */
//    if (myid == 0) {		/* root node already has data in place */
//      /* forces */
//      MPI_Gatherv(MPI_IN_PLACE, myatoms, MPI_VECTOR, forces, atom_len,
//        atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
//      /* energies */
//      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_DOUBLE, forces + natoms * 3,
//        conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      /* stresses */
//      MPI_Gatherv(MPI_IN_PLACE, myconf, MPI_STENS, forces + natoms * 3 + nconf,
//        conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
//    } else {
//      /* forces */
//      MPI_Gatherv(forces + firstatom * 3, myatoms, MPI_VECTOR, forces, atom_len,
//        atom_dist, MPI_VECTOR, 0, MPI_COMM_WORLD);
//      /* energies */
//      MPI_Gatherv(forces + natoms * 3 + firstconf, myconf, MPI_DOUBLE,
//        forces + natoms * 3, conf_len, conf_dist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      /* stresses */
//      MPI_Gatherv(forces + natoms * 3 + nconf + 6 * firstconf, myconf, MPI_STENS,
//        forces + natoms * 3 + nconf, conf_len, conf_dist, MPI_STENS, 0, MPI_COMM_WORLD);
//    }
//#endif /* MPI */

    // root process exits this function now
    if (settings->get_myid() == 0) {
      // Increase function call counter
      inc_fcalls();
      if (std::isnan(sum)) {
#ifdef DEBUG
	io->write("\n--> Force is nan! <--\n\n");
#endif // DEBUG
	set_error_sum(1e10);
	return 1e10;
      } else {
	set_error_sum(sum);
	return sum;
      }
    }
  } // end while

  // once a non-root process arrives here, all is done.
  return 0.;
}
