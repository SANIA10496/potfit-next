/****************************************************************
 *
 * structures.h:
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

#ifndef PTF_STRUCTURES_H
#define PTF_STRUCTURES_H

#include <vector>

#include "pointers.h"

#include "config/config.h"
#include "config/list_neighbors.h"

namespace POTFIT_NS {

  class Structures : protected Pointers {
  public:
    Structures(class POTFIT *);
    ~Structures();

    void init(void);

    void read_config(FILE *);
    void update_pointers(const int &);

    void print_mindist(void);

    const int get_num_contrib_atoms(void);
    const int get_num_contrib_energies(void);
    const int get_num_contrib_stresses(void);
    const int get_num_total_atoms(void);
    const int get_num_total_configs(void);

    void set_ntypes(const int &);
    int get_ntypes(void);

    std::vector<int> num_per_type; 	// total number of atoms per type
    int using_forces; 			// total number of configs using forces
    int using_stresses;			// total number of configs using stress
    // MPI communication related
    int firstconf; 			// number of first configuration
    int nconf; 				// number of configurations
    int myatoms; 			// number of atoms
    int firstatom; 			// number of first atom on this cpu
    int *atom_len;
    int *atom_dist;
    int *conf_len;
    int *conf_dist;

    double *min_dist; 			// minimal distances for all interactions

    std::vector<Config *> config; 	// configurations
    std::vector<Atom> atoms; 		// atoms
    std::vector<Neighbor_2> neigh_2; 	// two-body neighbors
    std::vector<Neighbor_3> neigh_3; 	// three-body neighbors

  private:
    int ntypes; 			// number of atom types
    int line;
    int total_num_conf; 		// total number of configurations
    int total_num_atoms;		// total number of atoms
  };

}

#endif /* PTF_SETTINGS_H */
