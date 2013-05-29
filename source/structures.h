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

namespace POTFIT_NS {

  class Structures : protected Pointers {
  public:
    Structures(class POTFIT *);
    ~Structures();

    void init(void);

    void read_config(FILE *);

    void print_mindist(void);

    int get_num_contrib_atoms(void);
    int get_num_contrib_energies(void);
    int get_num_contrib_stresses(void);
    int get_num_total_atoms(void);
    int get_num_total_configs(void);

    int ntypes; 			// number of atom types
    std::vector<int> num_per_type; 	// total number of atoms per type
    int using_forces; 			// total number of configs using forces
    int using_stresses;			// total number of configs using stress
    int firstconf; 			// number of first configuration (for MPI)
    int nconf; 				// number of configurations (for MPI)

    double *min_dist; 			// minimal distances for all interactions

    std::vector<Config *> config; 	// configurations
  private:
    int line;
    int total_num_conf; 		// total number of configurations
    int total_num_atoms;		// total number of atoms
  };

}

#endif /* PTF_SETTINGS_H */
