/****************************************************************
 *
 * potential.h:
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

#ifndef PTF_POTENTIAL_H
#define PTF_POTENTIAL_H

#include <cstring>
#include <vector>

#include "pointers.h"

#include "tables/table.h"
#include "tables/opt_table.h"
#include "tables/globals_table.h"
#include "tables/chempot_table.h"

namespace POTFIT_NS {

  class Potential : protected Pointers {
  public:
    Potential(class POTFIT *);
    ~Potential();

    void init(int);
    void read_globals(FILE *);
    void read_potentials(FILE *);

    void update_potentials(void);

    void init_opt_table(void);

    std::vector<char *> elements;

    int enable_cp; 		// chemical potential (only for pair)
    int num_globals; 		// global potentials

    int num_pots; 		// total number of potentials
    int num_free_pots; 		// total number of invariant potentials

    int num_params; 		// total number of free parameters
    int num_free_params; 	// total number of invariant parameters

    // invar_pot data from potential header
    int *invar_pot;

    int *idxpot; 		// indirect index for potentials
    int *idxparam; 		// indirect index for parameters

    double *rcut;
    double *rmin;
    double rcut_max;

    // actual potential data
    Table 		**pots; 		// individual potentials
    OptTable 		*opt; 			// optimization table
    GlobalsTable 	*global_params; 	// table for global parameters (apot)
    ChempotTable 	*chem_pot; 		// table for chemical potentials (pair)
  };

}

#endif /* PTF_POTENTIAL_H */
