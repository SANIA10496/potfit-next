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

#include <cstdio>

#include "pointers.h"

#include "tables/table.h"
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

    // no longer used, will be removed
    int format;

    // flags for
    int enable_cp; 		// chemical potential (only for pair)
    int enable_globals; 	// global potentials
    int have_grad; 		// ??
    int n_invar_pots; 		// number of invariant potentials

    int number; 		// total number of potentials
    int total_par; 		// total number of free parameters
    int total_ne_par; 		// total number of non-electrostatic parameters (??)

    // gradient and invar_pot data from potential header
    int *gradient;
    int *invar_pot;

    // actual potential data
    Table 		**pots; 		// individual potentials
    GlobalsTable 	*global_params; 	// table for global parameters (apot)
    ChempotTable 	*chem_pot; 		// table for chemical potentials (pair)
    Table 		*calc_pot; 		// table for calculations
  private:

  };

}

#endif /* PTF_POTENTIAL_H */
