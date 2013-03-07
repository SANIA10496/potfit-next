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
#include "table.h"

namespace POTFIT_NS {

  class Potential : protected Pointers {
  public:
    Potential(class POTFIT *);
    ~Potential();

    void init(int);
    void read_globals(FILE *);
    void read_potentials(FILE *);

    int enable_cp;
    int format;
    int have_grad;
    int n_invar_pots;

    int number; 		// total number of potentials
    int total_par; 		// total number of free parameters
    int total_ne_par; 		// total number of non-electrostatic parameters

    // global parameters for analytic potentials
    int have_globals;
    int number_globals;
    int *globals_usage;
    int ***globals_idx;

    // coulomb parameters
    double *ratio;
    double *charge;
    double last_charge;
    double *dp_kappa;
    int sw_kappa;

    // dipole parameters
    double *dp_alpha;
    double *dp_b;
    double *dp_c;

    // gradient and invar_pot data from potential header
    int *gradient;
    int *invar_pot;

    // actual potential data
    Table **pots;
    Table *global_params; 	// bare table for global parameters (apot)
    Table *chem_pot; 		// bare table for chemical potentials (pair)
    Table *calc_pot; 		// table for calculations
  private:

  };

}

#endif /* PTF_POTENTIAL_H */
