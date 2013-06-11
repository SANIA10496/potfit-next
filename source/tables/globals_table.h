/****************************************************************
 *
 * globals_table.h:
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

#ifndef PTF_GLOBALS_TABLE_H
#define PTF_GLOBALS_TABLE_H

#include <iostream>

#include "../pointers.h"

namespace POTFIT_NS {

  class GlobalsTable : protected Pointers {
  public:
    GlobalsTable(class POTFIT *, int);
    ~GlobalsTable();

    void add_param(int, const char *, double, double, double);

    int get_number_params(void);
    int get_number_free_params(void);

    void check_usage(void);

    void set_value(int, double);
    void set_opt_pot_start(const int &);

    int get_index(const char *);
    void add_usage(int, int, int);
    void get_value(int, double *);

    void update_potentials(void);

    void get_values(int *, double *);

    void write_potential(std::ofstream &);

  private:
    int num_globals; 		// number of global parameters
    int num_free_globals; 	// number of free global parameters
    int opt_pot_start;

    char **param_name; 		// parameter names
    double *values; 		// parameter values
    double *val_min; 		// parameter minimum
    double *val_max; 		// parameter maximum
    int *invar_par; 		// invariant parameters
    int *idx; 			// indirect index
    int *usage; 		// how often each parameter is used
    int ***global_idx; 		// where each parameter is used
  };
}

#endif // PTF_GLOBALS_TABLE_H
