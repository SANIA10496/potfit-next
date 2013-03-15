/****************************************************************
 *
 * chempot_table.h:
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

#ifndef PTF_CHEMPOT_TABLE_H
#define PTF_CHEMPOT_TABLE_H

#include "../pointers.h"

namespace POTFIT_NS {

  class ChempotTable : protected Pointers {
  public:
    ChempotTable(class POTFIT *, int);
    ~ChempotTable();

    void add_value(int, const char*, double, double, double);
    void set_value(int, double);

  private:
    int num_params;		// number of chemical potentials
    int num_invar_params; 	// number of invariant parameters

    double *values; 		// values
    char **param_name; 		// parameter names
    double *val_min; 		// parameter minimum
    double *val_max; 		// parameter maximum
    int *invar_par; 		// invariant parameters
  };
}

#endif // PTF_GLOBALS_TABLE_H
