/****************************************************************
 *
 * table_analytic.h:
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

#ifndef PTF_TABLE_ANALYTIC_H
#define PTF_TABLE_ANALYTIC_H

#include <iostream>

#include "table.h"

#include "../types.h"
#include "../functions/function.h"

namespace POTFIT_NS {

class TableAnalytic : public Table {
  public:
    TableAnalytic(class POTFIT *);
    ~TableAnalytic();

    // initialize with name and invar
    void init(const char *, int);
    void init_bare(const char *, int);

    void read_potential(FILE *);

    void set_cutoff(double);
    void set_value(int, double, double);
    void set_value(int, const char *, double, double, double);

    void update_value(int, double);
    double calc_value(double);

    const char *get_param_name(int);

    void write_potential(FILE *);
    void write_plot(FILE *);
    void write_plotpoint(FILE *);

  private:
    int bare;
    int smooth_pot;

    double *values; 	// parameter values
    char **param_name; 	// parameter names
    double *val_min; 	// parameter minimum
    double *val_max; 	// parameter maximum
    int *invar_par; 	// invariant parameters

    Function *function; 	// function pointer for analytic potentials
};
}

#endif /* PTF_TABLE_ANALYTIC_H */
