/****************************************************************
 *
 * table.h:
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

#ifndef PTF_TABLES_H
#define PTF_TABLES_H

#include <iostream>

#include "../pointers.h"

namespace POTFIT_NS {

  class Table : protected Pointers {
  public:
    Table(class POTFIT *);
    ~Table();

    // initialize with potential name
    virtual void init(const char *, int) = 0;
    // use raw potential for additional parameters
    virtual void init_bare(const char *, int) = 0;

    // read one potential
    virtual void read_potential(FILE *) = 0;

    // for tabulated potentials
    virtual void set_value(int, double, double) = 0;
    // for analytic potentials
    virtual void set_value(int, const char*, double, double, double) = 0;
    virtual const char *get_param_name(int) = 0;
  protected:
    char *name; 	// name of analytic function / potential type

    double begin; 	// starting position of potential = r_min
    double end; 	// end position of potential = cutoff radius

    int pot_number; 	// index of potential
    int n_par; 		// number of parameters
    int n_invar; 	// number of invariant parameters
    int init_done; 	// is the table initialized?
  };
}

#endif // PTF_TABLES_H
