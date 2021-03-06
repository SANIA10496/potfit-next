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

    // initialize with name and index
    void init(const std::string &, const int &);

    void read_potential(FILE *);

    int get_number_params(void);
    int get_number_free_params(void);
    double get_cutoff(void);
    double get_rmin(void);
    double get_val_min(const int &);
    double get_val_max(const int &);
    double get_plotmin(void);
    double get_value(const double &);

    void init_calc_table(void);
    void set_param(const int &, const double &);
    void update_potential(const int &);
    void update_slots(void);

    void write_potential(std::ofstream &);
    void write_plot(FILE *);
    void write_plotpoint(FILE *);

  private:
    void update_values(void);
    void update_calc_table(const int &);

    int smooth_pot;

    double cutoff(const double &, const double &, const double &);

    std::vector<std::string> param_name; 	// parameter names
    std::vector<double> val_min; 		// parameter minimum
    std::vector<double> val_max; 		// parameter maximum
    std::vector<double> stored_values; 		// values currently used in the calc_pot tables
    double fval; 				// temp variable

    Function *function; 	// function pointer for analytic potentials
  };
}

#endif /* PTF_TABLE_ANALYTIC_H */
