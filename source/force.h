/****************************************************************
 *
 * force.h:
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

#ifndef PTF_FORCE_H
#define PTF_FORCE_H

#include <iostream>

#include "pointers.h"

namespace POTFIT_NS {

  class Force : protected Pointers {
  public:
    Force(class POTFIT *);
    ~Force();

    virtual void read_additional_data(FILE *) = 0;

    virtual int num_slots(void) = 0;
    virtual int neigh_type(void) = 0;
    virtual int get_col(int, int, int) = 0;
    virtual void update_min_dist(double *) = 0;

    virtual double calc_forces(void) = 0;

    virtual int cols(void) = 0;

    void calc_pointers(void);

    int get_fcalls(void);
    void inc_fcalls(void);

    double *force_vect; 	// all deviations
    int energy_p; 		// pointer for energies
    int stress_p; 		// pointer for stresses
    int dummy_p; 		// pointer to dummy constraints
    int limit_p; 		// pointer to limiting constraints
    int punish_par_p; 		// pointer to parameter punishments
    int punish_pot_p; 		// pointer to potential punishments

  private:
    int fcalls;
  };

}

#endif // PTF_FORCE_H
