/****************************************************************
 *
 * error.h:
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

#ifndef PTF_ERROR_H
#define PTF_ERROR_H

#include <cstdio>

#include "pointers.h"

namespace POTFIT_NS {

  class Error : protected Pointers {
  public:
    Error(class POTFIT *);
    ~Error();

    void write_report(void);
  private:
    void calc_errors(void);

    double *force_vect;
    double total_sum;
    double force_sum, energy_sum, stress_sum, punish_sum;
    double rms_force, rms_energy, rms_stress;

    int total_contrib;
    int num_forces, num_energies, num_stresses;
    int fcalls;
  };
}

#endif /* PTF_ERROR_H */
