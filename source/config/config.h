/****************************************************************
 *
 * config.h:
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

#ifndef PTF_CONFIG_H
#define PTF_CONFIG_H

#include <iostream>
#include <vector>

#include "atom.h"

#include "../pointers.h"
#include "../types.h"

namespace POTFIT_NS {

  class Config : protected Pointers {
  public:
    Config(class POTFIT *);
    ~Config();

    void read(FILE *, int *);

    std::vector<Atom *> atoms;
    std::vector<int> num_per_type;
    int use_forces;
    int use_stresses;
    int num_atoms;

    double coh_energy;
    double conf_weight;
    double volume;

    sym_tens stress;
    vector box_x, box_y, box_z;

  private:
    void calc_neighbors(void);
    vector vec_prod(vector, vector);
  };

}

#endif /* PTF_CONFIG_H */
