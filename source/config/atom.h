/****************************************************************
 *
 * atom.h:
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

#ifndef PTF_ATOM_H
#define PTF_ATOM_H

#include <vector>

#include "neighbor.h"

#include "../pointers.h"
#include "../types.h"

namespace POTFIT_NS {

  class Atom : protected Pointers {
  public:
    Atom(class POTFIT *);
    ~Atom();

    int type;
    int num_neighbors;
    vector pos;
    vector force;
    double absforce;
    int contrib;

    std::vector<Neighbor *> neighs;

    // EAM
    double rho;
    double gradf;

    // ADP
    vector mu;
    sym_tens lambda;
    double nu;

    // dipole
    vector E_stat;
    vector p_sr;
    vector E_ind;
    vector p_ind;
    vector E_old;
    vector E_tot;
  };

}

#endif /* PTF_CONFIG_H */
