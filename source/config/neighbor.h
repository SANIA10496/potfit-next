/****************************************************************
 *
 * neighbor.h:
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

#ifndef PTF_NEIGHBOR_H
#define PTF_NEIGHBOR_H

#include <vector>

#include "../pointers.h"
#include "../types.h"

namespace POTFIT_NS {

  class Atom;

  class Neighbor : protected Pointers {
  public:
    Neighbor(class POTFIT *);
    ~Neighbor();

    virtual void init(class Config *, int, int, vector *) = 0;
    virtual void init(class Config *, int, int, int, vector *, vector *) = 0;

    int   type;
    int   nr;
    double r;
    vector dist;			/* distance divided by r */
    int   *slot;
    double *shift;
    double *step;
    int   *col;		/* coloumn of interaction for this neighbor */

    // ADP
    vector rdist;			/* real distance */
    sym_tens sqrdist;		/* real squared distance */
    double u_val, u_grad;		/* value and gradient of u(r) */
    double w_val, w_grad;		/* value and gradient of w(r) */

    // Coulomb
    double r2;			/* r^2 */
    double fnval_el;		/* stores tail of electrostatic potential */
    double grad_el;		/* stores tail of first derivative of electrostatic potential */
    double ggrad_el;		/* stores tail of second derivative of electrostatic potential */

  };
}

#endif /* PTF_NEIGHBOR_H */
