/****************************************************************
 *
 * neighbor.cpp:
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

#include "neighbor.h"

using namespace POTFIT_NS;

Neighbor::Neighbor(POTFIT *ptf) :
  Pointers(ptf),
  // ADP
  u_val(0.0),
  u_grad(0.0),
  w_val(0.0),
  w_grad(0.0),
  // Coulomb
  r2(0.0),
  fnval_el(0.0),
  grad_el(0.0),
  ggrad_el(0.0)
{
  type= -1;
  nr = -1;
  self = 0;
  r = 0.0;
  dist[0] = 0.0;
  dist[1] = 0.0;
  dist[2] = 0.0;

  for (int i=0;i<MAX_SLOTS;i++) {
    slot[i] = 0;
    shift[i] = 0.0;
    step[i] = 0.0;
    col[i] = 0;
  }

  // ADP
  rdist.x = 0.0;
  rdist.y = 0.0;
  rdist.z = 0.0;
  sqrdist.xx = 0.0;
  sqrdist.yy = 0.0;
  sqrdist.zz = 0.0;
  sqrdist.zx = 0.0;
  sqrdist.xy = 0.0;
  sqrdist.yz = 0.0;
}

Neighbor::~Neighbor() {}
