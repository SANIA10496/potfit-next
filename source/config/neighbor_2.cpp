/****************************************************************
 *
 * neighbor_2.cpp:
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

#include <cstdio>
#include <cmath>

#include "config.h"
#include "neighbor_2.h"

#include "../io.h"
#include "../interaction.h"
#include "../potential.h"
#include "../structures.h"
#include "../templates.h"
#include "../utils.h"

using namespace POTFIT_NS;

Neighbor_2::Neighbor_2(POTFIT *ptf) : Neighbor(ptf) {
  type = -1;
  nr = -1;
  r = 0.0;
  dist.x = 0.0;
  dist.y = 0.0;
  dist.z = 0.0;

  slot = new int[interaction->force->num_slots()];
  shift = new double[interaction->force->num_slots()];
  step = new double[interaction->force->num_slots()];
  col = new int[interaction->force->num_slots()];

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

  u_val = 0.0;
  u_grad = 0.0;
  w_val = 0.0;
  w_grad = 0.0;

  // Coulomb
  r2 = 0.0;
  fnval_el = 0.0;
  grad_el = 0.0;
  ggrad_el = 0.0;

  return;
}

Neighbor_2::~Neighbor_2() {
  delete [] slot;
  delete [] shift;
  delete [] step;
  delete [] col;

  return;
}

void Neighbor_2::init(Config *conf, int i, int j, vector *dd) {
  int type1 = conf->atoms[i]->type;
  int col_temp = 0;
  int k, l;

  r = sqrt(SPROD(*dd,*dd));

  dist.x = dd->x / r;
  dist.y = dd->y / r;
  dist.z = dd->z / r;

  type = conf->atoms[j]->type;
  nr = j;

  // Coulomb
  r2 = r * r;

  // ADP
  rdist.x = dd->x * r;
  rdist.y = dd->y * r;
  rdist.z = dd->z * r;
  sqrdist.xx = dd->x * dd->x * r * r;
  sqrdist.yy = dd->y * dd->y * r * r;
  sqrdist.zz = dd->z * dd->z * r * r;
  sqrdist.yz = dd->y * dd->z * r * r;
  sqrdist.zx = dd->z * dd->x * r * r;
  sqrdist.xy = dd->x * dd->y * r * r;

  conf->atoms[i]->num_neighbors++;

  col[0] = interaction->force->get_col(0, type1, type);
  structures->min_dist[col[0]] = MIN(structures->min_dist[col[0]], r);

  // loop over all slots
  for (k=0;k<interaction->force->num_slots();k++) {
    col_temp = interaction->force->get_col(k, type1, type);
    l = 0;
    while (potential->pots[col_temp]->xcoord[l] < r) {
      l++;
    }
    slot[k] = --l;
    step[k] = potential->pots[col_temp]->xcoord[l+1] - potential->pots[col_temp]->xcoord[l];
    shift[k] = (r - potential->pots[col_temp]->xcoord[l]) / step[k];
    col[k] = col_temp;
  }


  return;
}

void Neighbor_2::init(Config *conf, int i, int j, int k,  vector *dd_ij, vector *dd_ik) {
  io->error << "The three-body neighbor function cannot be called for two-body neighbor lists!" << std::endl;

  return;
}
