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

#include <algorithm>
#include <cstdio>

#include "config.h"
#include "neighbor_2.h"

#include "../io.h"
#include "../interaction.h"
#include "../potential.h"
#include "../structures.h"
#include "../templates.h"
#include "../utils.h"

using namespace POTFIT_NS;

Neighbor_2::Neighbor_2(POTFIT *ptf) :
  Neighbor(ptf)
{}

Neighbor_2::~Neighbor_2() {}

void Neighbor_2::init(Config *conf, int i, int j, vector *dd) {
  int type1 = conf->atoms[i].type;
  int col_temp = 0;
  int k, l;

  r = sqrt(sprod(*dd,*dd));

  dist[0] = dd->x / r;
  dist[1] = dd->y / r;
  dist[2] = dd->z / r;

  type = conf->atoms[j].type;
  nr = j;
  self = (j == i) ? 1 : 0;

  // Coulomb
  r2 = square(r);

  // ADP
  rdist.x = dd->x;
  rdist.y = dd->y;
  rdist.z = dd->z;
  sqrdist.xx = dd->x * dd->x * r * r;
  sqrdist.yy = dd->y * dd->y * r * r;
  sqrdist.zz = dd->z * dd->z * r * r;
  sqrdist.yz = dd->y * dd->z * r * r;
  sqrdist.zx = dd->z * dd->x * r * r;
  sqrdist.xy = dd->x * dd->y * r * r;

  conf->atoms[i].num_neighbors++;

  col[0] = interaction->force->get_col(0, type1, type);
  structures->min_dist[col[0]] = std::min(structures->min_dist[col[0]], r);

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
  io->pexit(EXIT_FAILURE);

  return;
}
