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
  type(-1),
  nr(-1),
  r(0.0),
  u_val(0.0),
  u_grad(0.0),
  w_val(0.0),
  w_grad(0.0),
  // Coulomb
  r2(0.0),
  fnval_el(0.0),
  grad_el(0.0),
  ggrad_el(0.0)
{}

Neighbor::~Neighbor() {}
