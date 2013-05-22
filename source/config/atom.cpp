/****************************************************************
 *
 * config.cpp:
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

#include "atom.h"

using namespace POTFIT_NS;

Atom::Atom(POTFIT *ptf) : Pointers(ptf) {
  type = 0;
  num_neighbors = 0;
  pos.x = 0.0;
  pos.y = 0.0;
  pos.z = 0.0;
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  absforce = 0.0;
  contrib = 0;

  // EAM
  rho = 0.0;
  gradf = 0.0;

  // ADP
  mu.x = 0.0;
  mu.y = 0.0;
  mu.z = 0.0;
  lambda.xx = 0.0;
  lambda.yy = 0.0;
  lambda.zz = 0.0;
  lambda.xx = 0.0;
  lambda.xx = 0.0;
  lambda.xx = 0.0;
  nu = 0.0;

  // dipole
  E_stat.x = 0.0;
  E_stat.y = 0.0;
  E_stat.z = 0.0;
  p_sr.x = 0.0;
  p_sr.y = 0.0;
  p_sr.z = 0.0;
  E_ind.x = 0.0;
  E_ind.y = 0.0;
  E_ind.z = 0.0;
  p_ind.x = 0.0;
  p_ind.y = 0.0;
  p_ind.z = 0.0;
  E_old.x = 0.0;
  E_old.y = 0.0;
  E_old.z = 0.0;
  E_tot.x = 0.0;
  E_tot.y = 0.0;
  E_tot.z = 0.0;

  return;
}

Atom::~Atom() {
  return;
}
