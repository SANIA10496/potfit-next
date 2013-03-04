/****************************************************************
 *
 * settings.cpp:
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

#include <mpi.h>

#include "settings.h"

using namespace POTFIT_NS;

Settings::Settings(POTFIT *ptf) : Pointers(ptf) {
  myid = -1;
  num_cpus = 0;
  opt = 1;

  eweight = -1.;
  sweight = -1.;
  extend = 0;
  evo_threshold = 0;
  anneal_temp = 0;
  apot_punish_value = 0;
  d_eps = 0;
  dp_cut = 0;
  dp_tol = 0;
  dp_mix = 0;
}

Settings::~Settings() {}
