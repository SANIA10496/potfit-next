/****************************************************************
 *
 * table.cpp:
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

#include <cstdlib>

#include "table.h"

using namespace POTFIT_NS;

Table::Table(POTFIT *ptf) :
  Pointers(ptf),
  splines(ptf),
  init_done(0),
  format(0),
  begin(0.0),
  end(0.0),
  plotmin(0.0),
  pot_number(0),
  num_params(0),
  num_free_params(0),
  opt_pot_start(-1),
  len(0),
  step(0.0),
  invstep(0.0),
  xcoord(NULL),
  table(NULL),
  d2tab(NULL)
{}

Table::~Table() {}
